/*
计算两个非简谐一维势能的时间关联函数
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>  // 计算程序运行时间
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

// using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

const int N = 1;            // 系统的自由度
const int P = 200;          // beads数量
const int steps = 20000;  // 演化的总步数
const int step_leap = 5;  // 跳过step_leap步记录一次位置和动量 计算时间关联函数的最小步长
const double m = 1;     // 质量矩阵的对角元
const double beta = 8;  // 温度的倒数
const double dt_PIMD = 0.01 * beta / sqrt(P); // PIMD中的时间步长
const double dt_HD = 0.001;
const double omega_P = sqrt(P) / beta;
const int n_x = 1000;  // 坐标空间中通过PIMD采样的点数目
const int n_p = 10;  // 固定x在动量空间中采样的点数目
const int PIMD_steps = 800000; // 坐标空间抽样时使用的PIMD步数

double M[N];      // 质量矩阵 1-D array
double M_inv[N];  // 质量矩阵的逆 1-D array

// 尽量不要动下面定义的变量
int t_steps_max = steps / step_leap - 200;
double x_array[P][N];     // 储存所有beads位置
double p_array[P][N];     // 储存所有beads动量
double nabla_V[N] = { 0 };  // 储存单个beads在Cartesian坐标下的势能梯度
double nabla_phi[P][N] = { 0 };  // 每个bead在staging坐标下的受力
double traj_x_init[n_p * n_x][N] = { 0 };  // 储存选定为初始值的平衡构型
double traj_p_init[n_p * n_x][N] = { 0 };  // 用作动力学的初始值
double x_list_traj[n_p * n_x][steps / step_leap][N];  // 所有轨线不同时刻的位置
double p_list_traj[n_p * n_x][steps / step_leap][N];  // 所有轨线不同时刻的动量
double t_list[steps / step_leap];  // 记录每个构型对应的演化时间

// 设置随机数种子
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> distribution(0.0, 1.0);

/*
staging变换，将Cartesian坐标变换成为staging坐标

x: double 2-D array 所有beads的Cartesian坐标
没有返回值，直接修改传入的数组
*/
void staging_transf(double(&x)[P][N]) {
  // 第一个bead的坐标不用动
  for (int i = 1; i < P; i++) {
    for (int j = 0; j < N; j++) {
      x[i][j] = x[i][j] - (double)(i * x[(i + 1) % P][j] + x[0][j]) / (i + 1);
    }
  };
}

/*
staging逆变换，将staging坐标转换成为Cartesian坐标
s: double 2-D array 所有beads的staging坐标
没有返回值
*/
void staging_transf_inv(double(&s)[P][N]) {
  // 第一个bead不用动
  for (int j = P - 1; j > 0; j--) {
    for (int i = 0; i < N; i++) {
      s[j][i] = s[j][i] + (double)(j * s[(j + 1) % P][i] + s[0][i]) / (j + 1);
    }
  }
}

/*
Cartesian坐标下真实的四次方势
x: doubel 1-D array Cartesian坐标
return: double 对应坐标下的势能
*/
double quartic_potential(double(&x)[N]) {
  return 0.25 * x[0] * x[0] * x[0] * x[0];
}
/*
不对称非简谐势能
x: double 1-D array Cartesian坐标
retrun: double 对应坐标下的势能
*/
double asy_aharmonic_potential(double(&x)[N]) {
  double x_2 = x[0] * x[0];
  return x_2 - 0.1 * x_2 * x[0] + 0.1 * x_2 * x_2;
}

/*
Cartesian坐标下四次方势的导数
x: double 1-D array Cartesian坐标
return: double 1-D array 势能的梯度
*/
void nabla_quartic_potential(double(&x)[N]) {
  nabla_V[0] = x[0] * x[0] * x[0];
}

/*
Cartesian坐标下非对称非简谐势能的导数
x: double 1-D array Cartesian坐标
return: double 1-D array 势能的梯度
*/
void nabla_asy_aharmonic_potential(double(&x)[N]) {
  nabla_V[0] = 2.0 * x[0] - 0.3 * x[0] * x[0] + 0.4 * x[0] * x[0] * x[0];
}

/*
四次方势能的质量加权Hessian矩阵（实际上就是一个数）
x: double 1-D array Cartesian坐标
return: double 1-D array 势能的梯度
*/
double Mhessian_quartic(double(&x)[N]) {
  return 3.0 * x[0] * x[0] / m;
}

/*
非简谐非对称势能的质量加权Hessian矩阵（实际上就是一个数）
x: double 1-D array Cartesian坐标
return: double 1-D array 势能的梯度
*/
double Mhessian_asy_aharmonic(double(&x)[N]) {
  return (2.0 - 0.6 * x[0] + 1.2 * x[0] * x[0]) / m;
}

/*
用于计算热质量矩阵
x: double 1-D array N 必须是cartesian坐标中的位置
*/
double M_therm(double(&x)[N], double (*Mhessian)(double(&)[N])) {
  double u = beta / 2.0;
  double w_2 = Mhessian(x);
  if (w_2 >= 0) {
    u = u * sqrt(w_2);
    return u / tanh(u);
  }
  else {
    u = u * sqrt(-w_2);
    return tanh(u) / u;
  }
}

/*
staging坐标下的所有beads的外势
s_chain: double 2-D array 所有beads的staging坐标
return: double 势能
*/
double phi_potential(double(&s_chain)[P][N],
  double (*potential)(double(&)[N])) {
  staging_transf_inv(s_chain);
  double phi = 0;
  for (int j = 0; j < P; j++) {
    phi += potential(s_chain[j]);
  }
  staging_transf(s_chain);  // 最后需要将Cartesian坐标变换为staging坐标
  return (double)phi / P;
}

/*
staging坐标下所有beads的势能梯度，在求解Langevin方程中直接使用
s_chain: double 2-D array 所有beads的staging坐标
return: douible 2-D array 所有beads在staging坐标下的势能梯度
注意这里没有修改传入的数组
*/
void nabla_phi_potential(double(&s_chain)[P][N],
  void (*nabla_potential)(double(&)[N])) {
  // 首先将数组初始化
  for (int j = 0; j < P; j++) {
    for (int k = 0; k < N; k++) {
      nabla_phi[j][k] = 0.0;
    }
  }
  staging_transf_inv(s_chain);
  // 将staging坐标变换为Cartesian坐标，方便计算势能
  {
    for (int j = 0; j < P; j++) {
      nabla_potential(s_chain[j]);  // 这一步修改了nabla_V中的值
      // 经过staging_transf_inv之后，s_chain表示对应的Cartesian坐标
      nabla_phi[0][0] += nabla_V[0];
    }
    nabla_phi[0][0] = nabla_phi[0][0] / P;
  }
  // 后面的每个bead的势能梯度都需要递推计算
  for (int j = 1; j < P; j++) {
    nabla_potential(s_chain[j]);  // 计算第j个cartisan坐标下的梯度
    nabla_phi[j][0] =
      (double)nabla_V[0] / P + (double)(j - 1) / j * nabla_phi[j - 1][0];
  }
  // 在程序结束前再将Cartesian坐标变换为staging坐标
  staging_transf(s_chain);
}

/*
动能的viral estimator, 直接输入staging坐标
x_chain: double 2-D array 所有beads的staging坐标
return: double 当前构型下的动能
*/
double kinetic_energy_viral(double(&s_chain)[P][N],
  double beta,
  void (*nabla_potential)(double(&)[N])) {
  staging_transf_inv(s_chain);
  double x_c[N] = { 0 };  // 构造所有beads质心的坐标
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < P; j++) {
      x_c[k] += s_chain[j][k];
    }
    x_c[k] = (double)x_c[k] / P;
  }
  double T = 0;
  for (int j = 0; j < P; j++) {
    nabla_potential(s_chain[j]);
    for (int k = 0; k < N; k++) {
      T += (s_chain[j][k] - x_c[k]) * nabla_V[k];
    }
  }
  staging_transf(s_chain);
  return (double)T / (2 * P) + (double)N / (2 * beta);
}

// 定义一些BAOAB函数常用的常数
const double c1_PIMD = exp(-omega_P * dt_PIMD);
const double c2_PIMD = sqrt((1 - c1_PIMD * c1_PIMD) / beta);
const double d1_PIMD = cos(omega_P * dt_PIMD / 2);  // 演化内势力
const double d2_PIMD = sin(omega_P * dt_PIMD / 2);
/*
以BAOAB的方式演化Langevin方程，均为staging坐标
x_chain: double 2-D array 表示前一个时刻所有beads的位置
p_chain: double 2-D array 表示前一个时刻所有beads的动量
没有返回值 直接修传入的数组
*/
void BAOAB(double(&s_chain)[P][N],
  double(&p_chain)[P][N],
  void (*nabla_potential)(double(&)[N]),
  double beta,
  double omega_P,
  double c1_PIMD,
  double c2_PIMD,
  double d1_PIMD,
  double d2_PIMD) {
  // 首先演化半步B动量项
  nabla_phi_potential(s_chain, nabla_potential);
  for (int j = 0; j < P; j++) {
    p_chain[j][0] += -0.5 * dt_PIMD * nabla_phi[j][0];
  }

  // 其次演化<半步>A内势力项
  s_chain[0][0] += +0.5 * dt_PIMD / M[0] * p_chain[0][0];
  double var = 0;  // 由于储存计算过程中的中间变量
  for (int j = 1; j < P; j++) {
    var = s_chain[j][0] * d1_PIMD +
      d2_PIMD * j / ((j + 1) * M[0] * omega_P) * p_chain[j][0];
    p_chain[j][0] = -omega_P * d2_PIMD * (j + 1) / j * M[0] * s_chain[j][0] +
      d1_PIMD * p_chain[j][0];
    s_chain[j][0] = var;
  }

  // 然后演化一步O控温项
  p_chain[0][0] =
    c1_PIMD * p_chain[0][0] + c2_PIMD * sqrt(M[0]) * distribution(generator);
  for (int j = 1; j < P; j++) {
    p_chain[j][0] =
      c1_PIMD * p_chain[j][0] +
      c2_PIMD * sqrt((double)(j + 1) / j * M[0]) * distribution(generator);
  }

  // 再演化半步A内势力项
  s_chain[0][0] += +0.5 * dt_PIMD / M[0] * p_chain[0][0];
  for (int j = 1; j < P; j++) {
    var = s_chain[j][0] * d1_PIMD +
      d2_PIMD * j / ((j + 1) * M[0] * omega_P) * p_chain[j][0];
    p_chain[j][0] = -omega_P * d2_PIMD * (j + 1) / j * M[0] * s_chain[j][0] +
      d1_PIMD * p_chain[j][0];
    s_chain[j][0] = var;
  }

  // 最后演化半步B动量项
  nabla_phi_potential(s_chain, nabla_potential);
  for (int j = 0; j < P; j++) {
    p_chain[j][0] += -0.5 * dt_PIMD * nabla_phi[j][0];
  }
}

/*
计算给定构型的温度
p_chain: double 2-D array 所有beads的动量
return: double 温度
*/
double temperature(double(&p_chain)[P][N]) {
  double T = 0;
  for (int k = 0; k < N; k++) {
    T += p_chain[0][k] * M_inv[k] * p_chain[0][k];
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      T += p_chain[j][k] * j * M_inv[k] * p_chain[j][k] / (j + 1);
    }
  }
  return T / (N * P);
}

/*
kubo形式的动量关联函数，只计算了两点的关联函数，没有进行时间平均
start: int
end: int
*/
double p_kubo_corr_func_(int start, int end, double (*Mhessian)(double(&)[N])) {
  double ave = 0;
  for (int j = 0; j < n_x * n_p; j++) {  // 对所有轨线平均
  // 对于高维系统应添加对维数的循环 for (int k = 0; k < n_atom * n_dim; k++){}
    ave += M[0] / M_therm(x_list_traj[j][start], Mhessian) * p_list_traj[j][start][0] *
      p_list_traj[j][end][0];
  }

  return ave / (n_x * n_p);
}
/*
kubo形式动量关联函数，考虑了在同一轨线不同时间的时间平均
*/
double p_kubo_corr_func(int t, double (*Mhessian)(double(&)[N]), int cnt = 50) {
  double ave = 0;
  for (int j = 0; j < cnt; j++) {
    ave += p_kubo_corr_func_(j, j + t, Mhessian);
  }
  return ave / cnt;
}

/*
标准形式位置pingfang关联函数的实部
*/
double x_square_std_corr_func_(int start, int end, double (*Mhessian)(double(&)[N])) {
  double ave = 0;
  for (int j = 0; j < n_x * n_p; j++) {
    double M_t = M_therm(x_list_traj[j][start], Mhessian);
    ave += (x_list_traj[j][start][0] * x_list_traj[j][start][0] +
      beta / (4.0 * M_t) - beta * beta * p_list_traj[j][start][0] * p_list_traj[j][start][0] /
      (4.0 * M_t * M_t)) * x_list_traj[j][end][0] * x_list_traj[j][end][0];
  }

  return ave / (n_x * n_p);
}

/*
时间平均的位置平方关联函数
*/
double x_square_std_corr_func(int t, double (*Mhessian)(double(&)[N]), int cnt = 50) {
  double ave = 0;
  for (int j = 0; j < cnt; j++) {
    ave += x_square_std_corr_func_(j, j + t, Mhessian);
  }
  return ave / cnt;
}

/*
哈密顿动力学
*/
void Velocity_Verlet(double(&x)[N],
  double(&p)[N],
  void (*nabla_V_potential)(double(&)[N]),
  double dt)
{
  nabla_V_potential(x);
  for (int i = 0; i < N; i++) // 半步动量
  {
    p[i] = p[i] - dt / 2.0 * nabla_V[i];
  }
  for (int i = 0; i < N; i++) // 一步位置
  {
    x[i] = x[i] + p[i] / M[i] * dt;
  }
  nabla_V_potential(x);
  for (int i = 0; i < N; i++)
  {
    p[i] = p[i] - dt / 2.0 * nabla_V[i];
  }
}
/*
最终的计算程序，输入势能梯度函数与质量加权hessian矩阵的函数
*/
void calculate(void (*nabla_V_potential)(double(&)[N]), double (*Mhessian)(double(&)[N])) {
  M[0] = m;
  M_inv[0] = 1 / m;

  time_t t_start, t_end;  // 记录程序运行时间
  t_start = time(NULL);

  // 构造PIMD的初始条件
  for (int j = 0; j < P; j++) {
    x_array[j][0] = 0;
    p_array[j][0] = 0;
  }
  staging_transf(x_array); // 将Cartesian坐标转化到staging坐标
  double T = 0; // 用于判断PIMD是否达到平衡?

  cout << "------------------START--PIMD------------------" << endl;
  cout << "the number of beads is " << P << endl;
  cout << "the setting temperature is " << beta << endl;
  cout << "total steps is " << PIMD_steps << endl;

  for (int i = 0; i < PIMD_steps; i++) {  // 用PIMD对量子正则分布进行抽样
    BAOAB(x_array, p_array, nabla_V_potential, beta, omega_P, c1_PIMD,
      c2_PIMD, d1_PIMD, d2_PIMD);
    T += temperature(p_array) / PIMD_steps;
    if (i % 50000 == 0) {
      t_end = time(NULL);
      printf("in PIMD, step = %d, time has been used = %.3f s\n", i,
        difftime(t_end, t_start));
    }
    // 
    if (i >= PIMD_steps - n_x) {  // 首先在位形空间采样，作为初始值
      for (int k = 0; k < N; k++) {
        for (int np = 0; np < n_p; np++) { // 每n_p个轨线采用相同的初始位置
          traj_x_init[(i - PIMD_steps + n_x) * n_p + np][k] = x_array[0][k];
        }
      }
      // 还需要对动量进行采样
      for (int k = 0; k < N; k++) {  // 空间维数
        for (int np = 0; np < n_p; np++) {  // 对于每一个位置，动量空间的采样
          traj_p_init[np + (i - PIMD_steps + n_x) * n_p][k] = sqrt(M_therm(x_array[0], Mhessian) / beta)
            * distribution(generator);
          // 希望这样采样可可以满足位置和动量的联合分布
        }
      }
    }
  }
  cout << endl;
  cout << "temperature of the system is " << T << ", ";
  cout << "and the setting temperature is " << 1 / beta << endl;
  /*
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos) <<
   T << "  "; OutFile1 << std::setiosflags(std::ios::scientific |
   std::ios::showpos) << (double)1 / beta << "  "; OutFile1 <<
   std::setiosflags(std::ios::scientific | std::ios::showpos) << K << "  ";
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos) <<
   V << "  "; OutFile1 << std::setiosflags(std::ios::scientific |
   std::ios::showpos) << K + V << endl;
    */
    // OutFile1.close();
    // fmt1.clear();
  cout << "-------------------END--PIMD-----------------" << endl;
  cout << "-----------BEIGN--REAL--TIME--DYNAMICS-------" << endl;
  for (int k = 0; k < n_x * n_p; k++) {
    // k: int 表示当前演化的轨线
    // 首先要初始化位置和动量
    if (k % 1000 == 0) {
      t_end = time(NULL);
      printf(
        "in classical HD, in %dth trajactory, time has been used = %.1f "
        "s\n",
        k, difftime(t_end, t_start));
    }
    for (int n = 0; n < N; n++) {
      x_array[0][n] = traj_x_init[k][n];
      p_array[0][n] = traj_p_init[k][n];
    }
    for (int j = 0; j < steps; j++) {
      if (j % step_leap == 0) {  // 每隔step_leap步储存一个构型
        for (int r = 0; r < N; r++) {
          x_list_traj[k][j / step_leap][r] = x_array[0][r];
          p_list_traj[k][j / step_leap][r] = p_array[0][r];
          t_list[j / step_leap] = j * dt_HD;
        }
      }
      Velocity_Verlet(x_array[0], p_array[0], nabla_V_potential, dt_HD);
    }
  }
  // 计算时间关联函数
  stringstream fmt1;  // 通过文件存储时间关联函数的数据
  fmt1 << "1-D_quartic_8_LSCIVR_x_square"
    << ".txt";
  ofstream OutFile1(fmt1.str());
  printf("need some time to calculate correlation function\n");
  for (int t = 0; t < t_steps_max; t++) {
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
      << t_list[t] << "  ";
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
      << x_square_std_corr_func_(0, t, Mhessian) << "  \n";
  }
  OutFile1.close();
  fmt1.clear();
  cout << "-------------------END--LSCIVR-----------------" << endl;
}

int main() {

  calculate(nabla_quartic_potential, Mhessian_quartic);
  int i;  // 让程序不要运行完就退出
  cin >> i;
  return 0;
}