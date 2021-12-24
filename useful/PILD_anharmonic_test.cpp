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

const int N = 1;           // 系统的自由度
const int P = 200;         // beads数量
const int steps_PILD = 2000000;  // 演化的总步数
const int step_leap = 10;  // 跳过step_leap步记录一次位置和动量 计算时间关联函数的最小步长
const double beta = 0.1;                         // 温度的倒数
const double dt_PILD = 0.0001 / sqrt(P);  // 时间步长 对于BAOAB_PILD算法，时间步长最长为2/omega_P
const double dt_PIMD = 0.1 / sqrt(P);  // PIMD中的时间步长
const double m = 1;                            // 质量矩阵的对角元
const double omega_P = sqrt(P) / beta;
const double gama_ad = 1e-4;  // 绝热参数
const double omega_ad = sqrt(P / gama_ad) / beta;  // 加入绝热参数后的内势力角频率
const int n_x = 40;
const int n_p = 3;
const int n_traj = n_x * n_p;  // 计算时间关联函数时所用轨线的数目
const int steps_PIMD = 2000000;  // 选取初始值时PIMD演化的步数

double Mtherm[N] = { 1.004 };
double M[N] = { 0 };      // 质量矩阵 1-D array
double M_inv[N] = { 0 };  // 质量矩阵的逆 1-D array

// 尽量不要动下面定义的变量
int t_steps_max = steps_PILD / step_leap - 100;
double x_array[P][N] = { 0 };     // 储存所有beads位置
double p_array[P][N] = { 0 };     // 储存所有beads动量
double nabla_V[N] = { 0 };  // 储存单个beads在Cartesian坐标下的势能梯度
double nabla_phi[P][N] = { 0 };       // 每个bead在staging坐标下的受力
double traj_x[n_traj][P][N] = { 0 };  // 储存选定为初始值的平衡构型
double traj_p[n_traj][P][N] = { 0 };
double x_list_traj[n_traj][steps_PILD / step_leap][N];  // 所有轨线不同时刻的位置
double p_list_traj[n_traj][steps_PILD / step_leap][N];  // 所有轨线不同时刻的动量
double t_list[steps_PILD / step_leap];  // 记录每个构型对应的演化时间
double p_init_array[n_traj][N] = { 0 };

// 设置随机数种子
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> distribution(0.0, 1.0);

void print_vector_2(double(&vec)[P][N]) {
  cout << "[";
  for (int i = 0; i < P; i++) {
    cout << "[";
    for (int j = 0; j < N; j++) {
      cout << vec[i][j] << ",";
    }
    cout << "], ";
  }
  cout << "]" << endl;
}

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
  if (abs(w_2) < 1e-12) {
    return 1.0;
  }
  else {
    if (w_2 >= 0) {
      u = u * sqrt(w_2);
      return u / tanh(u);
    }
    else {
      u = u * sqrt(-w_2);
      return tanh(u) / u;
    }
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
  // 首先将数组初始化 不清楚这一步对效率影响有多大？
  for (int j = 0; j < P; j++) {
    for (int k = 0; k < N; k++) {
      nabla_phi[j][k] = 0.0;
    }
  }
  staging_transf_inv(s_chain);
  // 将staging坐标变换为Cartesian坐标，方便计算势能
  {
    for (int j = 0; j < P; j++) {
      nabla_potential(s_chain[j]);
      for (int k = 0; k < N; k++) {
        // 经过staging_transf_inv之后，s_chain表示对应的Cartesian坐标
        nabla_phi[0][k] += nabla_V[k] / P;
      }
    }
  }
  // 后面的每个bead的势能梯度都需要递推计算
  for (int j = 1; j < P; j++) {
    nabla_potential(s_chain[j]);
    for (int k = 0; k < N; k++) {
      nabla_phi[j][k] =
        nabla_V[k] / P + (double)(j - 1) / j * nabla_phi[j - 1][k];
    }
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
  double var = 0;  // 由于储存计算过程中的中间变量
  // 首先演化半步B动量项
  nabla_phi_potential(s_chain, nabla_potential);
  for (int j = 0; j < P; j++) {
    p_chain[j][0] += -0.5 * dt_PIMD * nabla_phi[j][0];
  }

  // 其次演化<半步>A内势力项
  s_chain[0][0] += +0.5 * dt_PIMD / M[0] * p_chain[0][0];
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

// 定义一些BAOAB_PILD函数常用的常数
const double c1_PILD = exp(-omega_ad * dt_PILD);
const double c2_PILD = sqrt((1 - c1_PILD * c1_PILD) / beta);
const double d1_PILD = cos(omega_ad * dt_PILD / 2);  // 演化内势力
const double d2_PILD = sin(omega_ad * dt_PILD / 2);
/*
以BAOAB的方式演化PILD的方程
引入了绝热参数的概念，认为有效力只用一步采样就可以得到比较准确的结果？
x_chain: double 2-D array 表示前一个时刻所有beads的位置
p_chain: double 2-D array 表示前一个时刻所有beads的动量
nabla_potential: void 返回给定势能的梯度
之后是一些参数

使用平均的M_therm，在计算有效力时不考虑M_therm随位置的变化
*/
void BAOAB_PILD_average(double(&s_chain)[P][N],
  double(&p_chain)[P][N],
  void (*nabla_potential)(double(&)[N]),
  double beta,
  double omega_ad,
  double c1_PILD,
  double c2_PILD,
  double d1_PILD,
  double d2_PILD) {
  double var = 0;  // 由于储存计算过程中的中间变量

  // 首先演化半步B动量项 第一个bead要和其他bead分开演化
  nabla_phi_potential(s_chain, nabla_potential);
  p_chain[0][0] += -0.5 * dt_PILD * Mtherm[0] * M_inv[0] * nabla_phi[0][0];
  // M_therm(s_chain[0], Mhessian)

  for (int j = 1; j < P; j++) {
    p_chain[j][0] += -0.5 * dt_PILD * nabla_phi[j][0];
  }

  // 其次演化<半步>A内势力项
  s_chain[0][0] += +0.5 * dt_PILD / M[0] * p_chain[0][0];
  for (int j = 1; j < P; j++) {
    var = s_chain[j][0] * d1_PILD +
      d2_PILD * j / ((j + 1) * M[0] * omega_ad * gama_ad) * p_chain[j][0];
    p_chain[j][0] = -omega_ad * d2_PILD * (j + 1) / j * M[0] * gama_ad * s_chain[j][0] +
      d1_PILD * p_chain[j][0];
    s_chain[j][0] = var;
  }

  // 然后演化一步O控温项
  // 第一个bead不控温
  for (int j = 1; j < P; j++) {
    p_chain[j][0] =
      c1_PILD * p_chain[j][0] +
      c2_PILD * sqrt((double)(j + 1) / j * M[0] * gama_ad) * distribution(generator);
  }

  // 其次演化<半步>A内势力项 注意绝热参数出现的位置！！！
  s_chain[0][0] += +0.5 * dt_PILD / M[0] * p_chain[0][0];
  for (int j = 1; j < P; j++) {
    var = s_chain[j][0] * d1_PILD +
      d2_PILD * j / ((j + 1) * M[0] * omega_ad * gama_ad) * p_chain[j][0];
    p_chain[j][0] = -omega_ad * d2_PILD * (j + 1) / j * M[0] * gama_ad * s_chain[j][0] +
      d1_PILD * p_chain[j][0];
    s_chain[j][0] = var;
  }

  // 最后演化半步B动量项
  nabla_phi_potential(s_chain, nabla_potential);
  p_chain[0][0] += -0.5 * dt_PILD * Mtherm[0] * M_inv[0] * nabla_phi[0][0];
  for (int j = 1; j < P; j++) {
    p_chain[j][0] += -0.5 * dt_PILD * nabla_phi[j][0];
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

void calculate(void (*nabla_V_potential)(double(&)[N]), double (*Mhessain)(double(&)[N])) {
  for (int i = 0; i < N; i++)  // 初始化质量矩阵和其逆矩阵，假设均为对角矩阵
  {
    M[i] = m;
    M_inv[i] = 1 / m;
  }

  time_t t_start, t_end;  // 记录程序运行时间
  t_start = clock();

  for (int k = 0; k < N; k++)  // 构造PIMD的初始条件
  {
    for (int j = 0; j < P; j++) {
      x_array[j][k] = 0;
      p_array[j][k] = 0;
    }
  }
  staging_transf(x_array);
  double T = 0;  // 用于判断PIMD是否达到平衡
  cout << "------------------START--PIMD------------------" << endl;
  cout << "the number of beads is " << P << endl;
  cout << "the setting temperature is " << beta << endl;
  cout << "total steps is " << steps_PIMD << endl;

  for (int i = 0; i < steps_PIMD; i++) {  // 用PIMD对量子正则分布进行抽样
    BAOAB(x_array, p_array, nabla_V_potential, beta, omega_P, c1_PIMD,
      c2_PIMD, d1_PIMD, d2_PIMD);
    T += temperature(p_array) / steps_PIMD;
    if (i % 20000 == 0) {
      t_end = clock();
      printf("in PIMD, step = %d, time has been used = %.3f s\n", i,
        (double)(t_end - t_start) / CLOCKS_PER_SEC);
    }
    if (i >= steps_PIMD - n_x) {  // 首先在位形空间采样
      for (int k = 0; k < N; k++) {
        for (int j = 0; j < P; j++) {
        if (j == 0){
          for (int np = 0; np < n_p; np++) {
            traj_x[(i - steps_PIMD + n_x) * n_p + np][j][k] = x_array[j][k];
          }
        else {
        for (int np = 0; np < n_p; np++){
        traj_x[(i - steps_PIMD + n_x) * n_p + np][j][k] = 0;}
		 } }
        }
      }
      // 还需要对动量进行采样
      for (int k = 0; k < N; k++) {  // 空间维数
        for (int j = 0; j < P; j++) {
          if (j == 0) {
            for (int np = 0; np < n_p; np++) {  // 对于每一个位置，动量空间的采样
              traj_p[np + (i - steps_PIMD + n_x) * n_p][j][k] = sqrt(M_therm(x_array[0], Mhessain) / beta)
                * distribution(generator);
              // 希望这样采样可可以满足位置和动量的联合分布
            }
          }
          else {
            for (int np = 0; np < n_p; np++) {  // 对于每一个位置，动量空间的采样
              traj_p[np + (i - steps_PIMD + n_x) * n_p][j][k] = 0;
            }
          }
        }
      }
    }
  }
  cout << "temperature of the system is " << T << ", ";
  cout << "and the setting temperature is " << 1 / beta << endl;
  cout << "-------------------END--PIMD-----------------" << endl;
  cout << "-----------BEIGN--REAL--TIME--DYNAMICS-------" << endl;

  for (int k = 0; k < n_traj; k++) {
    // 这里尝试将演化到平衡的PIMD的构型拿过来直接使用
    // k: int 表示当前演化的轨线
    // 更新x_array p_array中的数据
    for (int r = 0; r < N; r++) {
      for (int s = 0; s < P; s++) {
        x_array[s][r] = traj_x[k][s][r];
        p_array[s][r] = traj_p[k][s][r];
      }
    }
    // 开始每一条轨线的演化
    for (int j = 0; j < steps_PILD; j++) {
      if (j % step_leap == 0) {  // 每隔step_leap步储存一个构型
        for (int r = 0; r < N; r++) {
          x_list_traj[k][j / step_leap][r] = x_array[0][r];
          p_list_traj[k][j / step_leap][r] = p_array[0][r];
        }
      }
      BAOAB_PILD_average(x_array, p_array, nabla_V_potential, beta,
        omega_ad, c1_PILD, c2_PILD, d1_PILD, d2_PILD);
      // BAOAB_PILD_1(x_array, p_array, nabla_quartic_potential,
      // Mhessian_quartic,
      //           beta, omega_ad, c1_PILD, c2_PILD, d1_PILD, d2_PILD);
      if (j % 100000 == 0) {
        t_end = clock();
        printf(
          "in PILD, in %dth trajactory ,step = %d, time has been used = %.1f "
          "s\n",
          k, j, (double)(t_end - t_start) / CLOCKS_PER_SEC);
      }
    }
  }
  // 计算时间关联函数
  stringstream fmt1;  // 通过文件存储时间关联函数的数据
  fmt1 << "1-D_quartic_0-1_average_compare"
    << ".txt";
    // compare: 表示将PILD的动量初始值全部设为0
  ofstream OutFile1(fmt1.str());
  for (int t = 0; t < t_steps_max; t++) {
    t_list[t] = t * dt_PILD * step_leap;
  }
  for (int t = 0; t < t_steps_max; t++) {
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
      << t_list[t] << "  ";
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
      << p_kubo_corr_func(t, Mhessain, 10) << "  \n";
  }
  OutFile1.close();
  fmt1.clear();
  cout << "-------------------END--PILD-----------------" << endl;
}

int main() {
  calculate(nabla_quartic_potential, Mhessian_quartic);
  // 更换测试体系时，更改caluclate中的参数、文件最开头的定义（温度），储存的文件名
  int i;  // 让程序不要运行完就退出
  cin >> i;
  return 0;
}
