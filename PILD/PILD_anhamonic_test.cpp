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
const int P = 500;         // beads数量
const int steps = 250000;  // 演化的总步数
const int step_leap =
    25;  // 跳过step_leap步记录一次位置和动量 计算时间关联函数的最小步长
const double dt =
    0.001 / sqrt(P);  // 时间步长 对于BAOAB_PILD算法，时间步长最长为2/omega_P
const double m = 1;     // 质量矩阵的对角元
const double beta = 8;  // 温度的倒数
const double omega_P = sqrt(P) / beta;
const double gama_ad = 10e-6;  // 绝热参数
const double omega_ad =
    sqrt(P / gama_ad) / beta;  // 加入绝热参数后的内势力角频率
const int n_traj = 100;  // 计算时间关联函数时所用轨线的数目
const int PIMD_steps = 500000;  // 选取初始值时PIMD演化的步数
const int n_atom = 1;

double M_thermal[N][N] = {{3.93462}};
double M_thermal_inv[N][N] = {{0.25415}};
double M_thermal_M_inv[N][N] = {{3.93462}};
double M_M_thermal_inv[N][N] = {{0.25415}};
double M[N];      // 质量矩阵 1-D array
double M_inv[N];  // 质量矩阵的逆 1-D array

// 尽量不要动下面定义的变量
int t_steps_max = steps / step_leap - 100;
double x_array[P][N];     // 储存所有beads位置
double p_array[P][N];     // 储存所有beads动量
double nabla_V[N] = {0};  // 储存单个beads在Cartesian坐标下的势能梯度
double nabla_phi[P][N] = {0};       // 每个bead在staging坐标下的受力
double traj_x[n_traj][P][N] = {0};  // 储存选定为初始值的平衡构型
double traj_p[n_traj][P][N] = {0};
double x_list_traj[n_traj][steps / step_leap][N];  // 所有轨线不同时刻的位置
double p_list_traj[n_traj][steps / step_leap][N];  // 所有轨线不同时刻的动量
double t_list[steps / step_leap];  // 记录每个构型对应的演化时间
double p_init_array[n_traj][N];

// 设置随机数种子
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> distribution(0.0, 1.0);

void print_vector_2(double (&vec)[P][N]) {
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
void staging_transf(double (&x)[P][N]) {
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
void staging_transf_inv(double (&s)[P][N]) {
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
double quartic_potential(double (&x)[N]) {
  return 0.25 * x[0] * x[0] * x[0] * x[0];
}
/*
不对称非简谐势能
x: double 1-D array Cartesian坐标
retrun: double 对应坐标下的势能
*/
double asy_aharmonic_potential(double (&x)[N]) {
  double x_2 = x[0] * x[0];
  return x_2 - 0.1 * x_2 * x[0] + 0.1 * x_2 * x_2;
}

/*
Cartesian坐标下四次方势的导数
x: double 1-D array Cartesian坐标
return: double 1-D array 势能的梯度
*/
void nabla_quartic_potential(double (&x)[N]) {
  nabla_V[0] = x[0] * x[0] * x[0];
}

/*
Cartesian坐标下非对称非简谐势能的导数
x: double 1-D array Cartesian坐标
return: double 1-D array 势能的梯度
*/
void nabla_asy_aharmonic_potential(double (&x)[N]) {
  nabla_V[0] = 2.0 * x[0] - 0.3 * x[0] * x[0] + 0.4 * x[0] * x[0] * x[0];
}

/*
staging坐标下的所有beads的外势
s_chain: double 2-D array 所有beads的staging坐标
return: double 势能
*/
double phi_potential(double (&s_chain)[P][N],
                     double (*potential)(double (&)[N])) {
  staging_transf_inv(s_chain);
  double phi = 0;
  for (int j = 0; j < P; j++) {
    phi += potential(s_chain[j]);
  }
  staging_transf(s_chain);
  return (double)phi / P;
}

/*
staging坐标下所有beads的势能梯度，在求解Langevin方程中直接使用
s_chain: double 2-D array 所有beads的staging坐标
return: douible 2-D array 所有beads在staging坐标下的势能梯度
注意这里没有修改传入的数组
*/
void nabla_phi_potential(double (&s_chain)[P][N],
                         void (*nabla_potential)(double (&)[N])) {
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
double kinetic_energy_viral(double (&s_chain)[P][N],
                            double beta,
                            void (*nabla_potential)(double (&)[N])) {
  staging_transf_inv(s_chain);
  double x_c[N] = {0};  // 构造所有beads质心的坐标
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
const double c1_PIMD = exp(-omega_P * dt);
const double c2_PIMD = sqrt((1 - c1_PIMD * c1_PIMD) / beta);
const double d1_PIMD = cos(omega_P * dt / 2);  // 演化内势力
const double d2_PIMD = sin(omega_P * dt / 2);
/*
以BAOAB的方式演化Langevin方程，均为staging坐标
x_chain: double 2-D array 表示前一个时刻所有beads的位置
p_chain: double 2-D array 表示前一个时刻所有beads的动量
没有返回值 直接修传入的数组
*/
void BAOAB(double (&s_chain)[P][N],
           double (&p_chain)[P][N],
           void (*nabla_potential)(double (&)[N]),
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
    for (int k = 0; k < N; k++) {
      p_chain[j][k] += -0.5 * dt * nabla_phi[j][k];
    }
  }

  // 其次演化<半步>A内势力项
  for (int k = 0; k < N; k++) {
    s_chain[0][k] += +0.5 * dt / M[k] * p_chain[0][k];
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      var = s_chain[j][k] * d1_PIMD +
            d2_PIMD * j / ((j + 1) * M[k] * omega_P) * p_chain[j][k];
      p_chain[j][k] = -omega_P * d2_PIMD * (j + 1) / j * M[k] * s_chain[j][k] +
                      d1_PIMD * p_chain[j][k];
      s_chain[j][k] = var;
    }
  }

  // 然后演化一步O控温项
  for (int k = 0; k < N; k++) {
    p_chain[0][k] = c1_PIMD * p_chain[0][k] +
                    c2_PIMD * sqrt(M[k]) * distribution(generator);
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      p_chain[j][k] =
          c1_PIMD * p_chain[j][k] +
          c2_PIMD * sqrt((double)(j + 1) / j * M[k]) * distribution(generator);
    }
  }

  // 再演化半步A内势力项
  for (int k = 0; k < N; k++) {
    s_chain[0][k] += +0.5 * dt / M[k] * p_chain[0][k];
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      var = s_chain[j][k] * d1_PIMD +
            d2_PIMD * j / ((j + 1) * M[k] * omega_P) * p_chain[j][k];
      p_chain[j][k] = -omega_P * d2_PIMD * (j + 1) / j * M[k] * s_chain[j][k] +
                      d1_PIMD * p_chain[j][k];
      s_chain[j][k] = var;
    }
  }

  // 最后演化半步B动量项
  nabla_phi_potential(s_chain, nabla_potential);
  for (int j = 0; j < P; j++) {
    for (int k = 0; k < N; k++) {
      p_chain[j][k] += -0.5 * dt * nabla_phi[j][k];
    }
  }
}

// 定义一些BAOAB_PILD函数常用的常数
const double c1_PILD = exp(-omega_ad * dt);
const double c2_PILD = sqrt((1 - c1_PILD * c1_PILD) / beta);
const double d1_PILD = cos(omega_ad * dt / 2);  // 演化内势力
const double d2_PILD = sin(omega_ad * dt / 2);
/*
以BAOAB的方式演化PILD的方程
引入了绝热参数的概念，认为有效力只用一步采样就可以得到比较准确的结果？
x_chain: double 2-D array 表示前一个时刻所有beads的位置
p_chain: double 2-D array 表示前一个时刻所有beads的动量
没有返回值 直接修传入的数组
*/
void BAOAB_PILD(double (&s_chain)[P][N],
                double (&p_chain)[P][N],
                void (*nabla_potential)(double (&)[N]),
                double beta,
                double omega_ad,
                double c1_PILD,
                double c2_PILD,
                double d1_PILD,
                double d2_PILD) {
  double var = 0;  // 由于储存计算过程中的中间变量

  // 首先演化半步B动量项 第一个bead要和其他bead分开演化
  nabla_phi_potential(s_chain, nabla_potential);
  for (int k = 0; k < N; k++) {
    double force = 0;
    for (int r = 0; r < N; r++) {
      force += M_thermal_M_inv[k][r] * nabla_phi[0][r];
    }
    p_chain[0][k] -= 0.5 * dt * force;
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      p_chain[j][k] -= -0.5 * dt * nabla_phi[j][k];
    }
  }

  // 其次演化<半步>A内势力项
  for (int k = 0; k < N; k++) {
    s_chain[0][k] += +0.5 * dt / M[k] * p_chain[0][k];
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      var = s_chain[j][k] * d1_PILD +
            d2_PILD * j / ((j + 1) * M[k] * omega_ad) * p_chain[j][k];
      p_chain[j][k] = -omega_ad * d2_PILD * (j + 1) / j * M[k] * s_chain[j][k] +
                      d1_PILD * p_chain[j][k];
      s_chain[j][k] = var;
    }
  }

  // 然后演化一步O控温项
  // 第一个bead不控温
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      p_chain[j][k] =
          c1_PILD * p_chain[j][k] +
          c2_PILD * sqrt((double)(j + 1) / j * M[k]) * distribution(generator);
    }
  }

  // 再演化半步A内势力项
  for (int k = 0; k < N; k++) {
    s_chain[0][k] += +0.5 * dt / M[k] * p_chain[0][k];
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      var = s_chain[j][k] * d1_PILD +
            d2_PILD * j / ((j + 1) * M[k] * omega_ad) * p_chain[j][k];
      p_chain[j][k] = -omega_ad * d2_PILD * (j + 1) / j * M[k] * s_chain[j][k] +
                      d1_PILD * p_chain[j][k];
      s_chain[j][k] = var;
    }
  }

  // 最后演化半步B动量项
  nabla_phi_potential(s_chain, nabla_potential);
  for (int k = 0; k < N; k++) {
    double force = 0;
    for (int r = 0; r < N; r++) {
      force += M_thermal_M_inv[k][r] * nabla_phi[0][r];
    }
    p_chain[0][k] -= 0.5 * dt * force;
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      p_chain[j][k] = p_chain[j][k] - 0.5 * dt * nabla_phi[j][k];
    }
  }
}

/*
计算给定构型的温度
p_chain: double 2-D array 所有beads的动量
return: double 温度
*/
double temperature(double (&p_chain)[P][N]) {
  double T = 0;
  for (int k = 0; k < N; k++) {
    T += p_chain[0][k] * M_inv[k] * p_chain[0][k];
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      T += p_chain[j][k] * j * M_inv[k] * p_chain[j][k] / (j + 1);
    }
  }
  return T / (N * P );
}

/*
kubo形式的动量关联函数
t_steps: int < steps / step_leap - cnt 计算t_steps时间节点上的平均值
cnt: int 每条轨线上计算平均值的次数
*/
double p_kubo_corr_func(int t_steps, int cnt = 50) {
  if (cnt + t_steps > steps / step_leap) {
    throw("error in p_kubo_func");
    return 0;
  }
  double ave = 0;
  for (int j = 0; j < n_traj; j++) {  // 对所有轨线平均
    for (int m = 0; m < cnt; m++) {   // 对同一条轨线作时间平均
      double tot = 0;
      for (int r = 0; r < N; r++) {
        double var = 0;
        for (int s = 0; s < N; s++) {
          var += M_M_thermal_inv[r][s] * p_list_traj[j][t_steps + m][s];
        }
        tot += p_list_traj[j][m][r] * var;
      }
      ave += tot;
    }
  }
  return ave / (cnt * n_traj);
}
/*
kubo形式的位置关联函数
t_steps: int < steps / step_leap - cnt 计算t_steps时间节点上的平均值
cnt: int 每条轨线上计算平均值的次数
有点问题，应该怎么写
*/
double x_kubo_corr_func(int t_steps, int cnt = 30) {
  if (cnt + t_steps > steps / step_leap) {
    throw("error in p_kubo_func");
    return 0;
  }
  return 0;
}
const double e[n_atom] = {1};  // 电荷
/*
kubo形式的偶极导数关联函数
*/
double dip_drv_kubo_corr_func(int t_steps, int cnt = 30) {
  if (cnt + t_steps > steps / step_leap) {
    throw("error in p_kubo_func");
    return 0;
  }
  double ave = 0;
  for (int j = 0; j < n_traj; j++) {  // 对所有轨线平均
    for (int m = 0; m < cnt; m++) {   // 对同一条轨线作时间平均
      double tot = 0;

      //有大问题，改日再写
      for (int r = 0; r < N; r++) {
        double var = 0;
        for (int s = 0; s < N; s++) {
          var += M_M_thermal_inv[r][s] * p_list_traj[j][t_steps + m][s];
        }
        tot += p_list_traj[j][m][r] * var;
      }
      ave += tot;
    }
  }
  return ave / (cnt * n_traj);
}

/*
根据M_therm生成对应的初始时刻动量的高斯分布
M_therm: double 2-D array N*N 热质量矩阵
beta: double
n_traj: int 轨线数目
p_init_array: double 2-D array n_traj*N 初始的动量分布
*/
void generate_p_init(double (&M_therm_inv)[N][N],
                     double (&p_init_array)[n_traj][N]) {
  // 这里要写入对角化M_therm与坐标变换的过程
  // mkl中应该有向量化的生成随机数的函数？
  double para = 1 / sqrt(beta);
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < n_traj; j++) {
      p_init_array[j][k] =
          para * distribution(generator) * sqrt(M_thermal[k][k]);
    }
  }
}

int main() {
  for (int i = 0; i < N; i++)  // 初始化质量矩阵和其逆矩阵，假设均为对角矩阵
  {
    M[i] = m;
    M_inv[i] = 1 / m;
  }

  time_t t_start, t_end;  // 记录程序运行时间
  t_start = time(NULL);

  for (int k = 0; k < N; k++)  //
  {
    for (int j = 0; j < P; j++) {
      x_array[j][k] = 0;
      p_array[j][k] = 0;
    }
  }
  staging_transf(x_array);
  // 用于判断PIMD是否达到平衡
  double T = 0;
  cout << "------------------START--PIMD------------------" << endl;
  cout << "the number of beads is " << P << endl;
  cout << "the setting temperature is " << beta << endl;
  cout << "total steps is " << PIMD_steps << endl;

  for (int i = 0; i < PIMD_steps; i++) {
    BAOAB(x_array, p_array, nabla_quartic_potential, beta, omega_P, c1_PILD,
          c2_PILD, d1_PILD, d2_PILD);
    T += temperature(p_array) / PIMD_steps;
    if (i % 20000 == 0) {
      t_end = time(NULL);
      printf("in PIMD, step = %d, time has been used = %.3f s\n", i,
             difftime(t_end, t_start));
    }
    if (i >= PIMD_steps - n_traj) {
      for (int k = 0; k < N; k++) {
        for (int j = 0; j < P; j++) {
          traj_x[i - PIMD_steps + n_traj][j][k] = x_array[j][k];
          traj_p[i - PIMD_steps + n_traj][j][k] = p_array[j][k];
        }
      }
    }  // 储存了空间构型，相当于对对积分中空间构型的采样
  }
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
  generate_p_init(M_thermal_inv, p_init_array);

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
    // 将第0个bead的动量换掉，使动量分布满足高斯分布
    for (int r = 0; r < N; r++) {
      p_array[0][r] = p_init_array[k][r];
    }
    // 开始每一条轨线的演化
    for (int j = 0; j < steps; j++) {
      if (j % step_leap == 0) {  // 每隔step_leap步储存一个构型
        for (int r = 0; r < N; r++) {
          x_list_traj[k][j / step_leap][r] = x_array[0][r];
          p_list_traj[k][j / step_leap][r] = p_array[0][r];
          t_list[j / step_leap] = j * dt;
        }
      }
      BAOAB_PILD(x_array, p_array, nabla_quartic_potential, beta, omega_ad,
                 c1_PILD, c2_PILD, d1_PILD, d2_PILD);
      if (j % 100000 == 0) {
        t_end = time(NULL);
        printf(
            "in PILD, in %dth trajactory ,step = %d, time has been used = %.1f "
            "s\n",
            k, j, difftime(t_end, t_start));
      }
    }
  }
  // 计算时间关联函数
  stringstream fmt1;  // 通过文件存储时间关联函数的数据
  fmt1 << "1-D_quartic"
       << ".txt";
  ofstream OutFile1(fmt1.str());
  for (int t = 0; t < t_steps_max; t++) {
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
             << t_list[t] << "  ";
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
             << p_kubo_corr_func(t) << "  \n";
  }
  OutFile1.close();
  fmt1.clear();
  cout << "-------------------END--PILD-----------------" << endl;

  int i;  // 让程序不要运行完就退出
  cin >> i;
  return 0;
}