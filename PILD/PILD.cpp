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

const int N = 3;            // 系统的自由度
const int P = 128;          // beads数量
const int steps = 200000;   // 演化的总步数
const int step_leap = 400;  // 跳过step_leap步记录一次位置和动量
const double dt =
    0.1 / sqrt(P);  // 时间步长 对于BAOAB_PILD算法，时间步长最长为2/omega_P
const double m = 1;       // 质量矩阵的对角元
const double beta = 0.5;  // 温度的倒数
const double omega_P = sqrt(P) / beta;
const double gama_ad = 0.1;  // 绝热参数
const double omega_ad =
    sqrt(P / gama_ad) / beta;  // 加入绝热参数后的内势力角频率
const int n_traj = 50;  // 计算时间关联函数时所用轨线的数目
const int PIMD_steps = 10000;  // 选取初始值时PIMD演化的步数

// 定义一些BAOAB_PILD函数常用的常数
const double c1_PILD = exp(-omega_ad * dt);
const double c2_PILD = sqrt((1 - c1_PILD * c1_PILD) / beta);
const double d1_PILD = cos(omega_ad * dt / 2);  // 演化内势力
const double d2_PILD = sin(omega_ad * dt / 2);

// 定义一些BAOAB函数常用的常数
const double c1_PIMD = exp(-omega_P * dt);
const double c2_PIMD = sqrt((1 - c1_PIMD * c1_PIMD) / beta);
const double d1_PIMD = cos(omega_P * dt / 2);  // 演化内势力
const double d2_PIMD = sin(omega_P * dt / 2);

// 定义一些与Morse势有关的常数
const double De = 20;      // 势阱深度
const double alpha = 2.0;  // 势阱宽度倒数
const double r_eq = 2.5;   // 平衡位置

double M_thermal[3][3] = {{1, 0, 0},
                          {0, 1, 0},
                          {0, 0, 1}};  // 热质量矩阵 暂时不会算，假设一个
double M_thermal_inv[3][3] = {{1, 0, 0},
                              {0, 1, 0},
                              {0, 0, 1}};  // 对于谐振子的情形
double M_thermal_M_inv[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
double M_M_thermal_inv[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
double M[N];      // 质量矩阵 1-D array
double M_inv[N];  // 质量矩阵的逆 1-D array

double x_array[P][N];     // 储存所有beads位置
double p_array[P][N];     // 储存所有beads动量
double nabla_V[N] = {0};  // 储存单个beads在Cartesian坐标下的势能梯度
double nabla_phi[P][N] = {0};       // 每个bead在staging坐标下的受力
double traj_x[n_traj][P][N] = {0};  // 储存选定为初始值的平衡构型
double traj_p[n_traj][P][N] = {0};
double config_list[n_traj][2][steps][N] = {0};  // 储存实时间动力学的构型
const int t_steps_max = steps / step_leap - 50;
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
Cartesian坐标下Morse势
x: doubel 1-D array Cartesian坐标
N: int 系统的自由度，在本例子中为3
return: double 对应坐标下的势能
*/
double Morse_potential(double (&x)[N]) {
  double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  double q = (1 - exp(-alpha * (r - r_eq)));
  return De * q * q;
}

/*
Cartesian坐标下真实势能的梯度
x: double 1-D array Cartesian坐标
N: int 系统的自由度
不返回值，修改之前定义的储存单个bead势能梯度的变量
*/
void nabla_Morse_potential(double (&x)[N]) {
  double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  double exp_a = exp(-alpha * (r - r_eq));
  double fractor = 2 * De * alpha * (1 - exp_a) * exp_a / r;
  for (int j = 0; j < N; j++) {
    nabla_V[j] = fractor * x[j];
  }
}

double Harmonic_potential(double (&x)[N]) {
  return 0.5 * (x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
}

void nabla_Harmonic_potential(double (&x)[N]) {
  for (int j = 0; j < N; j++) {
    nabla_V[j] = x[j];
  }
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
  return T / (N * P);
}

/*
kubo形式的动量关联函数
t_steps: int < steps / step_leap - cnt 计算t_steps时间节点上的平均值
cnt: int 每条轨线上计算平均值的次数
*/
double p_kubo_corr_func(int t_steps, int cnt = 30) {
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
    }
  }
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
      p_init_array[j][k] = para * distribution(generator);
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

  for (int k = 0; k < N; k++)  // 随机选取PIMD初始值
  {
    for (int j = 0; j < P; j++) {
      x_array[j][k] = distribution(generator);
      p_array[j][k] = distribution(generator);
    }
  }
  // 用于储存数据
  /*
stringstream fmt1;
fmt1 << "quartic_P_128"
     << ".txt";
ofstream OutFile1(fmt1.str());
*/
  // 用于判断PIMD是否达到平衡
  double T = 0;
  cout << "------------------START--PIMD------------------" << endl;
  cout << "the number of beads is " << P << endl;
  cout << "the setting temperature is " << (double)1 / beta << endl;
  cout << "total steps is " << steps << endl;
  for (int i = 0; i < PIMD_steps; i++) {
    BAOAB(x_array, p_array, nabla_Harmonic_potential, beta, omega_P, c1_PILD,
          c2_PILD, d1_PILD, d2_PILD);
    T += temperature(p_array) / steps;
    if (i % 2000 == 0) {
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
  cout << "and the setting temperature is " << (double)1 / beta << endl;
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
    // 更新x_array中的数据
    for (int r = 0; r < N; r++) {
      for (int s = 0; s < P; s++) {
        x_array[s][r] = traj_x[k][s][r];
        p_array[s][r] = traj_p[k][s][r];
      }
    }
    // 将第0个bead的动量换掉
    for (int r = 0; r < N; r++) {
      p_array[0][r] = p_init_array[k][r];
    }
    // 开始每一条轨线的演化
    for (int j = 0; j < steps; j++) {
      int cnt = 0;
      if (j % step_leap == 0) {
        cnt += 1;
        for (int r = 0; r < N; r++) {
          x_list_traj[k][cnt][r] = x_array[0][r];
          p_list_traj[k][cnt][r] = p_array[0][r];
          t_list[cnt] = j * dt;
        }
      }
      BAOAB_PILD(x_array, p_array, nabla_Harmonic_potential, beta, omega_ad,
                 c1_PILD, c2_PILD, d1_PILD, d2_PILD);
      if (j % 10000 == 0) {
        t_end = time(NULL);
        printf(
            "in PILD, in %dth trajactory ,step = %d, time has been used = %.3f "
            "s\n",
            k, j, difftime(t_end, t_start));
      }
    }
  }

  int i;  // 让程序不要运行完就退出
  cin >> i;
  return 0;
}