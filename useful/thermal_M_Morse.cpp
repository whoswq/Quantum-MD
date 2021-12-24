#include <math.h>
#include <mkl.h>
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
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

// 定义一些常数为全局变量，所有数组的定义都放在main函数中
const int N = 3;                  // 系统的自由度
const int P = 200;                // beads数量
const int steps = 10000000;       // 演化的总步数
const double dt = 0.2 / sqrt(P);  // 时间步长
const int step_leap = 25;         // 计算热质量矩阵的时间间隔
const double m = 1728.16;         // 质量矩阵的对角元
const double k_B = 3.16671e-6;    // Boltzman常数
const double tem = 100;           // 温度
const double beta = (double)1 / (k_B * tem);  // 温度的倒数
const double omega_P = sqrt(P) / beta;

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
void staging_transf(double(*x)) {
  // 第一个bead的坐标不用动
  for (int i = 1; i < P; i++) {
    for (int j = 0; j < N; j++) {
      x[i * N + j] = x[i * N + j] -
                     (double)(i * x[((i + 1) % P) * N + j] + x[j]) / (i + 1);
    }
  };
}

/*
staging逆变换，将staging坐标转换成为Cartesian坐标
s: double 2-D array 所有beads的staging坐标
没有返回值
*/
void staging_transf_inv(double(*s)) {
  // 第一个bead不用动
  for (int j = P - 1; j > 0; j--) {
    for (int i = 0; i < N; i++) {
      s[j * N + i] = s[j * N + i] +
                     (double)(j * s[((j + 1) % P) * N + i] + s[i]) / (j + 1);
    }
  }
}

// 定义一些与Morse势有关的常数
/* 参考2011年刘老师的文章
Insights in quantum dynamical effects in the infrared spectroscopy of liquid
water from a semiclassical study with an ab initio-based flexible and
polarizable force field J. Chem. Phys. 135, 244503 (2011);
https://doi.org/10.1063/1.3670960
*/
const double De = 0.1850;     // 势阱深度
const double alpha = 1.2102;  // 势阱宽度倒数
const double r_eq = 1.7799;   // 平衡位置
/*
Cartesian坐标下Morse势
x: doubel 1-D array Cartesian坐标
N: int 系统的自由度，在本例子中为3
return: double 对应坐标下的势能
*/
double Morse_potential(double(*x)) {
  double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  double q = (1 - exp(-alpha * (r - r_eq)));
  return De * q * q;
}
/*
Cartesian坐标下真实势能的梯度
x: double 1-D array Cartesian坐标
N: int 系统的自由度
*/
void nabla_Morse_potential(double(*x), double(*nabla_V)) {
  double r = sqrt(cblas_ddot(N, x, 1, x, 1));
  double exp_a = exp(-alpha * (r - r_eq));
  double fractor = 2 * De * alpha * (1 - exp_a) * exp_a / r;
  for (int j = 0; j < N; j++) {
    nabla_V[j] = fractor * x[j];
  }
}

/*
计算质量加权的hessian矩阵，通过中心差分的方法计算
x: pointer double 1-D array cartesian坐标
MHessian: pointer double 质量加权的hessian矩阵
M: pointer double 质量矩阵
d: double 有限差分时的步长
*/
void MHessian_Morse(double(*x),
                    double(*MHessian),
                    double(*M),
                    double d = 0.0002) {}

/*
staging坐标下的所有beads的外势
s_chain: double 2-D array 所有beads的staging坐标
return: double 势能
*/
double phi_potential(double (&s_chain)[P][N],
                     double (*potential)(double (&)[N])) {
  staging_transf_inv(s_chain);  // 這時候s_chain中儲存cartesian坐標
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
void nabla_Phi_potential(double (&s_chain)[P][N],
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
                            double (*nabla_potential)(double (&)[N])) {
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
const double c1 = exp(-omega_P * dt);
const double c2 = sqrt((1 - c1 * c1) / beta);
const double d1 = cos(omega_P * dt / 2);
const double d2 = sin(omega_P * dt / 2);
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
           double c1,
           double c2,
           double d1,
           double d2) {
  double var = 0;  // 由于储存计算过程中的中间变量
  // 首先演化半步B动量项
  nabla_Phi_potential(s_chain, nabla_potential);
  for (int j = 0; j < P; j++) {
    for (int k = 0; k < N; k++) {
      p_chain[j][k] = p_chain[j][k] - 0.5 * dt * nabla_phi[j][k];
    }
  }

  // 其次演化<半步>A内势力项
  for (int k = 0; k < N; k++) {
    s_chain[0][k] = s_chain[0][k] + 0.5 * dt / M[k] * p_chain[0][k];
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      var = s_chain[j][k] * d1 +
            d2 * j / ((j + 1) * M[k] * omega_P) * p_chain[j][k];
      p_chain[j][k] = -omega_P * d2 * (j + 1) / j * M[k] * s_chain[j][k] +
                      d1 * p_chain[j][k];
      s_chain[j][k] = var;
    }
  }

  // 然后演化一步O控温项
  for (int k = 0; k < N; k++) {
    p_chain[0][k] =
        c1 * p_chain[0][k] + c2 * sqrt(M[k]) * distribution(generator);
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      p_chain[j][k] =
          c1 * p_chain[j][k] +
          c2 * sqrt((double)(j + 1) / j * M[k]) * distribution(generator);
    }
  }

  // 再演化半步A内势力项
  for (int k = 0; k < N; k++) {
    s_chain[0][k] = s_chain[0][k] + 0.5 * dt / M[k] * p_chain[0][k];
  }
  for (int j = 1; j < P; j++) {
    for (int k = 0; k < N; k++) {
      var = s_chain[j][k] * d1 +
            d2 * j / ((j + 1) * M[k] * omega_P) * p_chain[j][k];
      p_chain[j][k] = -omega_P * d2 * (j + 1) / j * M[k] * s_chain[j][k] +
                      d1 * p_chain[j][k];
      s_chain[j][k] = var;
    }
  }

  // 最后演化半步B动量项
  nabla_Phi_potential(s_chain, nabla_potential);
  for (int j = 0; j < P; j++) {
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
用于计算热质量矩阵
x: double 1-D array N 必须是cartesian坐标中的位置
*/
double M_therm(double (&x)[N], double (*Mhessian)(double (&)[N])) {
  double u = beta / 2.0;
  double w_2 = Mhessian(x);
  if (w_2 >= 0) {
    u = u * sqrt(w_2);
    return u / tanh(u);
  } else {
    u = u * sqrt(-w_2);
    return tanh(u) / u;
  }
}

/*
经过测试之后发现如果只用M_thermal的primitiv estimator，收敛性较差
这里仿照计算势能的方法写一个计算在所有beads上平均的M_thermal
x: double 2-D array P*N 所有beads的直角坐标
*/
double M_therm_average(double (*Mhessian)(double (&)[N])) {
  // 在演化的过程中x_array中储存着beads的staging坐标，首先应变换回Cartesian坐标
  staging_transf_inv(x_array);
  double M_t = 0;
  for (int j = 0; j < P; j++) {
    M_t += Mhessian(x_array[j]);
  }
  // 最后还要变回去
  staging_transf(x_array);
  return M_t / P;
}

int main() {
  double* M;           // 质量矩阵 1-D array
  double* MHessian;    // 质量加权的Hessian矩阵
  double* M_inv;       // 质量矩阵的逆 1-D array
  double* M_inv_sqrt;  // 质量矩阵的平方根的逆矩阵
  double* x_array;     // 储存所有beads位置
  double* p_array;     // 储存所有beads动量
  double* nabla_V;     // 储存势能的梯度
  double* nabla_phi;   // 势能的梯度
  // 为指针申请内存
  M = (double*)mkl_malloc(N * N * sizeof(double), 128);
  M_inv_sqrt = (double*)mkl_malloc(N * N * sizeof(double), 128);
  MHessian = (double*)mkl_malloc(N * N * sizeof(double), 128);
  M_inv = (double*)mkl_malloc(N * N * sizeof(double), 128);
  x_array = (double*)mkl_malloc(P * N * sizeof(double), 128);
  p_array = (double*)mkl_malloc(P * N * sizeof(double), 128);
  nabla_V = (double*)mkl_malloc(1 * N * sizeof(double), 128);
  nabla_phi = (double*)mkl_malloc(P * N * sizeof(double), 128);

  for (int i = 0; i < N; i++) {  // 初始化质量矩阵和逆矩阵
    for (int j = 0; j < N; j++) {
      if (i == j) {
        M[i * N + j] = m;
        M_inv[i * N + j] = (double)1 / m;
        M_inv_sqrt[i * N + j] = sqrt((double)1 / m);
      } else {
        M[i * N + j] = 0;
        M_inv[i * N + j] = 0;
        M_inv_sqrt[i * N + j] = 0 ;
      }
    }
  }

  time_t start, end;  // 定义储存时间的变量
  start = time(NULL);
  // 初始化位置和动量，随机选取
  for (int j = 0; j < P; j++) {
    for (int k = 0; k < N; k++) {
      x_array[j * N + k] = distribution(generator);
      p_array[j * N + k] = distribution(generator);
    }
  }
  // 用于储存数据
  /*
  stringstream fmt1;
  fmt1 << "quartic_P_128"
       << ".txt";
  ofstream OutFile1(fmt1.str());
  */
  double T = 0;    // 判断是否达到平衡
  double M_t = 0;  // 保存热质量矩阵
  cout << "------------------START------------------" << endl;
  cout << "the number of beads is " << P << endl;
  cout << "the setting temperature is " << (double)1 / beta << endl;
  cout << "total steps is " << steps << endl;
  for (int i = 0; i < steps; i++) {
    BAOAB(x_array, p_array, nabla_asy_aharmonic_potential, beta, omega_P, c1,
          c2, d1, d2);
    T += temperature(p_array) / steps;
    if (i % 100000 == 0) {
      end = time(NULL);
      cout << "step " << i << ", time  " << (end - start) << " s" << endl;
    }
    if (i % step_leap == 0) {
      M_t += M_therm_average(Mhessian_quartic) / steps * step_leap;
    }
  }
  cout << "temperature of the system is " << T << ", ";
  cout << "the setting temperature is " << (double)1 / beta << endl;
  cout << "M_thermal = " << M_t << endl;
  cout << "-------------------END-------------------" << endl;

  // 为了不让程序运行结束后直接退出
  int i;
  cin >> i;
  return 0;
}