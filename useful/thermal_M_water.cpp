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
#include <mkl.h>
// using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

// 约定所有单位均为原子单位制
const int dim = 3; // 空间维数，均为3维
const int n_atom = 3;  // 原子个数
const int N = dim * n_atom;
const int P = 256;          // beads数量
const int step_leap = 5;  // 跳过step_leap步记录一次位置和动量 计算时间关联函数的最小步长
const double t_var = 300;                        // 设定温度
const double kb = 1.28065e-23 / 4.35974e-18; // Boltzman常数
const double beta = 1.0 / (kb * t_var);  // 温度的倒数
const double dt_PIMD = 0.25; // PIMD中的时间步长
const double omega_P = sqrt(P) / beta;
const int steps_PIMD = 500000; // 坐标空间抽样时使用的PIMD步数

double hessian[n_atom * dim * n_atom * dim] = { 0 };
double mhessian[n_atom * dim * n_atom * dim] = { 0 };
double Mtherm[n_atom * dim * n_atom * dim] = { 0 };
double M[dim * n_atom] = { 1837.7, 1837.7, 1837.7, 16 * 1837.7, 16 * 1837.7, 16 * 1837.7, 1837.7, 1837.7, 1837.7 };
double x_eq[N] = { 1.6774683, -0.04038825, 0., -0.1297673, 0.06031243, 0, -0.48502213, 1.83515775, 0 };

// 尽量不要动下面定义的变量
double x_array[P][n_atom * dim] = { 0 };     // 储存所有beads位置
double p_array[P][n_atom * dim] = { 0 };     // 储存所有beads动量
double nabla_V[dim * n_atom] = { 0 };  // 储存单个beads在Cartesian坐标下的势能梯度
double nabla_phi[P][dim * n_atom] = { 0 };  // 每个bead在staging坐标下的受力



unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> distribution(0.0, 1.0);
// 设置随机数种子

void print_matrix(double* x, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << x[i * m + j] << " ";
        }
        cout << endl;
    }
}

void print_vector(double* x, int n) {
    for (int k = 0; k < n; k++) {
        cout << x[k] << " ";
    }
    cout << endl;
}

/*
staging变换，将Cartesian坐标变换成为staging坐标

x: double 2-D array 所有beads的Cartesian坐标
没有返回值，直接修改传入的数组
*/
void staging_transf(double(&x)[P][dim * n_atom]) {
    // 第一个bead的坐标不用动
    for (int i = 1; i < P; i++) {
        for (int j = 0; j < dim * n_atom; j++) {
            x[i][j] = x[i][j] - (double)(i * x[(i + 1) % P][j] + x[0][j]) / (i + 1);
        }
    };
}

/*
staging逆变换，将staging坐标转换成为Cartesian坐标
s: double 2-D array 所有beads的staging坐标
没有返回值
*/
void staging_transf_inv(double(&s)[P][dim * n_atom]) {
    // 第一个bead不用动
    for (int j = P - 1; j > 0; j--) {
        for (int i = 0; i < dim * n_atom; i++) {
            s[j][i] = s[j][i] + (double)(j * s[(j + 1) % P][i] + s[0][i]) / (j + 1);
        }
    }
}

extern "C" {
    void h2opot_(double* v, double* dv, double* ddv, double* x, int* flag);
}

void nabla_V_Water(double* x, bool cstr = false) {
    double v = 0;
    int* nptr;
    int n = 1;
    nptr = &n;
    double ddv[9 * 9] = { 0 }; // 一定要初始化
    h2opot_(&v, nabla_V, ddv, x, nptr);
}


void Mhessian_Water(double* x) {
    double v = 0;
    double dv[n_atom * dim] = { 0 };
    int n = 2;
    h2opot_(&v, dv, hessian, x, &n);
    for (int i = 0; i < n_atom * dim; i++) {
        for (int j = 0; j < n_atom * dim; j++) {
            mhessian[i * n_atom * dim + j] = 1.0 /
                sqrt(M[i] * M[j]) * hessian[i * n_atom * dim + j];
        }
    }
}


void M_therm(double(&x)[dim * n_atom], void (*Mhessian)(double(*))) {
    double egv[n_atom * dim] = { 0 };
    int n = dim * n_atom;
    int supp[2 * n_atom * dim] = { 0 };
    double eg_vector[dim * n_atom * dim * n_atom] = { 0 };
    Mhessian_Water(x);
    double mhessian_cp[N * N] = { 0 };
    for (int k = 0; k < N * N; k++) {
        mhessian_cp[k] = mhessian[k];
    }
    LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', dim * n_atom, mhessian_cp,
        n_atom * dim, 0.0, 0.0, 0.0, 0.0, 1e-7, &n, egv, eg_vector,
        n_atom * dim, supp);
    // 特征值按升序排列，特征向量按列排列
    // 会修改传入的对角化的矩阵，应备份
    double eg_vector_cp[N * N] = { 0 };
    for (int k = 0; k < N * N; k++) {
        eg_vector_cp[k] = eg_vector[k];
    }
    // 构建Q(u)
    double Q_u[N] = { 0 };
    double u = beta / 2.0;
    for (int k = 0; k < N; k++) {
        if (abs(egv[k]) < 1e-12) {  // 小心数值误差111
            Q_u[k] = 1;
        }
        else {
            if (egv[k] >= 0) {
                u = u * sqrt(egv[k]);
                Q_u[k] = u / tanh(u);
            }
            else {
                u = u * sqrt(-egv[k]);
                Q_u[k] = tanh(u) / u;
            }
        }
    }
    double M_sup[N * N] = { 0 };
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double var = 0;
            for (int k = 0; k < N; k++) {
                var += eg_vector[i * N + k] * Q_u[k] * eg_vector[j * N + k];
            }
            M_sup[i * N + j] = var;
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Mtherm[i * n + j] = sqrt(M[i] * M[j]) * M_sup[i * N + j];
        }
    }
}

/*
staging坐标下的所有beads的外势
s_chain: double 2-D array 所有beads的staging坐标
return: double 势能
*/
double phi_potential(double(&s_chain)[P][dim * n_atom],
    double (*potential)(double(&)[dim * n_atom])) {
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
    void (*nabla_potential)(double(*), bool)) {
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
            nabla_potential(s_chain[j], false);  // 这一步修改了nabla_V中的值
            // 经过staging_transf_inv之后，s_chain表示对应的Cartesian坐标
            for (int k = 0; k < N; k++) {
                nabla_phi[0][k] += nabla_V[k];
            }
        }
    }
    // 后面的每个bead的势能梯度都需要递推计算
    for (int j = 1; j < P; j++) {
        nabla_potential(s_chain[j], false);  // 计算第j个cartisan坐标下的梯度
        for (int k = 0; k < N; k++) {
            nabla_phi[j][k] =
                (double)nabla_V[k] / P + (double)(j - 1) / j * nabla_phi[j - 1][k];
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
    void (*nabla_potential)(double(&)[N], bool)) {
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
        nabla_potential(s_chain[j], false);
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
    void (*nabla_potential)(double(*), bool),
    double beta,
    double omega_P,
    double c1_PIMD,
    double c2_PIMD,
    double d1_PIMD,
    double d2_PIMD) {
    // 首先演化半步B动量项
    nabla_phi_potential(s_chain, nabla_potential);
    for (int j = 0; j < P; j++) {
        for (int k = 0; k < N; k++) {
            p_chain[j][k] += -0.5 * dt_PIMD * nabla_phi[j][k];
        }
    }

    // 其次演化<半步>A内势力项
    for (int k = 0; k < N; k++) {
        s_chain[0][k] += +0.5 * dt_PIMD / M[k] * p_chain[0][k];
    }
    double var = 0;  // 储存计算过程中的中间变量
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
        p_chain[0][k] =
            c1_PIMD * p_chain[0][k] + c2_PIMD * sqrt(M[k]) * distribution(generator);
    }
    for (int j = 1; j < P; j++) {
        for (int k = 0; k < N; k++) {
            p_chain[j][k] =
                c1_PIMD * p_chain[j][k] +
                c2_PIMD * sqrt((double)(j + 1) / j * M[k]) * distribution(generator);
        }
    }

    // 演化<半步>A内势力项
    for (int k = 0; k < N; k++) {
        s_chain[0][k] += +0.5 * dt_PIMD / M[k] * p_chain[0][k];
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
            p_chain[j][k] += -0.5 * dt_PIMD * nabla_phi[j][k];
        }
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
        T += p_chain[0][k] / M[k] * p_chain[0][k];
    }
    for (int j = 1; j < P; j++) {
        for (int k = 0; k < N; k++) {
            T += p_chain[j][k] * j / M[k] * p_chain[j][k] / (j + 1);
        }
    }
    return T / (N * P * kb);
}

/*
热质量矩阵的estimator，只保留对角元
*/
double M_t[n_atom * dim] = { 0 };
void M_therm_average(void (*Mhessian)(double(*))) {
    staging_transf_inv(x_array);
    for (int k = 0; k < N; k++) {
        M_t[k] = 0;
    }
    for (int j = 0; j < P; j++) {
        M_therm(x_array[j], Mhessian);
        for (int k = 0; k < N; k++) {
            M_t[k] += Mtherm[k * N + k] / P;
        }
    }
    staging_transf(x_array);
}

/*
最终的计算程序，输入势能梯度函数与质量加权hessian矩阵的函数
*/
void calculate(void (*nabla_V_potential)(double(*), bool), void (*Mhessian)(double(*))) {

    time_t t_start, t_end;  // 记录程序运行时间
    t_start = clock();

    // 构造PIMD的初始条件
    for (int j = 0; j < P; j++) {
        for (int k = 0; k < N; k++) {
            x_array[j][k] = x_eq[k];
            p_array[j][k] = 0;
        }
    }
    staging_transf(x_array); // 将Cartesian坐标转化到staging坐标
    double T = 0; // 用于判断PIMD是否达到平衡?
    double M_t_out[N] = { 0 };
    cout << "------------------START--PIMD------------------" << endl;
    cout << "the number of beads is " << P << endl;
    cout << "the setting temperature is " << t_var << endl;
    cout << "total steps is " << steps_PIMD << endl;

    for (int i = 0; i < steps_PIMD; i++) {  // 用PIMD对量子正则分布进行抽样
        BAOAB(x_array, p_array, nabla_V_potential, beta, omega_P, c1_PIMD,
            c2_PIMD, d1_PIMD, d2_PIMD);
        T += temperature(p_array) / steps_PIMD;
        if (i % 10000 == 0) {
            t_end = clock();
            printf("in PIMD, step = %d, time has been used = %.3f s\n", i,
                (double)(t_end - t_start) / CLOCKS_PER_SEC);
        }
        if (i % step_leap == 0) {
            M_therm_average(Mhessian);
            for (int k = 0; k < N; k++) {
                M_t_out[k] += M_t[k] / steps_PIMD * step_leap;
            }
            //print_vector(M_t, N);
        }
    }
    cout << endl;
    cout << "temperature of the system is " << T << ", ";
    cout << "and the setting temperature is " << t_var << endl;
    print_vector(M_t_out, N);
    cout << "-------------------END--PIMD-----------------" << endl;
}

int main() {

    calculate(nabla_V_Water, Mhessian_Water);
    int i;  // 让程序不要运行完就退出
    cin >> i;
    return 0;
}
