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
const int steps_PILD = 1000000;  // 演化的总步数 
const int step_leap = 100;  // 跳过step_leap步记录一次位置和动量 计算时间关联函数的最小步长
const double t_var = 300;                        // 设定温度
const double kb = 1.28065e-23 / 4.35974e-18; // Boltzman常数
const double beta = 1.0 / (kb * t_var);  // 温度的倒数
const double dt_PIMD = 0.25; // PIMD中的时间步长
const double dt_PILD = 0.25;
const double omega_P = sqrt(P) / beta;
const double gama_ad = 1e-3;  // 绝热参数
const double omega_ad = sqrt(P / gama_ad) / beta;
const int n_x = 10;  // 坐标空间中通过PIMD采样的点数目
const int n_p = 3;  // 固定x在动量空间中采样的点数目
const int steps_PIMD = 400000; // 坐标空间抽样时使用的PIMD步数
const bool cstr_p = true;  // 是否约束hamilton动力学的初始位置和动量
const bool cstr_f = true;  // 是否在每一步中约束受力
const int cstr_f_step = 100; // 每隔cstr_f_step 步修正一次力

int step_now = 0;
// 这里认为热质量矩阵是对角的
double Mtherm[n_atom * dim] = { 1780.3,1785.01, 1796.11, 29287.2, 29297.3, 29317.3, 1779.2, 1784.46, 1794.71 };
double M[dim * n_atom] = { 1837.7, 1837.7, 1837.7, 16 * 1837.7, 16 * 1837.7, 16 * 1837.7, 1837.7, 1837.7, 1837.7 };
double x_eq[N] = { 1.6774683, -0.04038825, 0., -0.1297673, 0.06031243, 0, -0.48502213, 1.83515775, 0 };
double hessian[N * N] = { 0 };
double mhessian[N * N] = { 0 };

// 尽量不要动下面定义的变量
const int t_steps_max = steps_PILD / step_leap - 100;
double IDE[dim * dim] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
double x_array[P][n_atom * dim] = { 0 };     // 储存所有beads位置
double p_array[P][n_atom * dim] = { 0 };     // 储存所有beads动量
double nabla_V[dim * n_atom] = { 0 };  // 储存单个beads在Cartesian坐标下的势能梯度
double nabla_phi[P][dim * n_atom] = { 0 };  // 每个bead在staging坐标下的受力
double traj_x_init[n_p * n_x][P][n_atom * dim] = { 0 };  // 储存选定为初始值的平衡构型
double traj_p_init[n_p * n_x][P][n_atom * dim] = { 0 };  // 用作动力学的初始值
double x_list_traj[n_p * n_x][steps_PILD / step_leap][n_atom * dim] = { 0 };  // 所有轨线不同时刻的位置
double p_list_traj[n_p * n_x][steps_PILD / step_leap][n_atom * dim] = { 0 };  // 所有轨线不同时刻的动量
double t_list[steps_PILD / step_leap] = { 0 };  // 记录每个构型对应的演化时间


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

/*
对位置的约束
将水分子放在质心系
*/
void constrain_x(double* x) {
    double x_c[dim] = { 0 }; // 质心坐标
    double m_c = 0;  // 质心质量
    for (int n = 0; n < n_atom; n++) {
        m_c += M[n * dim];
    }
    for (int k = 0; k < dim; k++) {
        double var = 0;
        for (int n = 0; n < n_atom; n++) {
            var += x[dim * n + k] * M[dim * n + k];
        }
        x_c[k] = var / m_c;
    }
    for (int k = 0; k < dim; k++) {
        for (int n = 0; n < n_atom; n++) {
            x[dim * n + k] = x[dim * n + k] - x_c[k];
        }
    }
}

/*
计算惯量张量
计算结果储存在double指针I指向的位置
*/
void I_tensor(double* x, double* I) {
    double r_2[n_atom] = { 0 }; // 储存原子距原点的距离
    for (int n = 0; n < n_atom; n++) {
        double var = 0;
        for (int k = 0;k < dim; k++) {
            var += x[n * dim + k] * x[n * dim + k];
        }
        r_2[n] = var;
    }
    int d_ij = 0;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            double var = 0;
            if (i == j) {
                d_ij = 1;
            }
            else {
                d_ij = 0;
            }
            for (int n = 0; n < n_atom; n++) {
                var += M[n * dim] * (r_2[n] * d_ij - x[n * dim + i] * x[n * dim + j]);
            }
            I[i * dim + j] = var;
        }
    }
}

/*
计算惯量张量的逆
*/
void I_tensor_inv(double* I, double* I_inv) {
    for (int i = 0; i < dim * dim; i++) {
        I_inv[i] = IDE[i];  // 首先拷贝一份单位矩阵
    }
    double I_copy[dim * dim] = { 0 };
    for (int j = 0; j < dim * dim; j++) {
        I_copy[j] = I[j];
    }
    int ipiv[dim] = { 0 };  // 不知道lapack拿他干什么，但是必须有
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, dim, dim, I_copy, dim, ipiv, I_inv, dim);
    // 上述函数解线性方程组，当A*X=B中B为单位矩阵时解得X为A的逆
    // LAPACK_ROW_MAJOR表示行主序(C++)，fortran中数组为列主序
    // 上述函数会用解得的X覆盖B矩阵，如果需要使用B，应拷贝
    // 经过测试A矩阵也会被覆盖
}

/*
角动量
修改传入double指针L指向的内存
*/
void angular_momentum(double* x, double* v, double* L) {
    for (int k = 0; k < dim; k++) {
        double var = 0;
        for (int n = 0; n < n_atom; n++) {
            var += M[n * dim] * (x[n * dim + (k + 1) % dim]
                * v[n * dim + (k + 2) % dim] - v[n * dim + (k + 1) % dim]
                * x[n * dim + (k + 2) % dim]);
        }
        L[k] = var;
    }
}

/*
对速度的约束
经过反复测试，没有问题
测试方法，任意给定位置和动量，计算约束后的角动量
*/
void constrain_p(double* x, double* p)
{
    // 将位置修正到质心系
    constrain_x(x);
    double v[dim * n_atom] = { 0 }; // 构建每个原子的速度
    for (int i = 0; i < n_atom * dim; i++) {
        v[i] = p[i] / M[i];
    }
    // 将上述速度修正到质心参考系
    constrain_x(v);
    double L[dim] = { 0 }; // 体系的角动量
    angular_momentum(x, v, L);
    double I[dim * dim] = { 0 }; // 体系的惯量张量 用一维数组存储，对称，无论行或列优先
    I_tensor(x, I);
    double I_inv[dim * dim] = { 0 }; // 惯量张量的逆
    I_tensor_inv(I, I_inv);
    double omega[dim] = { 0 }; // 角速度
    cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, I_inv, dim, L, 1, 0.0, omega, 1);
    // CblasRowMajor: 行主序; CblasNoTrans: 不转置
    // 最后的计算结果储存在omega指针对应的地址
    double v_corr[dim * n_atom] = { 0 };  // 角速度与位置的叉乘
    for (int n = 0; n < n_atom; n++) {
        for (int k = 0; k < dim; k++) {
            v_corr[dim * n + k] = omega[(k + 1) % dim] * x[n * dim + (k + 2) % dim]
                - omega[(k + 2) % dim] * x[n * dim + (k + 1) % dim];
        }
    }
    // 最后将质心系中的速度减去v_corr, 计算动量
    for (int j = 0; j < n_atom * dim; j++) {
        p[j] = (v[j] - v_corr[j]) * M[j];
    }
}

/*
计算力矩
*/
void torque(double* x, double* f, double* T) {
    for (int k = 0; k < dim; k++) {
        double var = 0;
        for (int n = 0; n < n_atom; n++) {
            var += x[n * dim + (k + 1) % dim] * f[n * dim + (k + 2) % dim]
                - f[n * dim + (k + 1) % dim] * x[n * dim + (k + 2) % dim];
        }
        T[k] = var;
    }
}

/*
对力的约束，注意传入的是势能梯度不是力
要在函数内部将输入转换为力再约束最后转换为梯度
*/
void constrain_force(double* x, double* f) {
    double force[n_atom * dim] = { 0 };
    for (int k = 0; k < n_atom * dim; k++) {
        force[k] = -f[k];
    }
    // 首先定义质心受力
    double f_c[dim] = { 0 };
    for (int k = 0; k < dim; k++) {
        double var = 0;
        for (int n = 0; n < n_atom; n++) {
            var += force[dim * n + k];
        }
        f_c[k] = var / n_atom;
    }
    // 将质心受力扣除
    for (int k = 0; k < dim; k++) {
        for (int n = 0; n < n_atom; n++) {
            force[n * dim + k] = force[n * dim + k] - f_c[k];
        }
    }
    // 将坐标转换为质心系
    constrain_x(x);
    // 计算此时的合外力矩
    double N_tor[dim] = { 0 };
    torque(x, force, N_tor);
    // 计算角加速度
    double alpha[dim] = { 0 };
    //构造惯量张量和逆
    double I[dim * dim] = { 0 };
    double I_inv[dim * dim] = { 0 };
    I_tensor(x, I);
    I_tensor_inv(I, I_inv);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, I_inv, dim, N_tor, 1, 0.0, alpha, 1);
    // 计算每个原子受力的修正
    double f_corr[dim * n_atom] = { 0 };
    for (int k = 0; k < dim; k++) {
        for (int n = 0; n < n_atom; n++) {
            f_corr[n * dim + k] = M[n * dim] * (alpha[(k + 1) % dim] * x[n * dim + (k + 2) % dim]
                - alpha[(k + 2) % dim] * x[n * dim + (k + 1) % dim]);
        }
    }
    // 修正每个原子的受力
    for (int j = 0; j < n_atom * dim; j++) {
        f[j] = -force[j] + f_corr[j];
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
        Mtherm[i] = sqrt(M[i] * M[i]) * M_sup[i * N + i];
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
                nabla_phi[0][k] += nabla_V[k] / P;
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
    void (*nabla_potential)(double(*), bool),
    bool cstr,
    double beta,
    double omega_ad,
    double c1_PILD,
    double c2_PILD,
    double d1_PILD,
    double d2_PILD) {
    double var = 0;  // 由于储存计算过程中的中间变量

    // 首先演化半步B动量项 第一个bead要和其他bead分开演化
    nabla_phi_potential(s_chain, nabla_potential);
    if (cstr) {
        constrain_force(s_chain[0], nabla_phi[0]);
    }
    for (int k = 0; k < N; k++) {
        p_chain[0][k] += -0.5 * dt_PILD * Mtherm[k] / M[k] * nabla_phi[0][k];
        // M_therm(s_chain[0], Mhessian)
    }
    for (int j = 1; j < P; j++) {
        for (int k = 0; k < N; k++) {
            p_chain[j][k] += -0.5 * dt_PILD * nabla_phi[j][k];
        }
    }

    // 其次演化<半步>A内势力项
    for (int k = 0; k < N; k++) {
        s_chain[0][k] += +0.5 * dt_PILD / M[k] * p_chain[0][k];
    }
    for (int j = 1; j < P; j++) {
        for (int k = 0; k < N; k++) {
            var = s_chain[j][k] * d1_PILD +
                d2_PILD * j / ((j + 1) * M[k] * omega_ad * gama_ad) * p_chain[j][k];
            p_chain[j][k] = -omega_ad * d2_PILD * (j + 1) / j * M[k] * gama_ad * s_chain[j][k] +
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
                c2_PILD * sqrt((double)(j + 1) / j * M[k] * gama_ad) * distribution(generator);
        }
    }
    // 其次演化<半步>A内势力项
    for (int k = 0; k < N; k++) {
        s_chain[0][k] += +0.5 * dt_PILD / M[k] * p_chain[0][k];
    }
    for (int j = 1; j < P; j++) {
        for (int k = 0; k < N; k++) {
            var = s_chain[j][k] * d1_PILD +
                d2_PILD * j / ((j + 1) * M[k] * omega_ad * gama_ad) * p_chain[j][k];
            p_chain[j][k] = -omega_ad * d2_PILD * (j + 1) / j * M[k] * gama_ad * s_chain[j][k] +
                d1_PILD * p_chain[j][k];
            s_chain[j][k] = var;
        }
    }

    // 最后演化半步B动量项
    nabla_phi_potential(s_chain, nabla_potential);
    if (cstr) {
        constrain_force(s_chain[0], nabla_phi[0]);
    }
    for (int k = 0; k < N; k++) {
        p_chain[0][k] += -0.5 * dt_PILD * Mtherm[k] / M[k] * nabla_phi[0][k];
        // M_therm(s_chain[0], Mhessian)
    }
    for (int j = 1; j < P; j++) {
        for (int k = 0; k < N; k++) {
            p_chain[j][k] += -0.5 * dt_PILD * nabla_phi[j][k];
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
kubo形式的动量关联函数，只计算了两点的关联函数，没有进行时间平均
start: int
end: int
*/
double p_kubo_corr_func_(int start, int end) {
    double ave = 0;
    for (int j = 0; j < n_x * n_p; j++) {  // 对所有轨线平均
        M_therm(x_list_traj[j][start], Mhessian_Water);
        // 此时不再需要Mtherm来演化轨线，因此可以改动
        for (int k = 0; k < n_atom * dim; k++) {
            ave += M[k] / Mtherm[k] * p_list_traj[j][start][k] *
                p_list_traj[j][end][k] / (n_x * n_p);
        }
    }

    return ave;
}

/*
kubo形式动量关联函数，考虑了在同一轨线不同时间的时间平均
*/
double p_kubo_corr_func(int t, int cnt = 50) {
    double ave = 0;
    for (int j = 0; j < cnt; j++) {
        ave += p_kubo_corr_func_(j, j + t);
    }
    return ave / cnt;
}

/*
最终的计算程序，输入势能梯度函数与质量加权hessian矩阵的函数
*/
void calculate(void (*nabla_V_potential)(double(*), bool)) {

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
        // 
        if (i >= steps_PIMD - n_x) {  // 首先在位形空间采样
            for (int k = 0; k < N; k++) {
                for (int j = 0; j < P; j++) {
                    for (int np = 0; np < n_p; np++) {
                        traj_x_init[(i - steps_PIMD + n_x) * n_p + np][j][k] = x_array[j][k];
                    }
                }
            }
            // 还需要对动量进行采样
            for (int k = 0; k < N; k++) {  // 空间维数
                for (int j = 0; j < P; j++) {
                    if (j == 0) {
                        for (int np = 0; np < n_p; np++) {  // 对于每一个位置，动量空间的采样
                            M_therm(x_array[0], Mhessian_Water);
                            traj_p_init[np + (i - steps_PIMD + n_x) * n_p][j][k] = sqrt(Mtherm[k] / beta)
                                * distribution(generator);
                            // 希望这样采样可可以满足位置和动量的联合分布
                        }
                    }
                    else {
                        for (int np = 0; np < n_p; np++) {  // 对于每一个位置，动量空间的采样
                            traj_p_init[np + (i - steps_PIMD + n_x) * n_p][j][k] = 0;
                        }
                    }
                }
            }
        }
    }
    cout << endl;
    cout << "temperature of the system is " << T << ", ";
    cout << "and the setting temperature is " << t_var << endl;
    double Mtherm[n_atom * dim] = { 1948.34, 1932.72, 1896.91 ,29601.5, 29518.4, 29505.7, 1951.93, 1887.05, 1897.97 };
    cout << "-------------------END--PIMD-----------------" << endl;
    cout << "-----------BEIGN--REAL--TIME--DYNAMICS-------" << endl;
    double p_init_array[n_x * n_p][N] = { 0 };
    //generate_p_init(p_init_array);
    for (int k = 0; k < n_x * n_p; k++) {  // 演化每一根轨线
      // k: int 表示当前演化的轨线
        step_now = 0;
        for (int j = 0; j < P; j++) {
            for (int n = 0; n < N; n++) {  // 初始化位置和动量
                x_array[j][n] = traj_x_init[k][j][n];
                p_array[j][n] = traj_p_init[k][j][n];
            }
        }
        if (cstr_p) {
            constrain_p(x_array[0], p_array[0]);
        }
        for (int j = 0; j < steps_PILD; j++) {  // 演化每一根轨线
            if (j % 2000 == 0) {  // 输出当前轨线信息
                t_end = clock();
                printf(
                    "in PILD, in %dth trajactory, steps = %d,time has been used = %.1f "
                    "s\n",
                    k, j, (double)(t_end - t_start) / CLOCKS_PER_SEC);
            }
            if (j % step_leap == 0) {  // 每隔step_leap步储存一个构型
                for (int r = 0; r < N; r++) {
                    x_list_traj[k][j / step_leap][r] = x_array[0][r];
                    p_list_traj[k][j / step_leap][r] = p_array[0][r];
                    t_list[j / step_leap] = j * dt_PILD;
                }
            }
            BAOAB_PILD_average(x_array, p_array, nabla_V_potential, cstr_f, beta,
                omega_ad, c1_PILD, c2_PILD, d1_PILD, d2_PILD);
            if (step_now % cstr_f_step == 0) {
                constrain_p(x_array[0], p_array[0]);
            }
            step_now += 1;
        }
    }
    stringstream fmt1;  // 通过文件存储时间关联函数的数据
    fmt1 << "water_300_PILD_cstr_p"
        << ".txt";
    ofstream OutFile1(fmt1.str());
    printf("need some time to calculate correlation function\n");
    for (int t = 0; t < t_steps_max; t++) { // 计算时间关联函数
        if (t % 500 == 0) {
            t_end = clock();
            cout << "steps = " << t << ", time = " << (double)(t_end - t_start) / CLOCKS_PER_SEC << endl;
        }
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
            << t_list[t] << "  ";
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
            << p_kubo_corr_func(t, 50) << "  \n";

    }
    OutFile1.close();
    fmt1.clear();
    cout << "-------------------END--PILD-----------------" << endl;

}

int main() {
    calculate(nabla_V_Water);
    int i; // 让程序不要运行完就退出
    cin >> i;
    return 0;
}
