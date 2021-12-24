#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <mkl.h>

using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

// 约定所有单位均为原子单位制
const int N = 3;                             // 系统的总自由度
const double dt_MD = 1;                      // Langevin动力学的时间步长
const double dt_HD = 0.2;                    // Hamilton动力学的时间步长
const double T = 200;                        // 设定温度
const double kb = 1.28065e-23 / 4.35974e-18; // Boltzman常数
const double beta = 1.0 / (kb * T);
const double gama = 0.5;     // 摩擦系数 怎么选择比较合适？
const int steps_MD = 5000000; // 分子动力学（采样）的总步长
const int n_traj = 20000;     // Hamilton动力学时选取的轨线数目
const int steps_HD = 80000;    // Hamilton动力学总步长
const int step_leap = 5;      // 计算关联函数时跳过的步数
const bool constrain = true;
const double pi = 3.1415926535897932;
double x[N] = { 0 };
double p[N] = { 0 }; // 储存演化时的中间变量
double x_traj[N * n_traj] = { 0 };
double p_traj[N * n_traj] = { 0 };                         // 储存用于实时间动力学的轨线初始值
double nabla_V[N] = { 0 };                                 // 储存势能梯度
double f[N] = { 0 };                                       // 储存受力
double* x_array; // 用于计算关联函数的数据的列表
double* p_array; // 用mkl的函数申请内存可以避开数组长度的限制
double t_array[steps_HD / step_leap] = { 0 }; // 储存每一帧的时间
double p_corr[steps_HD / step_leap] = { 0 };
double x_corr[steps_HD / step_leap] = { 0 }; // 储存时间关联函数，暂时未使用

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
// 第一个参数为高斯分布的平均值，第二个参数为标准差
std::normal_distribution<double> distribution(0.0, 1.0);

// 定义一些与Morse势有关的常数
/* 参考2011年刘老师的文章
Insights in quantum dynamical effects in the infrared spectroscopy of liquid
water from a semiclassical study with an ab initio-based flexible and
polarizable force field J. Chem. Phys. 135, 244503 (2011);
https://doi.org/10.1063/1.3670960
*/
const double De = 0.1851;    // 势阱深度
const double alpha = 1.2102; // 势阱宽度倒数
const double r_eq = 1.7799;  // 平衡位置
const double m_H = 1837.7;
const double m_O = 16 * m_H;
const double mu = m_H * m_O / (m_H + m_O); // 约化质量
double M[3] = { mu, mu, mu };
double M_inv[3] = { 1.0 / mu, 1.0 / mu, 1.0 / mu };
/*
Cartesian坐标下Morse势
x: doubel 1-D array Cartesian坐标
N: int 系统的自由度，在本例子中为3
return: double 对应坐标下的势能
周期在400左右，MD与HD中时间步长选取可以参考
*/
double Morse_potential(double* x)
{ // 对于指针来说 p[i] = *(p + i)
    double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    double q = (1 - exp(-alpha * (r - r_eq)));
    return De * q * q;
}
/*
Cartesian坐标下真实势能的梯度
x: double 1-D array Cartesian坐标
N: int 系统的自由度
*/
void nabla_Morse_potential(double* x)
{
    double r = 0;
    for (int k = 0; k < N; k++)
    {
        r += x[k] * x[k];
    }
    r = sqrt(r);
    double exp_a = exp(-alpha * (r - r_eq));
    double fractor = 2.0 * De * alpha * (1.0 - exp_a) * exp_a / r;
    for (int j = 0; j < N; j++)
    {
        nabla_V[j] = fractor * x[j];
    }
}

/*
一个三维谐振子势能，用于检查MD与HD程序中的问题
*/
void nabla_Harmonic_potential(double* x)
{
    for (int k = 0; k < N; k++)
    {
        nabla_V[k] = 0.01 * pi * 0.01 * pi * M[k] * x[k];
    }
}

/*
对速度的约束，注意这里是约化的两体问题，不能用对多体问题使用的总角动量为0的约束
这里直接保留径向速度，完全扣除角向速度
*/
void constrain_p(double* x, double* p)
{
    double v[N] = { 0 }; // 速度
    for (int k = 0; k < N; k++)
    {
        v[k] = p[k] / M[k];
    }
    double r2 = 0; // 距离
    for (int k = 0; k < N; k++)
    {
        r2 += x[k] * x[k];
    }
    double v_r = 0;
    for (int k = 0; k < N; k++)
    {
        v_r += v[k] * x[k];
    }
    for (int k = 0; k < N; k++)
    { // 传递引用应该可以修改原有的数组
        p[k] = v_r * x[k] * M[k] / r2;
    }
}

const double c1 = exp(-gama * dt_MD);
const double c2 = sqrt(1.0 - c1 * c1);
/*
以BAOAB的方式演化Langevin方程
pre: 2*N维列表，pre[0]代表位置，pre[1]代表动量
dt: double 时间步长
gama: double Langevin动力学中的阻力常数
return: 2*N维列表，表示位置和动量
*/
void BAOAB(double* x,
    double* p,
    void (*nabla_V_potential)(double*),
    double dt,
    double c1,
    double c2)
{
    nabla_V_potential(x);
    for (int i = 0; i < N; i++) // 半步动量
    {
        p[i] = p[i] - dt / 2.0 * nabla_V[i];
    }
    for (int i = 0; i < N; i++) // 半步位置
    {
        x[i] = x[i] + p[i] / M[i] * dt / 2.0;
    }
    for (int i = 0; i < N; i++) // 控温
    {
        p[i] = c1 * p[i] + c2 * sqrt(M[i] / beta) * distribution(generator);
    }
    for (int i = 0; i < N; i++)
    {
        x[i] = x[i] + p[i] / M[i] * dt / 2.0;
    }
    nabla_V_potential(x);
    for (int i = 0; i < N; i++)
    {
        p[i] = p[i] - dt / 2.0 * nabla_V[i];
    }
}

/*
哈密顿动力学
*/
void Velocity_Verlet(double* x,
    double* p,
    void (*nabla_V_potential)(double*),
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
计算指定构型的温度
*/
double temperature(double* p)
{
    double T = 0;
    for (int k = 0; k < N; k++)
    {
        T += p[k] * p[k] / M[k];
    }
    return T / (N * kb);
}

/*
计算给定轨线集合的位置时间关联函数，未考虑在轨线上的时间平均
*x_list: 为储存所有帧位置的数组的头指针
idx: 需要计算时间关联函数的帧数
*/
double x_corr_pmt(double* x_list, int idx)
{
    double corr = 0;
    for (int j = 0; j < n_traj; j++)
    {
        double var = 0;
        for (int k = 0; k < N; k++)
        {
            var += x_list[N * steps_HD / step_leap * j + k] * x_list[N * steps_HD / step_leap * j + idx * N + k];
        }
        corr += var / n_traj;
    }
    return corr;
}

/*
计算给定轨线集合的动量时间关联函数，未考虑在轨线上的时间平均
*p_list: 为储存所有帧位置的数组的头指针
idx: 需要计算时间关联函数的帧数
*/
double p_corr_pmt(double* p_list, int idx)
{
    double corr = 0;
    for (int j = 0; j < n_traj; j++)
    {
        double var = 0;
        for (int k = 0; k < N; k++)
        {
            var += p_list[N * steps_HD / step_leap * j + k] * p_list[N * steps_HD / step_leap * j + idx * N + k];
        }
        corr += var / n_traj;
    }
    return corr;
}

int main()
{
    x[0] = r_eq;
    x[1] = 0.1;
    x[2] = 0.1;
    cout << "===================="
        << "START LANGEVIN DYNAMICS"
        << "====================" << endl;
    x_array = (double *)MKL_malloc((steps_HD / step_leap) * N * n_traj * sizeof(double), 128);
    p_array = (double *)MKL_malloc((steps_HD / step_leap) * N * n_traj * sizeof(double), 128);
    time_t begin, end;
    begin = clock();
    double tem = 0;
    for (int step = 0; step < steps_MD; step++)
    {
        BAOAB(x, p, nabla_Morse_potential, dt_MD, c1, c2);
        tem += temperature(p);
        if (step % 1000000 == 0)
        {
            end = clock();
            cout << "In steps " << step << ", time = " << (double)(end - begin) / CLOCKS_PER_SEC << ", temperature = " << double(tem / step) << endl;
        }
        if (step >= steps_MD - n_traj)
        {
            for (int k = 0; k < N; k++)
            {
                x_traj[(step + n_traj - steps_MD) * N + k] = x[k];
                p_traj[(step + n_traj - steps_MD) * N + k] = p[k];
            }
        }
    }
    cout << "position is " << x[0] << " " << x[1] << " " << x[2] << endl;
    cout << "====================="
        << "END LANGEVIN DYNAMICS"
        << "=====================" << endl;
    cout << "===================="
        << "START HAMILTON DYNAMICS"
        << "====================" << endl;
    // 构建每一帧的时间
    for (int n = 0; n < steps_HD / step_leap; n++)
    {
        t_array[n] = n * dt_HD * step_leap;
    }

    for (int j = 0; j < n_traj; j++)
    {
        // j表示在第j根轨线，首先修改x，p为初始值
        if (j % 1000 == 0)
        {
            end = clock();
            cout << "in trajectory " << j << ", time = " << (double)(end - begin) / CLOCKS_PER_SEC << endl;
        }
        for (int k = 0; k < N; k++)
        {
            x[k] = x_traj[N * j + k];
            p[k] = p_traj[N * j + k];
        }
        // 此时对初始动量进行约束
        if (constrain)
        {
            constrain_p(x, p);
        }
        for (int step = 0; step < steps_HD; step++)
        {
            // 首先应该在x_array, p_array中储存轨线，按每条轨线的顺序
            // 存储，即[N * steps_HD / step_leap * j]为第j条轨线的头指针
            if (step % step_leap == 0)
            {
                for (int k = 0; k < N; k++)
                {
                    x_array[steps_HD / step_leap * N * j + N * (step / step_leap) + k] = x[k];
                    p_array[steps_HD / step_leap * N * j + N * (step / step_leap) + k] = p[k];
                }
            }
            Velocity_Verlet(x, p, nabla_Morse_potential, dt_HD);
        }
    }
    // 计算时间关联函数
    stringstream fmt1; // 通过文件存储时间关联函数的数据
    fmt1 << "3-D_Morse_p_cstr"
        << ".txt";
    stringstream fmt2; // 通过文件存储时间关联函数的数据
    fmt2 << "3-D_Morse_x_cstr"
        << ".txt";
    ofstream OutFile1(fmt1.str());
    ofstream OutFile2(fmt2.str());
    for (int t = 0; t < steps_HD / step_leap; t++)
    {
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
            << t_array[t] << "  ";
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
            << p_corr_pmt(p_array, t) << "  \n"; // p_array[N * steps_HD / step_leap * 200 +  N * t + 2]
        OutFile2 << std::setiosflags(std::ios::scientific | std::ios::showpos)
            << t_array[t] << "  ";
        OutFile2 << std::setiosflags(std::ios::scientific | std::ios::showpos)
            << x_corr_pmt(x_array, t) << "  \n";
    }
    OutFile1.close();
    fmt1.clear();
    OutFile2.close();
    fmt2.clear();
    cout << "====================="
        << "END HAMILTON DYNAMICS"
        << "=====================" << endl;
    int i;
    cin >> i;
    return 0;
}