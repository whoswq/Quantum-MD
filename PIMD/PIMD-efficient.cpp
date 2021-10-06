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
#include <windows.h>

//using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

const int N = 1;                 // 系统的自由度
const int P = 128;                 // beads数量
const int steps = 1000000;       // 演化的总步数
const double dt = 0.1 / sqrt(P); // 时间步长 对于BAOAB算法，时间步长最长为2/omega_P
// double h_bar = 1; // 约化Planck常数
// M = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; // 质量矩阵
const double m = 1; // 质量矩阵的对角元
double M[N];        // 质量矩阵 1-D array
double M_inv[N];    // 质量矩阵的逆 1-D array

double x_array[P][N];    // 储存所有beads位置
double p_array[P][N];    // 储存所有beads动量
double nabla_V[N] = {0}; // 储存单个beads在Cartesian坐标下的势能梯度
double nabla_phi[P][N] = {0};

// 设置随机数种子
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> distribution(0.0, 1.0);

void print_vector_2(double (&vec)[P][N])
{
    cout << "[";
    for (int i = 0; i < P; i++)
    {
        cout << "[";
        for (int j = 0; j < N; j++)
        {
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
void staging_transf(double (&x)[P][N])
{
    // 第一个bead的坐标不用动
    for (int i = 1; i < P; i++)
    {
        for (int j = 0; j < N; j++)
        {
            x[i][j] = x[i][j] - (double)(i * x[(i + 1) % P][j] + x[0][j]) / (i + 1);
        }
    };
}

/*
staging逆变换，将staging坐标转换成为Cartesian坐标
s: double 2-D array 所有beads的staging坐标
没有返回值
*/
void staging_transf_inv(double (&s)[P][N])
{
    // 第一个bead不用动
    for (int j = P - 1; j > 0; j--)
    {
        for (int i = 0; i < N; i++)
        {
            s[j][i] = s[j][i] + (double)(j * s[(j + 1) % P][i] + s[0][i]) / (j + 1);
        }
    }
}

/*
Cartesian坐标下真实的势能函数
x: doubel 1-D array Cartesian坐标
return: double 对应坐标下的势能
*/
double V_potential(double (&x)[N])
{
    double v = 0;
    for (int j = 0; j < N; j++)
    {
        v += 0.25 * x[j] * x[j] * x[j] * x[j];
    }
    return v;
}

/*
staging坐标下的所有beads的外势
s_chain: double 2-D array 所有beads的staging坐标
return: double 势能
*/
double phi_potential(double (&s_chain)[P][N])
{
    staging_transf_inv(s_chain);
    double phi = 0;
    for (int j = 0; j < P; j++)
    {
        phi += V_potential(s_chain[j]);
    }
    staging_transf(s_chain);
    return (double)phi / P;
}

/*
Cartesian坐标下真实势能的梯度
x: double 1-D array Cartesian坐标
return: double 1-D array 势能的梯度
*/
void nabla_V_potential(double (&x)[N])
{

    for (int j = 0; j < N; j++)
    {
        nabla_V[j] = x[j] * x[j] * x[j];
    }
}

/*
staging坐标下所有beads的势能梯度，在求解Langevin方程中直接使用
s_chain: double 2-D array 所有beads的staging坐标
return: douible 2-D array 所有beads在staging坐标下的势能梯度
注意这里没有修改传入的数组
*/
void nabla_Phi_potential(double (&s_chain)[P][N])
{
    // 首先将数组初始化 不清楚这一步对效率影响有多大？
    for (int j = 0; j < P; j++)
    {
        for (int k = 0; k < N; k++)
        {
            nabla_phi[j][k] = 0.0;
        }
    }
    staging_transf_inv(s_chain);
    // 将staging坐标变换为Cartesian坐标，方便计算势能
    {
        for (int j = 0; j < P; j++)
        {
            nabla_V_potential(s_chain[j]);
            for (int k = 0; k < N; k++)
            {
                // 经过staging_transf_inv之后，s_chain表示对应的Cartesian坐标
                nabla_phi[0][k] += nabla_V[k] / P;
            }
        }
    }
    // 后面的每个bead的势能梯度都需要递推计算
    for (int j = 1; j < P; j++)
    {
        nabla_V_potential(s_chain[j]);
        for (int k = 0; k < N; k++)
        {
            nabla_phi[j][k] = nabla_V[k] / P + (double)(j - 1) / j * nabla_phi[j - 1][k];
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
double kinetic_energy_viral(double (&s_chain)[P][N], double beta)
{
    staging_transf_inv(s_chain);
    double x_c[N] = {0}; // 构造所有beads质心的坐标
    for (int k = 0; k < N; k++)
    {
        for (int j = 0; j < P; j++)
        {
            x_c[k] += s_chain[j][k];
        }
        x_c[k] = (double)x_c[k] / P;
    }
    double T = 0;
    for (int j = 0; j < P; j++)
    {
        nabla_V_potential(s_chain[j]);
        for (int k = 0; k < N; k++)
        {
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
void BAOAB(double (&s_chain)[P][N], double (&p_chain)[P][N], double beta, double omega_P, double c1, double c2, double d1, double d2)
{
    double var = 0; // 由于储存计算过程中的中间变量
    // 首先演化半步B动量项
    nabla_Phi_potential(s_chain);
    for (int j = 0; j < P; j++)
    {
        for (int k = 0; k < N; k++)
        {
            p_chain[j][k] = p_chain[j][k] - 0.5 * dt * nabla_phi[j][k];
        }
    }

    // 其次演化<半步>A内势力项
    for (int k = 0; k < N; k++)
    {
        s_chain[0][k] = s_chain[0][k] + 0.5 * dt / M[k] * p_chain[0][k];
    }
    for (int j = 1; j < P; j++)
    {
        for (int k = 0; k < N; k++)
        {
            var = s_chain[j][k] * d1 + d2 * j / ((j + 1) * M[k] * omega_P) * p_chain[j][k];
            p_chain[j][k] = -omega_P * d2 * (j + 1) / j * M[k] * s_chain[j][k] + d1 * p_chain[j][k];
            s_chain[j][k] = var;
        }
    }

    // 然后演化一步O控温项
    for (int k = 0; k < N; k++)
    {
        p_chain[0][k] = c1 * p_chain[0][k] + c2 * sqrt(M[k]) * distribution(generator);
    }
    for (int j = 1; j < P; j++)
    {
        for (int k = 0; k < N; k++)
        {
            p_chain[j][k] = c1 * p_chain[j][k] + c2 * sqrt((double)(j + 1) / j * M[k]) * distribution(generator);
        }
    }

    // 再演化半步A内势力项
    for (int k = 0; k < N; k++)
    {
        s_chain[0][k] = s_chain[0][k] + 0.5 * dt / M[k] * p_chain[0][k];
    }
    for (int j = 1; j < P; j++)
    {
        for (int k = 0; k < N; k++)
        {
            var = s_chain[j][k] * d1 + d2 * j / ((j + 1) * M[k] * omega_P) * p_chain[j][k];
            p_chain[j][k] = -omega_P * d2 * (j + 1) / j * M[k] * s_chain[j][k] + d1 * p_chain[j][k];
            s_chain[j][k] = var;
        }
    }

    // 最后演化半步B动量项
    nabla_Phi_potential(s_chain);
    for (int j = 0; j < P; j++)
    {
        for (int k = 0; k < N; k++)
        {
            p_chain[j][k] = p_chain[j][k] - 0.5 * dt * nabla_phi[j][k];
        }
    }
}

/*
计算给定构型的温度
p_chain: double 2-D array 所有beads的动量
return: double 温度
*/
double temperature(double (&p_chain)[P][N])
{
    double T = 0;
    for (int k = 0; k < N; k++)
    {
        T += p_chain[0][k] * M_inv[k] * p_chain[0][k];
    }
    for (int j = 1; j < P; j++)
    {
        for (int k = 0; k < N; k++)
        {
            T += p_chain[j][k] * j * M_inv[k] * p_chain[j][k] / (j + 1);
        }
    }
    return T / (N * P);
}

int main()
{
    for (int i = 0; i < N; i++)
    {
        M[i] = m;
        M_inv[i] = 1 / m;
    }
    double start = GetTickCount();
    // 初始化位置和动量，随机取点
    for (int j = 0; j < P; j++)
    {
        for (int k = 0; k < N; k++)
        {
            x_array[j][k] = distribution(generator);
            p_array[j][k] = distribution(generator);
        }
    }
    // 用于储存数据
    stringstream fmt1;
    fmt1 << "quartic_P_128"
         << ".txt";
    ofstream OutFile1(fmt1.str());
    // 用于计算平均值
    double K;
    double T;
    double V;
    double temp = 0.1;
    double beta;    // 温度的倒数
    double omega_P; // 内势力角频率
    double c1;
    double c2;
    double d1;
    double d2;
    cout << "------------------START------------------" << endl;
    for (int n = 0; n < 50; n++)
    {
        K = 0;
        T = 0;
        V = 0;
        beta = 1 / (0.1 + 0.2 * n);
        omega_P = sqrt(P) / beta;
        // 定义一些BAOAB函数常用的常数
        c1 = exp(-omega_P * dt);
        c2 = sqrt((1 - c1 * c1) / beta);
        d1 = cos(omega_P * dt / 2);
        d2 = sin(omega_P * dt / 2);
        cout << "the number of beads is " << P << endl;
        cout << "the setting temperature is " << (double)1 / beta << endl;
        cout << "total steps is " << steps << endl;
        for (int i = 0; i < steps; i++)
        {
            BAOAB(x_array, p_array, beta, omega_P, c1, c2, d1, d2);
            T += temperature(p_array) / steps;
            K += kinetic_energy_viral(x_array, beta) / steps;
            V += phi_potential(x_array) / steps;
            if (i % 100000 == 0)
            {
                cout << "step" << i << ", time has been used is " << (GetTickCount() - start) / 1000 << " s" << endl;
            }
        }
        cout << "temperature of the system is " << T << ", ";
        cout << "and the setting temperature is " << (double)1 / beta << endl;
        cout << "kinetic energy is " << K << endl;
        cout << "potential energy is " << V << endl;
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos) << T << "  ";
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos) << (double)1 / beta << "  ";
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos) << K << "  ";
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos) << V << "  ";
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos) << K + V << endl;
    }
    OutFile1.close();
    fmt1.clear();
    cout << "-------------------END-------------------" << endl;
    int i;
    cin >> i;
    return 0;
}