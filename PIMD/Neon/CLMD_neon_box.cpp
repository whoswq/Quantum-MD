#include <chrono>
#include <cmath>
#include <fstream>
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

const int N_mol = 108;       // 系统的分子数
const int N = N_mol * 3;     // 系统的总自由度
const double M = 1;          // 质量矩阵的对角元
const double M_inv = 1;      // 质量矩阵逆的对角元
const double dt = 0.005;     // 时间步长
const double L = 5.4155;     // 模拟所用的盒子的边长（正方形盒子）
const double beta = 1.1173;  // 测试选取的温度
const double r_m = 4;        // 势能函数的截断半径
double x_array[N] = {0};     // 储存某时刻下所有分子的位置
double p_array[N] = {0};     // 储存某时刻下所有分子的动量
double force_ij[3] = {0};    // 储存某两个分子之间的势能梯度
double force_i[3] = {0};     // 储存单个分子收到的力
double force_array[N] = {0}; // 储存所有分子的受力情况 都是一维数组，注意处理细节
double x_single[3] = {0};    // 储存单个分子的坐标
double y_single[3] = {0};    // 储存单个分子的坐标（中间变量）

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
// 第一个参数为高斯分布的平均值，第二个参数为标准差
std::normal_distribution<double> distribution(0.0, 1.0);

void print_vector_1(double vec)
{
    cout << "[";
    for (int i = 0; i < vec.size(); i++)
    {
        cout << vec[i] << " ,";
    }
    cout << "]" << endl;
}

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
Cartesian坐标下两个分子之间的势能
x: double 1-D array
y: double 2-D array 两个分子之间的坐标
return: double 两个分子之间产生的势能
*/
double pair_potential(double (&x)[3], double (&y)[3])
{
    double r_2 = (x[0] - y[0]) * (x[0] - y[0]) +
                 (x[1] - y[1]) * (x[1] - y[1]) +
                 (x[2] - y[2]) * (x[2] - y[2]); // 计算两个分子之间距离的平方
    double r_6 = r_2 * r_2 * r_2;               // 这样可以少算乘法
    return (double)4 * ((double)1 / (r_6 * r_6) - (double)1 / r_6);
}

/*
Cartesian坐标下两个分子之间的作用力
x: double 1-D array
y: double 1-D array 位于x的分子受到y分子的作用力
不返回任何值，修改force_ij中的变量
*/
void force_pair(double (&x)[3], double (&y)[3])
{
    double r_2 = (x[0] - y[0]) * (x[0] - y[0]) +
                 (x[1] - y[1]) * (x[1] - y[1]) +
                 (x[2] - y[2]) * (x[2] - y[2]); // 计算两个分子之间距离的平方
    double r_6 = r_2 * r_2 * r_2;
    double r_8 = r_6 * r_2;
    for (int k = 0; k < 3; k++)
    {
        force_ij[k] = (double)24 * (x[k] - y[k]) * (2(double)1 / (r_8 * r_6) - (double)1 / r_8);
    }
}

/*
计算一个分子感受到的势能
j: int 第j个分子
x_array: double 1-D array N 所有分子的当前构型
return: double 该分子感受到的势能
*/
double single_potential(int j, double (&x_array)[N])
{
    // 首先遍历所有的分子所有的自由度
    for (int k = 0; k < N_mol; k++)
    {
        pass;
    }
}

/*
计算所有分子感受到的势能
x_array: double 1-D N 所有分子当前的构型
return: double 所有分子的总势能
*/
double total_potential(double (&x_array)[N])
{
    double potential = 0;
    for (int j = 0; j < N_mol; j++)
    {
    }
}

/*
计算一个分子受到的力
j: int 第j个分子
x_array: double 1-D array N　所有分子的当前构型
*/
void　force_single(int j; double (&x_array)[N])
{
    // 首先对force_i进行初始化 这样的代码段会不会影响效率？
    for (int k = 0; k < 3; k++)
    {
        force_i[k] = 0;
    }
    // 计算受力
}

/*
生成所有分子受力形成的列表，直接用于下一步langevin方程的演化
x_array: double 1-D array N 所有分子的当前构型
没有返回值，直接修改force_array这个数组
*/
void force(double (&x_array)[N])
{
}

/*
以BAOAB的方式演化Langevin方程
需要好好修改！！
*/
double BAOAB(double pre, double dt, double gama = 0.5, double beta = 1)
{
    // 首先定义一些量
    double nxt;
    int n = pre[0].size(); // 空间维数
    // 定义了在O过程中出现的两个常数
    double c1 = exp(-gama * dt);
    double c2 = sqrt(1 - c1 * c1);

    double x = pre[0];
    double p = pre[1]; // 分别提取位置和动量信息

    // 首先演化半步B动量项
    for (int i = 0; i < n; i++)
    {
        p[i] = p[i] - dt / 2 * nabla_V_potential(x)[i];
    }

    // 其次演化<半步>A内势力项
    for (int i = 0; i < n; i++)
    {
        x[i] = x[i] + 1 / M[i][i] * p[i] * dt / 2;
    }

    // 然后演化一步O控温项
    for (int i = 0; i < n; i++)
    {
        p[i] = c1 * p[i] + c2 * sqrt((double)1 / beta) * sqrt(M[i][i]) * distribution(generator);
    }

    // 再演化半步A位置
    for (int i = 0; i < n; i++)
    {
        x[i] = x[i] + 1 / M[i][i] * p[i] * dt / 2;
    }

    // 最后演化半步B动量项
    for (int i = 0; i < n; i++)
    {
        p[i] = p[i] - dt / 2 * nabla_V_potential(x)[i];
    }
    nxt.push_back(x);
    nxt.push_back(p);
    return nxt;
}

double kinetic_energy(double p)
{
    double T = 0;
    for (int i = 0; i < p.size(); i++)
    {
        T += p[i] * p[i] / M[i][i] / 2;
    }
    return T;
}

int main()
{
    double i;
    cin >> i;
    return 0;
}