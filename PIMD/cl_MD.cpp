#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <windows.h>

//using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

// 关于时间关联函数的一些模拟，暂时局限于一维体系

const int N = 1;                 // 系统的总自由度
const double M = 1;              // 质量矩阵的对角元
const double M_inv = 1;          // 质量矩阵逆的对角元
const double dt = 0.02;          // 时间间隔
const double gamma = 0.5;        // 摩擦系数 怎么选择比较合适？
const int total_steps = 1000000; // 总步长
const double x_0 = 1.5;
double x = 0;
double p = 0;
double nabla_V[N] = {0}; // 储存中间变量

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
// 第一个参数为高斯分布的平均值，第二个参数为标准差
std::normal_distribution<double> distribution(0.0, 1.0);

/*
void print_vector_2(double vec)
{
    cout << "[";
    for (int i = 0; i < vec.size(); i++)
    {
        cout << "[";
        for (int j = 0; j < vec[i].size(); j++)
        {
            cout << vec[i][j] << ",";
        }
        cout << "], ";
    }
    cout << "]" << endl;
}
*/

/*
简谐势能
x: vector<double> Cartesian坐标
return: double
*/
double V_harmonic(double(&x))
{
    return 0.5 * (x - x_0) * (x - x_0);
}

/*
Cartesian坐标下真实势能的梯度
x: vector<double> Cartesian坐标
return: vector<double> 势能的梯度
*/
void nabla_V_potential(double(&x))
{
    nabla_V = x;
}

/*


*/
/*
以BAOAB的方式演化Langevin方程
pre: 2*N维列表，pre[0]代表位置，pre[1]代表动量
dt: double 时间步长
gama: double Langevin动力学中的阻力常数
return: 2*N维列表，表示位置和动量
*/
double BAOAB(double pre, double dt, double gama = 0.5, double beta = 1)
{
    // 首先定义一些量
    double nxt;
    int n = pre[0].size(); // 空间维数
    // 定义了在O过程中出现的两个常数
    double c1 = exp(-gama * dt);
    double c2 = sqrt(1 - c1 * c1);

    vector<double> x = pre[0];
    vector<double> p = pre[1]; // 分别提取位置和动量信息

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

double kinetic_energy(vector<double> p)
{
    double T = 0;
    for (int i = 0; i < p.size(); i++)
    {
        T += p[i] * p[i] / M[i][i] / 2;
    }
    return T;
}

/*
系综平均 
c_list: 所有构型的列表
estimator: 估计量
return: 平均值
*/
double ensemble_average(double(&c_list), double (&estimator)(vector<double>))
{
    int tot = c_list.size();
    double average = 0;
    for (int i = 0; i < tot; i++)
    {
        average += estimator(c_list[i]);
    }
    return average / tot;
}

int main()
{
    double i;
    cin >> i;
    return 0;
}