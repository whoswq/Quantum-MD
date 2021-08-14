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
using std::vector;

// 为了找出PIMD程序的问题，写一个md程序

int N = 1; // 系统的自由度
//vector<vector<double>> M = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; // 质量矩阵
vector<vector<double>> M = {{1}};
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
// 第一个参数为高斯分布的平均值，第二个参数为标准差
std::normal_distribution<double> distribution(0.0, 1.0);

void print_vector_1(vector<double> vec)
{
    cout << "[";
    for (int i = 0; i < vec.size(); i++)
    {
        cout << vec[i] << " ,";
    }
    cout << "]" << endl;
}

void print_vector_2(vector<vector<double>> vec)
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

/*
Cartesian坐标下真实的势能函数
x: vector<double> Cartesian坐标
return: double
*/
double V_potential(vector<double> x)
{
    double x_2 = 0;
    for (int j = 0; j < x.size(); j++)
    {
        x_2 += x[j] * x[j];
    }
    return 0.5 * x_2;
}

/*
Cartesian坐标下真实势能的梯度
x: vector<double> Cartesian坐标
return: vector<double> 势能的梯度
*/
vector<double> nabla_V_potential(vector<double> x)
{
    vector<double> f;

    for (int j = 0; j < x.size(); j++)
    {
        f.push_back(x[j]);
    }
    return f;
}

/*
以BAOAB的方式演化Langevin方程
pre: 2*N维列表，pre[0]代表位置，pre[1]代表动量
dt: double 时间步长
gama: double Langevin动力学中的阻力常数
return: 2*N维列表，表示位置和动量
*/
vector<vector<double>> BAOAB(vector<vector<double>> pre, double dt, double gama = 0.5, double beta = 1)
{
    // 首先定义一些量
    vector<vector<double>> nxt;
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
double ensemble_average(vector<vector<double>>(&c_list), double (&estimator)(vector<double>))
{
    int tot = c_list.size();
    double average = 0;
    for (int i = 0; i < tot; i++)
    {
        average += estimator(c_list[i]);
    }
    return average / tot;
}

/*
平均动能
不明原因导致ensemble_average()运行很慢，重写对应的函数
*/
/*
double average_kinetic_energy(vector<vector<double>> c_list){
    
}
*/

int main()
{
    cout << "----------the program begins----------" << endl;
    vector<double> x_init = {1};
    vector<double> p_init = {0};
    vector<vector<double>> x_p = {x_init, p_init};
    vector<vector<double>> x_list;
    vector<vector<double>> p_list;
    cout << "the initial condition is: " << endl;
    print_vector_2(x_p);
    /*  将轨线储存为文件
    stringstream fmt1;
    fmt1 << "cl_MD_test"
         << ".txt";
    ofstream OutFile1(fmt1.str());
    */
    int total_steps; // 演化的步数
    cout << "total steps = ?" << endl;
    cin >> total_steps;
    double beta; // 系统设定的温度
    cout << "beta = ?" << endl;
    cin >> beta;
    double gama; // 系统的阻力系数
    cout << "gamma = ?" << endl;
    cin >> gama;
    double start = GetTickCount(); // 开始计时
    for (int i = 0; i < total_steps; i++)
    {
        x_list.push_back(x_p[0]);
        p_list.push_back(x_p[0]);
        //OutFile1 << x_p[0][0] << ", " << x_p[1][0]
        //         << "\n";
        if (i % 1000 == 0)
        {
            cout << i << " steps, " << GetTickCount() - start << "ms has been used, total steps are " << total_steps << endl;
        }
        x_p = BAOAB(x_p, 0.01, gama, beta);
    }
    double k_e = ensemble_average(p_list, kinetic_energy);
    double v_e = ensemble_average(x_list, V_potential);
    cout
        << "---------------RESULT---------------" << endl;
    cout << "the kinetic enenrgy of the system is " << k_e << endl;
    cout << "the potential energy of the system is " << v_e << endl;
    cout << "total time is " << GetTickCount() - start << " ms" << endl;
    cout << "---------------END---------------" << endl;

    double i;
    cin >> i;
    return 0;
}