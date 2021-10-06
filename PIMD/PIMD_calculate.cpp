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

int N = 1;        // 系统的自由度
double h_bar = 1; // 约化Planck常数
//vector<vector<double>> M = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; // 质量矩阵
vector<vector<double>> M = {{1}};

// 设置随机数种子
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> distribution(0.0, 1.0);

/*
staging变换，将Cartesian坐标变换成为staging坐标

x: vector<vector<double>> 输入的Cartesian坐标，嵌套列表，可以支持高维体系的变换
return: vector<<vector<double>> 返回staging坐标中的位置
经过测试无误
*/
vector<vector<double>> staging_transformation(vector<vector<double>> x)
{
    vector<vector<double>> staging;
    for (int i = 0; i < x.size(); i++) // x.size() 表示了beads的数量
    {
        if (i == 0)
        {
            staging.push_back(x[i]);
        }
        else if (0 < i and i < x.size() - 1)
        {
            vector<double> vec; // 构建某个bead的坐标
            for (int j = 0; j < x[i].size(); j++)
            {
                vec.push_back(x[i][j] - (double)i / (i + 1) * x[i + 1][j] - x[0][j] / (i + 1));
            } // 注意这里要将i转换为double
            staging.push_back(vec);
        }
        else
        {
            vector<double> vec;
            for (int j = 0; j < x[i].size(); j++)
            {
                vec.push_back(x[i][j] - x[0][j]);
            }
            staging.push_back(vec);
        }
    }
    return staging;
}

/*
staging逆变换，将staging坐标转换成为Cartesian坐标
s: vector<vector<double>> 输入的staging坐标
return: vector<vector<double>> 返回的Cartesian坐标
经过测试无误
*/
vector<vector<double>> staging_transformation_inv(vector<vector<double>> s)
{
    vector<vector<double>> x;
    for (int i = 0; i < s.size(); i++)
    {
        if (i == 0)
        {
            x.push_back(s[i]);
        }
        else
        {
            vector<double> vec;
            for (int j = 0; j < s[i].size(); j++)
            {
                double cd = s[0][j];
                for (int l = i; l < s.size(); l++)
                {
                    cd = cd + (double)i / l * s[l][j];
                }
                vec.push_back(cd);
            }
            x.push_back(vec);
        }
    }
    return x;
}

/*
在终端中打印 vector<vector<<double>>中的元素
*/
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
    double x_4 = 0;
    for (int j = 0; j < x.size(); j++)
    {
        x_4 += x[j] * x[j] * x[j] * x[j];
    }
    return 0.25 * x_4;
}

/*
staging坐标下的所有beads的外势
s_chain: vector<vector<double>> staging坐标下所有beads的坐标
return: double 
*/
double Phi_potential(vector<vector<double>> s_chain)
{
    int p = s_chain.size();
    vector<vector<double>> x_chain = staging_transformation_inv(s_chain);
    double phi = 0;
    for (int j = 0; j < p; j++)
    {
        phi += V_potential(x_chain[j]);
    }
    return (double)1 / p * phi;
}

/*
Cartesian坐标下的势能
x_chain: vector<vector<double>> Cartesian坐标
return: double 势能
*/
double potential(vector<vector<double>> x_chain, double beta)
{
    int p = x_chain.size();
    double potential = 0;
    for (int j = 0; j < p; j++)
    {
        potential += V_potential(x_chain[j]);
    }
    return (double)1 / p * potential;
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
        f.push_back(x[j] * x[j] * x[j]);
    }
    return f;
}

/*
staging坐标下所有beads对应的势能梯度，在求解Langevin方程中直接使用
s_chain: vector<vector<double>> 所有beads的staging坐标
return: vector<vector<double>>
*/
vector<vector<double>> nabla_Phi_potential(vector<vector<double>> s_chain)
{
    vector<vector<double>> phi_drv;
    vector<vector<double>> x_chain = staging_transformation_inv(s_chain);
    // 将staging坐标变换为Cartesian坐标，方便计算势能
    {
        vector<double> drv; // 计算第一个bead的梯度，为后面递推做准备
        for (int k = 0; k < s_chain[0].size(); k++)
        {
            double v = 0;
            for (int j = 0; j < s_chain.size(); j++)
            {
                v += nabla_V_potential(x_chain[j])[k];
            }
            drv.push_back((double)1 / s_chain.size() * v);
        }
        phi_drv.push_back(drv);
    }
    // 后面的每个bead的势能梯度都需要递推计算

    for (int i = 1; i < s_chain.size(); i++)
    {
        vector<double> drv;
        for (int j = 0; j < s_chain[0].size(); j++)
        {
            drv.push_back((double)1 / s_chain.size() * nabla_V_potential(x_chain[i])[j] + (double)(i - 1) / i * phi_drv[i - 1][j]);
        }
        phi_drv.push_back(drv);
    }
    return phi_drv;
}

/*
动能的primitive estimator，Cartesian坐标下使用
x_chain: vector<vector<double>> 所有beads的Cartesian坐标
return: double 动能
*/
double kinetic_energy_primitive(vector<vector<double>> x_chain, double beta)
{
    int n = x_chain[0].size();
    int p = x_chain.size();
    double T = 0;
    for (int j = 0; j < p; j++)
    {
        double v = 0;
        for (int k = 0; k < n; k++)
        {
            v += (x_chain[(j + 1) % p][k] - x_chain[j][k]) * M[k][k] * (x_chain[(j + 1) % p][k] - x_chain[j][k]);
        }
        T += v;
    }

    return (double)-T * p / (2 * beta * beta * h_bar * h_bar) + (double)n * p / (2 * beta);
}

/*
动能的viral estimator, Cartesian坐标下使用
x_chain: vector<vector<double>> 所有beads的Cartesian坐标
return: double 动能
*/
double kinetic_energy_viral(vector<vector<double>> x_chain, double beta)
{
    int p = x_chain.size();    // beads数量
    int n = x_chain[0].size(); // 系统的维数
    double T = 0;
    vector<double> x_c; // 构造beads质心的坐标
    for (int k = 0; k < n; k++)
    {
        double crd = 0;
        for (int j = 0; j < p; j++)
        {
            crd += x_chain[j][k];
        }
        x_c.push_back((double)1 / p * crd);
    }
    for (int j = 0; j < p; j++)
    {
        for (int k = 0; k < n; k++)
        {
            T += (x_chain[j][k] - x_c[k]) * nabla_V_potential(x_chain[j])[k];
        }
    }
    return (double)1 / (2 * p) * T + (double)n / (2 * beta);
}

/*
给定估计量的系综平均值
estimator: func 估计量
config_list: K * P * N 给定的一定数量的构型分布，只包含各个beads的位置信息，不包含动量，
             staging坐标，但是estimator往往是定义在Cartesian坐标下，注意转换
return: double

// 注意这里传递参数尽量传递引用
*/
double ensemble_average(vector<vector<vector<double>>> &config_list, double (&estimator)(vector<vector<double>>, double beta), double beta)
{
    double average = 0;
    int total = config_list.size(); // 所选取的构型的数量
    for (int i = 0; i < total; i++)
    {
        average += estimator(staging_transformation_inv(config_list[i]), beta);
    }
    return (double)average / total;
}

/*
以OBABO的方式演化Langevin方程
pre: 2*P*N维列表，pre[0]代表所有beads的位置，pre[1]代表所有beads的动量
dt: double 代表时间间隔
return: 2*P*N维列表，表示位置和动量
*/
vector<vector<vector<double>>> OBABO(vector<vector<vector<double>>> pre, double dt, double beta)
{
    // 首先定义一些量
    vector<vector<vector<double>>> nxt;
    int n = pre[0][0].size();                  // 空间维数
    int p = pre[0].size();                     // beads数量
    double omega_p = sqrt(p) / (beta * h_bar); // 内势力常数
    // 在staging变换下，所有每一个bead的所有自由度的最优阻力系数相同(ωp)
    // 定义了在O过程中出现的两个常数
    double c1 = exp(-omega_p * dt / 2);
    double c2 = sqrt(1 - c1 * c1);

    vector<vector<double>> s_chain = pre[0];
    vector<vector<double>> p_chain = pre[1]; // 分别提取位置和动量信息

    // 首先演化半步(dt/2) O项，即热库项
    for (int j = 0; j < p; j++)
    {
        if (j == 0)
        {
            for (int k = 0; k < n; k++)
            {
                p_chain[0][k] = c1 * p_chain[0][k] + c2 * sqrt(1 / beta) * sqrt(M[k][k]) * distribution(generator);
            }
        }
        else
        {
            for (int k = 0; k < n; k++)
            {
                p_chain[j][k] = c1 * p_chain[j][k] + c2 * sqrt(1 / beta) * sqrt((double)(j + 1) / j * M[k][k]) * distribution(generator);
            }
        }
    }

    // 其次演化半步B项（动量）
    for (int j = 0; j < p; j++)
    {
        for (int k = 0; k < n; k++)
        {
            p_chain[j][k] = p_chain[j][k] - dt / 2 * nabla_Phi_potential(s_chain)[j][k];
        }
    }

    // 然后演化一步内势力A项
    for (int j = 0; j < p; j++)
    {
        if (j == 0)
        {
            for (int k = 0; k < n; k++)
            {
                s_chain[0][k] = s_chain[0][k] + 1 / M[k][k] * p_chain[0][k] * dt;
            }
        }
        else
        {
            for (int k = 0; k < n; k++)
            {
                s_chain[j][k] = s_chain[j][k] * cos(omega_p * dt) + 1 / omega_p * sin(omega_p * dt) * j / ((j + 1) * M[k][k]) * p_chain[j][k];
                p_chain[j][k] = -omega_p * sin(omega_p * dt) * (j + 1) / j * M[k][k] * s_chain[j][k] + cos(omega_p * dt) * p_chain[j][k];
            }
        }
    }

    // 再演化半步B动量项
    for (int j = 0; j < p; j++)
    {
        for (int k = 0; k < n; k++)
        {
            p_chain[j][k] = p_chain[j][k] - dt / 2 * nabla_Phi_potential(s_chain)[j][k];
        }
    }

    // 最后演化半步O项
    for (int j = 0; j < p; j++)
    {
        if (j == 0)
        {
            for (int k = 0; k < n; k++)
            {
                p_chain[0][k] = c1 * p_chain[0][k] + c2 * sqrt(1 / beta) * sqrt(M[k][k]) * distribution(generator);
            }
        }
        else
        {
            for (int k = 0; k < n; k++)
            {
                p_chain[j][k] = c1 * p_chain[j][k] + c2 * sqrt(1 / beta) * sqrt((double)(j + 1) / j * M[k][k]) * distribution(generator);
            }
        }
    }
    nxt.push_back(s_chain);
    nxt.push_back(p_chain);
    return nxt;
}

/*
以BAOAB的方式演化Langevin方程
pre: 2*P*N维列表，pre[0]代表所有beads的位置，pre[1]代表所有beads的动量
dt: double 时间步长
return: 2*P*N维列表，表示位置和动量
*/
vector<vector<vector<double>>> BAOAB(vector<vector<vector<double>>> pre, double dt, double beta)
{
    // 首先定义一些量
    vector<vector<vector<double>>> nxt;
    int n = pre[0][0].size();                          // 空间维数
    int p = pre[0].size();                             // beads数量
    double omega_p = (double)sqrt(p) / (beta * h_bar); // 内势力常数
    // 在staging变换下，所有每一个bead的所有自由度的最优阻力系数相同(ωp)
    // 定义了在O过程中出现的两个常数 注意这里的时间步长与OBABO不一样
    double c1 = exp(-omega_p * dt);
    double c2 = sqrt(1 - c1 * c1);

    vector<vector<double>> s_chain = pre[0];
    vector<vector<double>> p_chain = pre[1]; // 分别提取位置和动量信息

    // 首先演化半步B动量项
    for (int j = 0; j < p; j++)
    {
        for (int k = 0; k < n; k++)
        {
            p_chain[j][k] = p_chain[j][k] - dt / 2 * nabla_Phi_potential(s_chain)[j][k];
        }
    }

    // 其次演化<半步>A内势力项
    for (int j = 0; j < p; j++)
    {
        if (j == 0)
        {
            for (int k = 0; k < n; k++)
            {
                s_chain[0][k] = s_chain[0][k] + 1 / M[k][k] * p_chain[0][k] * dt / 2;
            }
        }
        else
        {
            for (int k = 0; k < n; k++)
            {
                s_chain[j][k] = s_chain[j][k] * cos(omega_p * dt / 2) + 1 / omega_p * sin(omega_p * dt / 2) * j / ((j + 1) * M[k][k]) * p_chain[j][k];
                p_chain[j][k] = -omega_p * sin(omega_p * dt / 2) * (j + 1) / j * M[k][k] * s_chain[j][k] + cos(omega_p * dt / 2) * p_chain[j][k];
            }
        }
    }

    // 然后演化一步O控温项

    for (int j = 0; j < p; j++)
    {
        if (j == 0)
        {
            for (int k = 0; k < n; k++)
            {
                p_chain[0][k] = c1 * p_chain[0][k] + c2 * sqrt(1 / beta) * sqrt(M[k][k]) * distribution(generator);
            }
        }
        else
        {
            for (int k = 0; k < n; k++)
            {
                p_chain[j][k] = c1 * p_chain[j][k] + c2 * sqrt(1 / beta) * sqrt((double)(j + 1) / j * M[k][k]) * distribution(generator);
            }
        }
    }

    // 再演化半步A内势力项
    for (int j = 0; j < p; j++)
    {
        if (j == 0)
        {
            for (int k = 0; k < n; k++)
            {
                s_chain[0][k] = s_chain[0][k] + 1 / M[k][k] * p_chain[0][k] * dt / 2;
            }
        }
        else
        {
            for (int k = 0; k < n; k++)
            {
                s_chain[j][k] = s_chain[j][k] * cos(omega_p * dt / 2) + 1 / omega_p * sin(omega_p * dt / 2) * j / ((j + 1) * M[k][k]) * p_chain[j][k];
                p_chain[j][k] = -omega_p * sin(omega_p * dt / 2) * (j + 1) / j * M[k][k] * s_chain[j][k] + cos(omega_p * dt / 2) * p_chain[j][k];
            }
        }
    }

    // 最后演化半步B动量项
    for (int j = 0; j < p; j++)
    {
        for (int k = 0; k < n; k++)
        {
            p_chain[j][k] = p_chain[j][k] - dt / 2 * nabla_Phi_potential(s_chain)[j][k];
        }
    }
    nxt.push_back(s_chain);
    nxt.push_back(p_chain);
    return nxt;
}

/*
为了判断是否达到平衡，还需要一个计算某段时间内平均温度的函数
p_list: K * P * N 给定数目的动量构型
return: double 对应温度的β值
*/
double temperature(vector<vector<vector<double>>> &p_list)
{
    double e = 0;
    int n = p_list[0][0].size();
    int p = p_list[0].size();
    int K = p_list.size();
    for (int k = 0; k < K; k++)
    {
        for (int j = 0; j < p; j++)
        {
            if (j == 0)
            {
                for (int l = 0; l < n; l++)
                {
                    e += p_list[k][j][l] * 1 / M[l][l] * p_list[k][j][l];
                }
            }
            else
            {
                for (int l = 0; l < n; l++)
                {
                    e += p_list[k][j][l] * j / (M[l][l] * j + 1) * p_list[k][j][l];
                }
            }
        }
    }
    e = e / (n * p * K);
    return e;
}

int main()
{
    int P; // PIMD中的bead数量
    cout << "the nunber of beads = ?" << endl;
    cin >> P;
    double beta; // 系统设定的温度
    cout << "beta = ?" << endl;
    cin >> beta;
    int total_steps; // 演化的步数
    cout << "the number of total steps = ?" << endl;
    cin >> total_steps;

    stringstream fmt1;
    fmt1 << "quartic_P_64"
         << ".txt";
    ofstream OutFile1(fmt1.str());

    stringstream fmt2;
    fmt2 << "quartic_P_64_temperature"
         << ".txt";
    ofstream OutFile2(fmt2.str());

    double start = GetTickCount(); // 开始计时
    for (int t = 0; t < 100; t++)
    {
        double T = 0.1 + (double)t / 10;
        beta = 1 / T;
        vector<double> x_i(N, 0.05);
        vector<vector<double>> x(P, x_i);
        vector<double> p_i(N, 0);
        vector<vector<double>> p(P, p_i); // 构造初始位置和动量
        cout << "------------STEP " << t << "--------------" << endl;
        cout << "-----------BEGIN STEP-----------" << endl;
        cout << "the dimension of the system " << N << endl;
        cout << "the number of beads " << P << endl;
        cout << "the temperature of the system is " << T << endl;
        cout << "initial positon in Cartesian coordinate:" << endl;
        print_vector_2(x);
        cout << "initial momentum in Cartesian coordinate: " << endl;
        print_vector_2(p);
        vector<vector<double>> s = staging_transformation(x);

        vector<vector<vector<double>>> s_p = {x, p};
        vector<vector<vector<double>>> p_list;
        vector<vector<vector<double>>> s_list; // 储存构型分布
        vector<double> beta_list;
        for (int i = 0; i < total_steps; i++)
        {
            s_list.push_back(s_p[0]);
            p_list.push_back(s_p[1]);
            if (i % 1000 == 0)
            {
                cout << i << " steps, " << GetTickCount() - start << "ms has been used, total steps are " << total_steps << endl;
            }
            if (i % 500 == 0)
            {
                OutFile2 << std::setiosflags(std::ios::scientific | std::ios::showpos) << temperature(p_list) << " ";
            }
            s_p = BAOAB(s_p, 0.1*sqrt(P), beta);
        }
        OutFile2 << endl;

        cout << "---------------REASULT---------------" << endl;
        OutFile1 << beta << " ";

        double KE_primitive = ensemble_average(s_list, kinetic_energy_primitive, beta);
        cout << "the kinetic energy(primitive) at " << T << " is: ";
        cout << std::setiosflags(std::ios::scientific | std::ios::showpos) << KE_primitive << endl;
        OutFile1 << KE_primitive << " ";

        double KE_viral = ensemble_average(s_list, kinetic_energy_viral, beta);
        cout << "the kinetic energy(viral) at " << T << " is: ";
        cout << KE_viral << endl;
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos) << KE_viral << " ";

        double POT = ensemble_average(s_list, potential, beta);
        cout << "the potential energy is: ";
        cout << POT << endl;
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos) << POT << " ";
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos) << POT + KE_viral << " ";
        OutFile1 << endl;

        double TEP = temperature(p_list);
        cout << "the temperature is: ";
        cout << TEP << endl;

        cout << "total time is: " << GetTickCount() - start << " ms" << endl;
        cout << "---------------END STEP---------------" << endl;
    }
    OutFile1.close();
    fmt1.clear();
    cout << "-------------END  ALL-------------" << endl;

    int i;
    cin >> i;
    return 0;
}