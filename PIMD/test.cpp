#include <iostream>
using std::cin;
using std::cout;
using std::endl;

const int P = 8;
const int N = 1;
double nabla_V[N];
double nabla_phi[P][N];

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

void nabla_V_potential(double (&x)[N])
{

    for (int j = 0; j < N; j++)
    {
        nabla_V[j] = x[j];
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
    print_vector_2(s_chain);
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

int main()
{
    double var[P][N] = {1,1, 2, 2, 3, 3, 4, 4};
    nabla_Phi_potential(var);
    print_vector_2(nabla_phi);
    print_vector_2(var);
    int i;
    cin >> i;
    return 0;
}