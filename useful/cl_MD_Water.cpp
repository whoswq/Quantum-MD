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
const int dim = 3; // 空间维数，均为3维
const int n_atom = 3;  // 原子个数
const double dt_MD = 2;                      // Langevin动力学的时间步长
const double dt_HD = 2;                    // Hamilton动力学的时间步长
const double T = 300;                        // 设定温度
const double kb = 1.28065e-23 / 4.35974e-18; // Boltzman常数
const double beta = 1.0 / (kb * T);
const double gama = 0.2;     // 摩擦系数 怎么选择比较合适？
const int steps_MD = 500000; // 分子动力学（采样）的总步长
const int n_traj = 100;     // Hamilton动力学时选取的轨线数目
const int steps_HD = 200000;    // Hamilton动力学总步长
const int step_leap = 5;      // 计算关联函数时跳过的步数
const bool cstr_p = true;
const bool cstr_f = true;
const int cstr_f_step = 100;  // 表示每100步约束一次力
const int max_ave = 150; // 时间关联函数的最大平均步长
int step_now = 0; // 记录当前HD演化的步数，为了实现每隔cstr_f_step对力约束
double IDE[dim * dim] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
double x[dim * n_atom] = { 0 };
double p[dim * n_atom] = { 0 }; // 储存演化时的中间变量
double x_traj[dim * n_atom * n_traj] = { 0 };
double p_traj[dim * n_atom * n_traj] = { 0 };                         // 储存用于实时间动力学的轨线初始值
double nabla_V[dim * n_atom] = { 0 };                                 // 储存势能梯度
double f[dim * n_atom] = { 0 };                                       // 储存受力
double* x_list_traj; // 用于计算关联函数的数据的列表
double* p_list_traj; // 用mkl的函数申请内存可以避开数组长度的限制
double t_array[steps_HD / step_leap] = { 0 }; // 储存每一帧的时间
double p_corr[steps_HD / step_leap] = { 0 };
double x_corr[steps_HD / step_leap] = { 0 }; // 储存时间关联函数，暂时未使用
double M[dim * n_atom] = { 1837.7, 1837.7, 1837.7, 16 * 1837.7, 16 * 1837.7, 16 * 1837.7, 1837.7, 1837.7, 1837.7 };

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
// 第一个参数为高斯分布的平均值，第二个参数为标准差
std::normal_distribution<double> distribution(0.0, 1.0);

void print_vector(double* x, int n) {
    for (int k = 0; k < n; k++) {
        cout << x[k] << " ";
    }
    cout << endl;
}

/*
对位置的约束
将水分子放在质心系
*/
void constrain_x(double* x) {
    double x_c[3] = { 0 }; // 质心坐标
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
用于检查constrain_x
*/
double center(double* x) {
    double var = 0;
    for (int k = 0; k < n_atom * dim; k++) {
        var += M[k] * x[k];
    }
    return var;
}

/*
计算惯量张量
计算结果储存在double指针I指向的位置
*/
void I_tensor(double* x, double* I) {
    double r_2[n_atom] = { 0 }; // 储存原子距原点的距离
    for (int n = 0; n < n_atom; n++) {
        double var = 0;
        for (int k = 0; k < dim; k++) {
            var += x[n * dim + k] * x[n * dim + k];
        }
        r_2[n] = var;
    }
    int d_ij = 0;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if (i == j) {
                d_ij = 1;
            }
            else {
                d_ij = 0;
            }
            double var = 0;
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
对力的约束
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

/*
水分子的势函数
*/
extern "C" {
    void h2opot_dv_(double* v, double* dvdx, double* x, double* xd, int* n);
}
void nabla_V_Water(double* x, bool cstr = false) {
    double v = 0;
    int n = 9;
    double xd[9 * 9] = { 0 };
    for (int k = 0; k < n_atom * dim; k++) {
        xd[9 * k + k] = 1;
    }
    h2opot_dv_(&v, nabla_V, x, xd, &n);
    if (cstr and (step_now / cstr_f_step == 0)) {
        constrain_force(x, nabla_V);
    }
}

const double c1 = exp(-gama * dt_MD);
const double c2 = sqrt(1.0 - c1 * c1);
/*
以BAOAB的方式演化Langevin方程
pre: 2*N维列表，pre[0]代表位置，pre[1]代表动量
dt: double 时间步长
gama: double Langevin动力学中的阻力常数
*/
void BAOAB(double* x,
    double* p,
    void (*nabla_V_potential)(double*, bool),
    double dt,
    double c1,
    double c2,
    bool cstr)
{
    nabla_V_potential(x, cstr);
    for (int i = 0; i < dim * n_atom; i++) // 半步动量
    {
        p[i] = p[i] - dt * 0.5 * nabla_V[i];
    }
    for (int i = 0; i < dim * n_atom; i++) // 半步位置
    {
        x[i] = x[i] + p[i] / M[i] * dt * 0.5;
    }
    for (int i = 0; i < dim * n_atom; i++) // 控温
    {
        p[i] = c1 * p[i] + c2 * sqrt(M[i] / beta) * distribution(generator);
    }
    for (int i = 0; i < dim * n_atom; i++)
    {
        x[i] = x[i] + p[i] / M[i] * dt * 0.5;
    }
    nabla_V_potential(x, cstr);
    for (int i = 0; i < dim * n_atom; i++)
    {
        p[i] = p[i] - dt * 0.5 * nabla_V[i];
    }
}

/*
哈密顿动力学
*/
void Velocity_Verlet(double* x,
    double* p,
    void (*nabla_V_potential)(double*, bool),
    double dt,
    bool cstr)
{
    nabla_V_potential(x, cstr);
    for (int i = 0; i < dim * n_atom; i++) // 半步动量
    {
        p[i] = p[i] - dt * 0.5 * nabla_V[i];
    }
    for (int i = 0; i < dim * n_atom; i++) // 一步位置
    {
        x[i] = x[i] + p[i] / M[i] * dt;
    }
    nabla_V_potential(x, cstr);
    for (int i = 0; i < dim * n_atom; i++)
    {
        p[i] = p[i] - dt * 0.5 * nabla_V[i];
    }
}

/*
计算指定构型的温度
*/
double temperature(double* p)
{
    double T = 0;
    for (int k = 0; k < dim * n_atom; k++)
    {
        T += p[k] * p[k] / M[k];
    }
    return T / (dim * n_atom * kb);
}

/*
计算给定轨线集合的位置时间关联函数，未考虑在轨线上的时间平均
*/
double x_corr_pmt(double* x_list, int begin, int end)
{
    double corr = 0;
    for (int j = 0; j < n_traj; j++)
    {
        double var = 0;
        for (int k = 0; k < dim * n_atom; k++)
        {
            var += x_list[n_atom * dim * steps_HD / step_leap * j + begin * n_atom * dim + k] *
                x_list[n_atom * dim * steps_HD / step_leap * j + end * n_atom * dim + k];
        }
        corr += var;
    }
    return corr / n_traj;
}
/*
计算给定轨线集合的位置时间关联函数，考虑了轨线上的时间平均
*/
double x_corr_func(double* x_list, int t, int cnt = 50) {
    double ave = 0;
    for (int i = 0; i < cnt; i++) {
        ave += x_corr_pmt(x_list, i, i + t);
    }
    return ave / cnt;
}

/*
计算给定轨线集合的动量时间关联函数，未考虑在轨线上的时间平均
*/
double p_corr_pmt(double* p_list, int begin, int end)
{
    double corr = 0;
    for (int j = 0; j < n_traj; j++)
    {
        double var = 0;
        for (int k = 0; k < dim * n_atom; k++)
        {
            var += p_list[n_atom * dim * steps_HD / step_leap * j + begin * n_atom * dim + k] *
                p_list[n_atom * dim * steps_HD / step_leap * j + end * n_atom * dim + k];
        }
        corr += var;
    }
    return corr / n_traj;
}

/*
计算给定轨线集合的动量时间关联函数，考虑时间平均
*/
double p_corr_func(double* p_list, int t, int cnt = 50) {
    double ave = 0;
    for (int i = 0; i < cnt; i++) {
        ave += p_corr_pmt(p_list, i, i + t);
    }
    return ave / cnt;
}

int main()
{   // 势能面的极小值点是从作业中，通过共轭梯度法搜索出来的
    double x[n_atom * dim] = { 1.6774683, -0.04038825, 0, -0.1297673, 0.06031243, 0, -0.48502213, 1.83515775, 0 };
    cout << "===================="
        << "START LANGEVIN DYNAMICS"
        << "====================" << endl;
    x_list_traj = (double*)MKL_malloc((steps_HD / step_leap) * n_atom * dim * n_traj * sizeof(double), 128);
    p_list_traj = (double*)MKL_malloc((steps_HD / step_leap) * n_atom * dim * n_traj * sizeof(double), 128);
    time_t begin, end;
    begin = clock();
    double tem = 0;
    for (int step = 0; step < steps_MD; step++)
    {
        BAOAB(x, p, nabla_V_Water, dt_MD, c1, c2, false);
        tem += temperature(p);
        if (step % 100000 == 0)
        {
            end = clock();
            cout << "In steps " << step << ", time = " << (double)(end - begin) / CLOCKS_PER_SEC
                << ", temperature = " << double(tem / (step + 1)) << endl;
        }
        if (step >= steps_MD - n_traj)
        {
            for (int k = 0; k < dim * n_atom; k++)
            {
                x_traj[(step + n_traj - steps_MD) * dim * n_atom + k] = x[k];
                p_traj[(step + n_traj - steps_MD) * dim * n_atom + k] = p[k];
            }
        }
    }
    cout << endl;
    cout << "====================="
        << "END LANGEVIN DYNAMICS"
        << "=====================" << endl;
    cout << "===================="
        << "START HAMILTON DYNAMICS"
        << "====================" << endl;
    // 构建每一帧的时间
    for (int n = 0; n < steps_HD / step_leap - max_ave; n++)
    {
        t_array[n] = n * dt_HD * step_leap;
    }

    for (int j = 0; j < n_traj; j++)
    {
        // j表示在第j根轨线，首先修改x，p为初始值
        step_now = 0;
        if (j % 5 == 0)
        {
            end = clock();
            cout << "in trajectory " << j << ", time = " << (double)(end - begin) / CLOCKS_PER_SEC << endl;
        }
        for (int k = 0; k < dim * n_atom; k++)
        {
            x[k] = x_traj[dim * n_atom * j + k];
            p[k] = p_traj[dim * n_atom * j + k];
        }
        // 此时对初始动量进行约束
        if (cstr_p)
        {
            constrain_p(x, p);
        }
        for (int step = 0; step < steps_HD; step++)
        {
            // 首先应该在x_list_traj, p_list_traj中储存轨线，按每条轨线的顺序
            // 存储，即[N * steps_HD / step_leap * j]为第j条轨线的头指针
            if (step % step_leap == 0)
            {
                for (int k = 0; k < dim * n_atom; k++)
                {
                    x_list_traj[steps_HD / step_leap * dim * n_atom * j + dim * n_atom * (step / step_leap) + k] = x[k];
                    p_list_traj[steps_HD / step_leap * dim * n_atom * j + dim * n_atom * (step / step_leap) + k] = p[k];
                }
            }
            if (step_now % cstr_f_step == 0){
                constrain_p(x, p);
            }
            Velocity_Verlet(x, p, nabla_V_Water, dt_HD, cstr_f);
            step_now += 1;
        }
    }
    cout << "====================="
        << "END HAMILTON DYNAMICS"
        << "=====================" << endl;
    // 计算时间关联函数
    stringstream fmt1; // 通过文件存储时间关联函数的数据
    fmt1 << "cl_water_p_cstr_f_cstr"
        << ".txt";
    stringstream fmt2; // 通过文件存储时间关联函数的数据
    fmt2 << "cl_water_x_cstr_f_cstr"
        << ".txt";
    ofstream OutFile1(fmt1.str());
    ofstream OutFile2(fmt2.str());
    cout << "need some time to calculate correlation function" << endl;
    for (int t = 0; t < steps_HD / step_leap; t++)
    {
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
            << t_array[t] << "  ";
        OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
            << p_corr_pmt(p_list_traj, 0, t) << "  \n"; // p_list_traj[N * steps_HD / step_leap * 200 +  N * t + 2]
        OutFile2 << std::setiosflags(std::ios::scientific | std::ios::showpos)
            << t_array[t] << "  ";
        OutFile2 << std::setiosflags(std::ios::scientific | std::ios::showpos)
            << x_corr_func(x_list_traj, t) << "  \n";
        if (t % 1000 == 0) {
            end = clock();
            cout << "steps = " << t << ", time = " << double(end - begin) / CLOCKS_PER_SEC << endl;
        }
    }
    OutFile1.close();
    fmt1.clear();
    OutFile2.close();
    fmt2.clear();
    cout << "==================="
        << "END ALL"
        << "====================" << endl;
    int i;
    cin >> i;
    return 0;
}
