import numpy as np
import time  # 计算程序运行时间
"""
start = time.time()
t = time.time() - start
"""
# 定义一些常数
P = 1  # beads的数量
N = 1  # 系统的维数
M = np.array([1])  # 系统的质量矩阵，由于是对角矩阵，只储存了对角元
M_inv = 1 / M  # 质量矩阵的逆，同样只储存了对角元，矩阵乘法时要格外注意
beta = 0.01  # 温度的倒数，不要超过15
steps = 200000  # 演化的总步数
hbar = 1  # Planck常数
omega_P = np.sqrt(P) / (beta * hbar)  # 内势力的角频率
dt = 0.1 / omega_P  # 时间步长


def staging_transf(x_array):
    """
    x_array: numpy 2-D array 所有beads的Cartesian坐标
    retrun: numpy 2-D array 所有beads的staging坐标
    """
    s_array = np.zeros((P, N))
    s_array[0] += x_array[0]
    for j in range(1, P):
        s_array[j] += x_array[j] - (j * x_array[
            (j + 1) % P] + x_array[0]) / (j + 1)
    return s_array


def staging_transf_inv(s_array):
    """
    s_array: numpy 2_D array 所有beads的staging坐标
    return: numpy 2-D array 所有beads的Cartesian坐标
    """
    x_array = np.zeros((P, N))
    x_array[0] += s_array[0]
    for j in range(P - 1, 0, -1):
        x_array[j] += s_array[j] + (j * x_array[
            (j + 1) % P] + s_array[0]) / (j + 1)
    return x_array


def V_potential(x):
    """
    x: numpy 1-D array Cartesian坐标下单个bead的坐标
    return: float 单个bead的势能
    """
    return 0.5 * x.dot(x)


def nabla_V_potential(x):
    """
    x: numpy 1-D array Cartesian坐标下单个bead的坐标
    return: numpy 1-D array 受力
    """
    return x


def phi_potential(s_array):
    """
    x_array: numpy 2-D array staging坐标下所有beads的位置
    return: float 势能
    """
    x_array = staging_transf_inv(s_array)
    phi = 0
    for j in range(P):
        phi += V_potential(x_array[j])
    return phi / P


def nabla_phi_potential(s_array):
    """
    s_array: numpy 2-D array 所有beads的staging坐标
    return: numpy 2-D array 所有beads在staging坐标下的受力
    """
    x_array = staging_transf_inv(s_array)
    s_force_array = np.zeros((P, N))
    force0 = nabla_V_potential(x_array[0])
    # 构造第一个bead的受力
    for j in range(1, P):
        force0 += nabla_V_potential(x_array[j])
    s_force_array[0] += 1 / P * force0
    # 递推其他bead的受力
    for j in range(1, P):
        s_force_array[j] += 1 / P * nabla_V_potential(
            x_array[j]) + (j - 1) / j * s_force_array[j - 1]
    return s_force_array


def kinetic_energy_viral(x_array):
    """
    x_array: numpy 2-D array 所有beads的Cartesian坐标
    return: float 系统的动能（通过viral定理估计）
    """
    x_c = np.zeros(N)
    for j in range(P):
        x_c += x_array[j]
    x_c = x_c / P  # 构造质心
    K = 0
    for j in range(P):
        K += (x_array[j] - x_c).dot(nabla_V_potential(x_array[j]))
    return N / (2 * beta) + K / (2 * P)


# 记录一些演化时会用到的变量
c1 = np.exp(-omega_P * dt)
c2 = np.sqrt((1 - c1 * c1) / beta)
d1 = np.cos(omega_P * dt * 0.5)
d2 = np.sin(omega_P * dt * 0.5)


def BAOAB(pre):
    """
    使用BAOAB的算法演化Lagevin方程
    pre: 前一个时刻的位置和动量 staging坐标下
    return: 下一个时刻的位置和动量 staging坐标下
    """
    s_array = pre[0]
    p_array = pre[1]
    # 首先演化半步动量
    p_array = p_array - nabla_phi_potential(s_array) * dt * 0.5
    # 其次演化半步内势力
    s_array[0] = s_array[0] + M_inv * p_array[0] * dt * 0.5  # 注意质量矩阵与动量对应元素相乘
    for j in range(1, P):
        s_j = d1 * s_array[j] + d2 / omega_P * j / (j + 1) * M_inv * p_array[j]
        p_array[j] = -omega_P * d2 * (j + 1) / j * M * s_array[j] + np.cos(
            omega_P * dt * 0.5) * p_array[j]
        s_array[j] = s_j
    # 演化一步控温
    p_array[
        0] = c1 * p_array[0] + c2 * np.sqrt(M) * np.random.standard_normal(N)
    for j in range(1, P):
        p_array[j] = c1 * p_array[j] + c2 * np.sqrt(
            (j + 1) / j * M) * np.random.standard_normal(N)
    # 再演化半步内势力
    s_array[0] = s_array[0] + M_inv * p_array[0] * dt * 0.5
    for j in range(1, P):
        s_j = d1 * s_array[j] + d2 / omega_P * j / (j + 1) * M_inv * p_array[j]
        p_array[j] = -omega_P * d2 * (j + 1) / j * M * s_array[j] + np.cos(
            omega_P * dt * 0.5) * p_array[j]
        s_array[j] = s_j
    # 最后演化半步动量
    p_array = p_array - nabla_phi_potential(s_array) * dt * 0.5

    return s_array, p_array


# 记录一些演化时会用到的变量
o1 = np.exp(-omega_P * dt * 0.5)
o2 = np.sqrt((1 - c1 * c1) / beta)


def OABAO(pre):
    """
    使用OABAO的算法演化Lagevin方程
    pre: numpy 3-D array 前一个时刻的位置和动量 staging坐标下
    return: numpy 3-D array 下一个时刻的位置和动量 staging坐标下
    """
    s_array = pre[0]
    p_array = pre[1]

    # 首先演化半步控温
    p_array[
        0] = o1 * p_array[0] + o2 * np.sqrt(M) * np.random.standard_normal(N)
    for j in range(1, P):
        p_array[j] = o1 * p_array[j] + o2 * np.sqrt(
            (j + 1) / j * M) * np.random.standard_normal(N)
    # 再演化半步内势力
    s_array[0] = s_array[0] + M_inv * p_array[0] * dt * 0.5
    for j in range(1, P):
        s_j = d1 * s_array[j] + d2 / omega_P * j / (j + 1) * M_inv * p_array[j]
        p_array[j] = -omega_P * d2 * (j + 1) / j * M * s_array[j] + np.cos(
            omega_P * dt * 0.5) * p_array[j]
        s_array[j] = s_j
    # 再演化一步动量
    p_array = p_array - nabla_phi_potential(s_array) * dt
    # 再演化半步内势力
    s_array[0] = s_array[0] + M_inv * p_array[0] * dt * 0.5
    for j in range(1, P):
        s_j = d1 * s_array[j] + d2 / omega_P * j / (j + 1) * M_inv * p_array[j]
        p_array[j] = -omega_P * d2 * (j + 1) / j * M * s_array[j] + np.cos(
            omega_P * dt * 0.5) * p_array[j]
        s_array[j] = s_j
    # 最后演化半步控温
    p_array[
        0] = o1 * p_array[0] + o2 * np.sqrt(M) * np.random.standard_normal(N)
    for j in range(1, P):
        p_array[j] = o1 * p_array[j] + o2 * np.sqrt(
            (j + 1) / j * M) * np.random.standard_normal(N)

    return s_array, p_array


def temperature(p_array):
    """
    用于计算当前构型的温度
    p_array: numpy 2-D array 当前轨迹的动量
    retrun: float 当前构型的温度
    """
    T = 0
    T += p_array[0].dot(M_inv * p_array[0])
    for j in range(1, P):
        T += p_array[j].dot(j / (j + 1) * M_inv * p_array[j])
    return T / (N * P)


def main():
    """
    主程序
    """
    # 定义轨线初始位置，按照标准正态分布随机在相空间内取点
    x_init = np.array([np.random.standard_normal(N) for j in range(P)])
    p_init = np.array([np.random.standard_normal(N) for j in range(P)])
    s_init = staging_transf(x_init)
    init = (s_init, p_init)
    # 定义储存对应构型动能、势能、温度的变量
    potential_array = []
    temperature_array = []
    kinetic_array = []

    print("the number of beads = %d" % P)
    print("the temperature you set = %.2f" % (1 / beta))
    print("total steps = %d" % steps)
    print("---------------START---------------")
    start = time.time()
    for i in range(steps):
        init = BAOAB(init)
        potential_array.append(phi_potential(init[0]))
        temperature_array.append(temperature(init[1]))
        kinetic_array.append(kinetic_energy_viral(staging_transf_inv(init[0])))
        if i % 5000 == 0:
            print("step = %d, time = %.6f" % (i, time.time() - start))
    T = sum(temperature_array) / steps
    KE = sum(kinetic_array) / steps
    V = sum(potential_array) / steps
    print("the temperature of the trajectory is %.4f,"
          "the setting temperature is %.4f" % (T, 1 / beta))
    print("the kinetic energy is %.4f" % KE)
    print("the potential energy is %.4f" % V)
    print("----------------END----------------")


if __name__ == "__main__":
    main()
