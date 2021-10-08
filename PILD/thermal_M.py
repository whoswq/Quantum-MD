"""
calculate the M_thermal
"""
import numpy as np
import numpy.linalg as liag

N = 3  # the freedom of the system
P = 200  # the number of beads
steps = 200000  # total steps of evolution
hbar = 1  # reduced Planck constant
beta = 0.5  # the inverse of temperature
omega_P = np.sqrt(P) / (beta * hbar)  # the frequency of internal force
dt = 0.1 / omega_P  # time step
M = np.array([])  # mass matrix
M_inv = np.array([])
# params of morse potential
De = 20
alpha = 2.0
r_eq = 2.5
# 选定一个标准的平衡构型


def Morse_potential(x):
    """
    Morse potential at x
    x: 1d array length=N coordinate
    return: double
    """
    # Morse势的参数
    q = 1 - np.exp(-alpha * (liag.norm(x) - r_eq))
    return De * q * q


def nabla_Morse_potential(x):
    """
    the gradient of Morse potential at x
    x: 1d array length=N coordinate
    return: 1d array length=N
    """
    r = np.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    exp_a = np.exp(-alpha * (r - r_eq))
    fractor = 2 * De * alpha * (1 - exp_a) * exp_a / r
    return np.array([fractor * x[0], fractor * x[1], fractor * x[2]])


def Hessian_Morse_potential(x, dx=0.01):
    """
    the hessian matrix of Morse potential at x
    由于计算导数的解析表达式过于复杂，考虑使用中心差分计算高阶导数
    x: 1d array length=N coordinate
    return 2d array N*N
    """
    a0 = nabla_Morse_potential((x[0] + dx, x[1], x[2]))
    b0 = nabla_Morse_potential((x[0] - dx, x[1], x[2]))
    a1 = nabla_Morse_potential((x[0], x[1] + dx, x[2]))
    b1 = nabla_Morse_potential((x[0], x[1] - dx, x[2]))
    a2 = nabla_Morse_potential((x[0], x[1], x[2] + dx))
    b2 = nabla_Morse_potential((x[0], x[1], x[2] - dx))
    return np.array([
        [(a0[i] - b0[i]) / (2*dx) for i in range(3)],
        [(a1[i] - b1[i]) / (2*dx) for i in range(3)],
        [(a2[i] - b2[i]) / (2*dx) for i in range(3)]
    ])


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


def Harmonic_potential(x):
    """
    x: numpy 1-D array Cartesian坐标下单个bead的坐标
    return: float 单个bead的势能
    """
    return 0.5 * x.dot(x)


def nabla_Harmonic_potential(x):
    """
    x: numpy 1-D array Cartesian坐标下单个bead的坐标
    return: numpy 1-D array 受力
    """
    return x


def phi_potential(s_array, potential):
    """
    x_array: numpy 2-D array staging坐标下所有beads的位置
    potential: the potential energy of the system
    return: float 势能
    """
    x_array = staging_transf_inv(s_array)
    phi = 0
    for j in range(P):
        phi += potential(x_array[j])
    return phi / P


def nabla_phi_potential(s_array, nabla_potential):
    """
    s_array: numpy 2-D array 所有beads的staging坐标
    return: numpy 2-D array 所有beads在staging坐标下的受力
    """
    x_array = staging_transf_inv(s_array)
    s_force_array = np.zeros((P, N))
    force0 = nabla_potential(x_array[0])
    # 构造第一个bead的受力
    for j in range(1, P):
        force0 += nabla_potential(x_array[j])
    s_force_array[0] += 1 / P * force0
    # 递推其他bead的受力
    for j in range(1, P):
        s_force_array[j] += 1 / P * nabla_potential(
            x_array[j]) + (j - 1) / j * s_force_array[j - 1]
    return s_force_array


def kinetic_energy_viral(x_array, nabla_potential):
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
        K += (x_array[j] - x_c).dot(nabla_potential(x_array[j]))
    return N / (2 * beta) + K / (2 * P)


# 记录一些演化时会用到的变量
c1 = np.exp(-omega_P * dt)
c2 = np.sqrt((1 - c1 * c1) / beta)
d1 = np.cos(omega_P * dt * 0.5)
d2 = np.sin(omega_P * dt * 0.5)


def BAOAB(s_array, p_array, nabla_potential):
    """
    使用BAOAB的算法演化Lagevin方程
    pre: 前一个时刻的位置和动量 staging坐标下
    return: 下一个时刻的位置和动量 staging坐标下
    """
    # 首先演化半步动量
    p_array = p_array - \
        nabla_phi_potential(s_array, nabla_potential) * dt * 0.5
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
    p_array = p_array - \
        nabla_phi_potential(s_array, nabla_potential) * dt * 0.5

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


def M_thermal_digonal_estimator_p(x_0, T_0):
    """
    the diagonal term of M_thermal(ensemble average)
    x_0: 1-d array length=N
    T_0: 2-d array N*N the eigenvector of Hessian matrix
    """

    pass


def M_thermal_estimator_p(x_array, T_0):
    """
    the ensemble average of M_thermal
    x_0: 1-d array length=N
    T_0: 2-d array N*N the eigenvector of Hessian matrix
    """
    pass


