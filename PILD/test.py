import numpy as np

m_H = np.array([1837, 1837, 1837])
m_O = np.array([1823, 1823, 1823]) * 16
M = m_H * m_O / (m_H + m_O)
M_inv = 1 / M
M = np.diag(M)
M_inv = np.diag(M_inv)
M_t = np.array([[48006.248658, 119.025604, -164.144211],
                [119.025604, 2458.158386, -0.675177],
                [-164.144211, -0.675177, 2458.392372]])

M_t_inv = np.linalg.inv(M_t)

print(np.linalg.eigh(M_t_inv)[0])
