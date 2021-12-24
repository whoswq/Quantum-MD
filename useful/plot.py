import numpy as np
import matplotlib.pyplot as plt
from numpy import fft

data_1 = np.loadtxt("./1-D_quartic_0-1_average_compare.txt")
data_2 = np.loadtxt("./1-D_quartic_0-1_LSCIVR_p.txt")
t_1 = data_1.T[0]
p_corr_1 = data_1.T[1]
t_2 = data_2.T[0]
p_corr_2 = data_2.T[1]
"""
N = len(t)
dt = t[1] - t[0]
omega = np.array([2 * np.pi * k / (N * dt) for k in range(N)])
amp = np.abs(fft.fft(p_corr)) / N  # 归一化的振幅谱
plt.plot(omega[0:N//2], amp[0:N//2])"""
plt.plot(t_1, p_corr_1, label="PILD")
plt.plot(t_2[0:-1000], p_corr_2[0:-1000], label="LSC-IVR")
plt.title("quartic poential, $\\beta = 0.1$")
plt.ylabel("$\\left<p(0)p(t)\\right>_{\\mathrm{Kubo}}$")
plt.xlabel("time")
plt.legend()
plt.savefig("1-D_quartic_0-1_compare.png", dpi=800)
plt.show()
