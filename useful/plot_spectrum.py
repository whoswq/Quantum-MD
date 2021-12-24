import numpy as np
import matplotlib.pyplot as plt
from numpy import fft

data_1 = np.loadtxt("/home/wcb/liugroup/my attemps/PILD/water_300_LSCIVR_cstr_p_f_cstr.txt")

t = data_1.T[0]
p_corr = data_1.T[1]

N = len(t)
dt = t[1] - t[0]
factor = 6.57968e15 / 2.998e8 / 100
omega = np.array([2 * np.pi * k / (N * dt) for k in range(N)]) * factor
amp = np.abs(fft.fft(p_corr - sum(p_corr) / len(p_corr))) / N  # 归一化的振幅谱
plt.plot(omega[0:N // 30], amp[0:N // 30], label="LSC-IVR")
plt.vlines(1594.78, 0, 0.7, color="red", label="exp data")
plt.vlines(3657.04, 0, 0.7, color="red")
plt.vlines(3755.96, 0, 0.7, color="red")
plt.legend()
plt.savefig("/home/wcb/liugroup/my attemps/PILD/water_300_LSCIVR_cstr_p_f_cstr.pdf")
plt.show()
