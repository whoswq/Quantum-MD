import numpy as np
import matplotlib.pyplot as plt
from numpy import fft

data = np.loadtxt("./3-D_Harmonic.txt")

t = data.T[0]
p_corr = data.T[1]
N = len(t)
dt = t[1] - t[0]
omega = np.array([2 * np.pi * k / (N * dt) for k in range(N)])
amp = np.abs(fft.fft(p_corr)) / N  # 归一化的振幅谱
plt.plot(omega[0:N//2], amp[0:N//2])
plt.show()
