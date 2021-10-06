氖原子之间的势能采用Lennard-Jones势：
$$
V_{LJ}(r) := 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]
$$
其中参数为：
$$
\sigma = 2.75A,\, \frac{\epsilon}{k_{B}} = 35.8K
$$
平衡位置为：$r_{eq} = 2^{1/6}\sigma$

考虑经典MD，为了处理方便，将方程转化为无量纲的方程（这是数值计算时常用的一个策略，事实上**原子单位制**就是由氢原子的Schrodinger方程无量纲化而来）。为了简单起见，考虑如下Hamilton量：
$$
H(r, p) = \frac{p^2}{2m} + V_{LJ}(r)
$$
令$H = \epsilon H_0,\,p = \beta_{p}p_0,\,m = \beta_m m_0,\,r = \sigma r_0$,，其中下标为0的代表没有量纲的常数，量纲全部吸收进系数中，将它们带入Hamilton量：
$$
H_0 = \frac{\beta_p p_0^2}{2\beta_m\epsilon m_0} + 4\left[\frac{1}{r_0^{12}} - \frac{1}{r_0^{6}}\right]
$$
如果为了保持Hamilton量的形式，可以使：
$$
\beta_p^2 = \beta_m \epsilon
$$
同时对Hamilton正则方程进行无量纲化，希望Hamilton正则方程的形式不变：
$$
\frac{\sigma\mathrm{d}r_0}{\beta_t\mathrm{d}t_0} = \frac{\epsilon\partial H_0}{\beta_p \partial p_0}
$$
（对另一个方程无量纲化可以得到相同的结果）这就要求：
$$
\sigma \beta_p = \epsilon\beta_t
$$
联立方程(5)(7)可以得到：
$$
\beta_t^2 = \frac{\beta_m\sigma^2}{\epsilon}
$$


这时只要给$\beta_m$一个合适的定义就可以使得在Lennard-Jones势下的Hamilton方程转化为无量纲的形式，可以有两种选择，第一选择氖原子的质量为单位质量即$m_{neon} = \beta_m$，这时（将国际单位制下对应带量纲系数的真实数值代入，可以得到在我们定义的这套“单位制”下时间与国际单位制中时间的换算关系）:
$$
\beta_t = \sqrt{\frac{m_{neon}\sigma^2}{\epsilon}} =
2.2643\times 10^{-12}\mathrm{s}
$$
也可以将原子量的单位道尔顿选为标准质量，得到另一套单位制。

如果要应用Langevin动力学，我们还必须讨论$k_B$的单位换算，但是考虑到$k_B T$具有能量量纲，直接选用$\epsilon$为单位表示即可

第一次测试选取$T = 40K,\,\rho\sigma^3 = 0.68$