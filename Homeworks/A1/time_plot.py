import numpy as np
import matplotlib.pyplot as plt

m_list = np.array([10, 100, 1000, 1000, 10000])
n_list = np.array([10, 100, 1000, 1000, 10000])

mn = m_list*n_list

times = [0.000118017196655, 0.00748491287231, 0.849985122681, 0.971807956696, 115.979768038]

x = np.array([10, 10000])
x = x*x
first_order = x*1e-4
second_order = x*x*1e-7

plt.figure()
plt.loglog(mn, times, '-o',label='timing')
plt.plot(x,first_order,label='first order')
plt.plot(x,second_order,label='second order')
plt.xlabel('m*n')
plt.ylabel('running time')
plt.legend()
plt.savefig('timeing.png')
plt.show()
