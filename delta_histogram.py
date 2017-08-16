import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import seaborn as sns;sns.set_style('darkgrid')

data0 = np.loadtxt('0 corr_arrays.txt')
data1 = np.loadtxt('1 corr_arrays.txt')
data2 = np.loadtxt('2 corr_arrays.txt')
ra0,dec0,g10,g20 = data0[0],data0[1],data0[2],data0[3]
ra1,dec1,g11,g21 = data1[0],data1[1],data1[2],data1[3]
ra2,dec2,g12,g22 = data2[0],data2[1],data2[2],data2[3]
g1_avg = 0.5*(g10 + g11)
g2_avg = 0.5*(g20 + g21)
delta_g1 = g11-g10
delta_g2 = g21-g20

# Now make a histogram

plt.figure()
plt.histogram(delta_g1,label=r'$\frac{\detla g_1}{<g1>}$')
plt.histogram(delta_g2,label=r'$\frac{\delta g_2$}{<g2>}')
plt.xlabel('g')
#plt.ylabel(r'$\xi$')
plt.xscale('log')
plt.yscale('log')
plt.legend([r'$\xi_+$',r'$\xi_-$'],bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0.)
plt.savefig('fractional_delta_shear.pdf')
plt.close()
