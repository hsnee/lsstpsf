import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import seaborn as sns;sns.set_style('darkgrid')

data0 = np.loadtxt('0 corr_arrays.txt')
data1 = np.loadtxt('1 corr_arrays.txt')

ra0,dec0,g10,g20,size0 = data0[0],data0[1],data0[2],data0[3],data0[4]
ra1,dec1,g11,g21,size1 = data1[0],data1[1],data1[2],data1[3],data1[4]

size_avg = 0.5*(size0 + size1)
delta_size = size1-size0

# Now make a histogram

plt.figure()
plt.hist(delta_size,label=r'$\frac{\Delta \sigma}{< \sigma >}$')

# plt.xlabel(r'\')
# plt.ylabel(r'$\xi$')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0.)
plt.savefig('delta_size_histogram.pdf')
plt.close()
