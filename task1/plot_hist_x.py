# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np
import math
import matplotlib.mlab as mlab
from scipy.stats import norm
# input file
filename = 'theta_dist.dat'
# import data
data = np.loadtxt(filename)
# initial size of plot window
plt.figure(figsize=(8,6))
h = []
# plot


# labels
plt.xlabel(r'$x = cos(\theta)$', fontsize=20)
plt.ylabel('P(x)', fontsize=20)

x_vec = []
prob_fun = []

for i in range(0, len(data)):
    x_vec.append(-1.0+i*2.0/ len(data))
    prob_fun.append(np.cos(data[i]))

#plt.plot(x_vec, prob_fun, linewidth = 2, label=r'$0.5*sin(\theta)$')

plt.hist(prob_fun, bins=50, density=True, alpha=0.7, rwidth=0.85)



# legend
#plt.legend(loc='upper right')
#leg = plt.gca().get_legend()
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize=12)


# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Histogram of P(x), where ' + r'$x=cos(\theta)$')
plt.savefig('1_px.png')
# display the plot

plt.show()
