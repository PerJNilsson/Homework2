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

"""
(mu4, sigma4)  = norm.fit(data[:,3])
(mu3, sigma3)  = norm.fit(data[:,2])
(mu2, sigma2)  = norm.fit(data[:,1])
(mu1, sigma1)  = norm.fit(data[:,0])

x = np.linspace(mu1 - 3*sigma1, mu1 + 3*sigma1, 100)
plt.plot(x,mlab.normpdf(x, mu1, sigma1), '--', color='Black', linewidth=3)
x1 = np.linspace(mu2 - 3*sigma2, mu2 + 3*sigma2, 100)
plt.plot(x1,mlab.normpdf(x1, mu2, sigma2), '--', color='Black', linewidth=3)
x2 = np.linspace(mu3 - 3*sigma3, mu3 + 3*sigma3, 100)
plt.plot(x2,mlab.normpdf(x2, mu3, sigma3), '--', color='Black', linewidth=3)
x3 = np.linspace(mu4 - 3*sigma4, mu4 + 3*sigma4, 100)
plt.plot(x3,mlab.normpdf(x3, mu4, sigma4), '--', color='Black', linewidth=3)
"""
plt.show()
