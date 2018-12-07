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

plt.hist(data, bins=50, density=True, alpha=0.7, rwidth=0.85, label=r'$\theta$')

# labels
plt.xlabel('Radians', fontsize=20)
plt.ylabel(r'$acos(\frac{r_1r_2}{|r_1||r_2|})$', fontsize=20)

# legend
plt.legend(loc='upper left')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits
# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Historgram of '+ r'$acos(\frac{r_1r_2}{|r_1||r_2|})$')

# display the plot

"""
timestep = 1E-7
radius = (2.79*1E-6)/2.0;
density = 2.65 * 1E3;
m = 4.0*math.pi*radius*radius*radius *density*3.0; # assuming particle is sphere-ish
temperature = 297;
tau = 48.5*1E-6; # 48.5 | 147.3;
eta = 1.0/tau;
c_0 = math.exp(-2*eta*timestep);
c_1 = math.exp(-eta*timestep);

k_b = 1.380*1E-23;
v0 = 2*1E-3

const_part = m /(k_b*temperature*(1.0-c_0))
exp_part = -m*((0-v0*c_1)**2)/(2.0*k_b*temperature*(1.0-c_0))
f = (math.sqrt(const_part)*math.exp(exp_part))

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
