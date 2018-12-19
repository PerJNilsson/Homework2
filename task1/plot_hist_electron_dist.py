# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np
import math
import matplotlib.mlab as mlab

# input file
filename = 'electron_dist.dat'
# import data
data = np.loadtxt(filename)
# initial size of plot window
plt.figure(figsize=(8,6))
#h = []

# plot histogram, data[:,0]=electron1, data[:1]=electron2
plt.hist(data[:,0], bins=100, density=True, alpha=0.7, rwidth=0.85, label='simulated distances')

# plot theoretical curve
nbr_points = 100
r = np.linspace(0, 4, nbr_points)
p = np.zeros(nbr_points)

Z = 2
for i in range(0,nbr_points):
    p[i] = Z**3*4*r[i]*r[i]*np.exp(-2*Z*r[i])
plt.plot(r, p, '--', color='grey',linewidth=2, label=r'$\rho(r)=Z^34r^2e^{-2Zr}$' + ' with Z = 2')

Z = 27/16
for i in range(0,nbr_points):
    p[i] = Z**3*4*r[i]*r[i]*np.exp(-2*Z*r[i])
plt.plot(r, p, '--', color='black',linewidth=2, label=r'$\rho(r)=Z^34r^2e^{-2Zr}$' + ' with Z = 27/16')


# labels
plt.xlabel('distance [a.u.]', fontsize=20)
plt.ylabel('Probability', fontsize=20)

# legend
plt.legend(loc='upper right')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12)

# axis limits
# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Historgram of electron distance from nucleus')
plt.savefig('1_electron_distance.png')
# display the plot

plt.show()
