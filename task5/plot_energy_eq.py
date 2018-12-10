# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'energy_data_eq.dat'

# import data
data = np.loadtxt(filename)
cols = np.size(data,1)

# initial size of plot window
plt.figure(figsize=(8,6))

# plot
plt.plot(data[:,0], data[:,1],'-')


# labels
plt.xlabel('Iterations', fontsize=20)
plt.ylabel('Energy$_{local}$ / [a.u.]', fontsize=20)

# legend
#plt.legend(loc='upper right')
#leg = plt.gca().get_legend()
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize=12)

# axis limits
# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Energy plot from system with an initial configuration of '+ r'$\vec r_1 = \vec r_2=[50,50,50]$'+' a.u.')
plt.savefig('2_equilibration.png')
# display the plot


plt.show()
