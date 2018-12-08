# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'alpha_energies.dat'

# import data
data = np.loadtxt(filename)
cols = np.size(data,1)

# initial size of plot window
plt.figure(figsize=(8,6))

x = np.linspace(0, len(data), len(data))
# plot

plt.plot(data[:,0], data[:,1],'b-')


# labels
plt.xlabel('alpha', fontsize=20)
plt.ylabel('Energy / [a.u.]', fontsize=20)

# legend
plt.legend(loc='upper right')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits
# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Energy of Helium atom with respect to alpha')

# display the plot


plt.show()
