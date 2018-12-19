# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np
#plt.style.use('seaborn-whitegrid')

# input file
filename = 'alpha_energies.dat'

# import data
data = np.loadtxt(filename)
cols = np.size(data,1)

# initial size of plot window
plt.figure(figsize=(8.5,6.2))

x = np.linspace(0, len(data), len(data))
# plot
xerr = 0
plt.errorbar(data[:,0], data[:,1], yerr=2*data[:,2],fmt='o', ecolor='r', capthick=2, label='error bars = 2'+r'$\sigma$')


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
plt.title('Average energy of the Helium atom with respect to alpha with error bars at two standard deviations')
plt.savefig('3_alpha_min.png')
# display the plot


plt.show()
