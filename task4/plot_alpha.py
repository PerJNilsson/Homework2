# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'beta_alpha.dat'

# import data
data = np.loadtxt(filename)


# initial size of plot window
plt.figure(figsize=(8,6))

x = np.linspace(0, len(data), len(data))
# plot
plt.plot(data[:,0],'-', label=r'$\beta=$'+str(0.51))
for i in range(1,6):
    plt.plot(data[:,i],'-', label=r'$\beta=$'+str(0.5+i/10))


# labels
plt.xlabel('# gradient steps', fontsize=20)
plt.ylabel(r'$\alpha$', fontsize=20)

# legend
plt.legend(loc='upper right')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits)
plt.ylim(0.05, 0.25)
# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Convergence of ' + r'$\alpha$' + ' for different '+r'$\beta$' '')
plt.savefig('4_alpha_convergence.png')
# display the plot


plt.show()
