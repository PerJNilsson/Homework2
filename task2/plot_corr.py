1# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'values.dat'

# import data
data = np.loadtxt(filename)

# initial size of plot window
plt.figure(figsize=(8,6))

# plot
dasploot=plt.plot(data[:,0], data[:,1],'-')
plt.plot(10.40, np.exp(-2), 'rx', markersize=15, label=r'$\phi=e^{-2}$'+' at k=10.40')

# labels
plt.xlabel('k', fontsize=20)
plt.ylabel('$\phi$(k)', fontsize=20)
plt.title('Plot of autocorrelation function with respect to the lag k')
# legend
plt.legend(loc='upper right')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12)

# axis limits
#plt.xlim([0,data[-1,0]])
#plt.ylim([min(data[:,3])-0.002,max(data[:,1])+0.002])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# display the plot
plt.savefig('afc.png')
plt.show()
