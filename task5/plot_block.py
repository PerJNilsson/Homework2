# plot the energies
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np
import math

# input file
filename = 'block_value.dat'

# import data
data = np.loadtxt(filename)

# initial size of plot window
plt.figure(figsize=(8,6))


sum = 0
start_avg = int(len(data[:,1]) / 5); 
for i in range(start_avg, len(data[:,1])):
    sum = sum + data[i,1] / (len(data[:,1])-start_avg)
g = float("{0:.2f}".format(sum))
sum2 = []
for i in range(0, len(data[:,0])):
    sum2.append(sum)

# plot
plt.plot(data[:,0], data[:,1],'X')
plt.plot(data[:,0], sum2,'-r', label=r'$s=\lim_{B\rightarrow large}\frac{B*Var(F)}{Var(f)}$='+str(g))

# labels
plt.xlabel('Blocks', fontsize=20)
plt.ylabel('Statistical inefficiency', fontsize=20)
plt.title('Statistical inefficiency='+r'$B\frac{Var(E_{L,block})}{Var(E_L)}$')
# legend
plt.legend(loc='lower right')
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
plt.savefig('block.png')

plt.show()
