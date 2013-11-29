#!/usr/bin/python
from lmtest13_setup import *

print 'Plotting...'
rmsenp = np.load(v0Name + 'rmse.npy')

for ra in [0, 1]:
    for var in range(3):
        plt.close()
        legvals = []
        for v1 in range(len(v1vector)):
            for v2 in range(len(v2vector)):
                plt.plot(v0vector, rmsenp[ra, :, v1, v2, var], v1Lines[v1] + v2Lines[v2]);
                legvals.append(v1Titles[v1] + v2Titles[v2])
    
        plt.ylabel(varTitles[var] + ylabel[ra])
        plt.xlabel(xlabel)
        
        plt.legend(legvals)
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1, x2 + x2 * 0.4, y1, y2))
        plt.savefig(v0Name + varTitles[var] + raTitles[ra] + '.png')
        plt.close()
    
