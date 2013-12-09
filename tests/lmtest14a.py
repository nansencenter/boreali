#!/usr/bin/python
from boreali import lm, Boreali

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm

cmap = cm.get_cmap('jet')


# PLOT INSITU SPECTRA
a = np.loadtxt('conc_albe.txt')

ylabels = [
'Ed(-0), w m-2',
'Ed(-1), w m-2',
'Lu(-0), w m-2 sr-1',
'Lu(-0), w m-2 sr-1',
'Rrsw(-0), sr-1',
'Rrsw(-1), sr-1',
'Ed(-0) - Ed(-1), w m-2',
'Rrsw(-0) - Rrsw(-1), w m-2',
]

fnames = [
'ed0',
'ed1',
'lu0',
'lu1',
'rr0',
'rr1',
'ed01',
'rr01',
]

for var in range(8):
    leg = []
    for i in [0,1,2,3,4,5,6,7,8]:
        if var < 6:
            vals = a[:, 1 + i * 6 + var]
        elif var == 6:
            vals = a[:, 1 + i * 6 + 0] - a[:, 1 + i * 6 + 1]
        elif var == 7:
            vals = a[:, 1 + i * 6 + 4] - a[:, 1 + i * 6 + 5]

        plt.plot(a[:, 0], vals, '-', c=cmap(i/8.))
        leg.append(i+1)
    plt.legend(leg)
    plt.xlabel('wavelength, nm')
    plt.ylabel(ylabels[var])
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1, x2 + x2 * 0.05, y1, y2))
    plt.savefig(fnames[var])
    plt.close()
    

