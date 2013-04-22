#!/usr/bin/python
import lm

import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate
# Interpolate chara/cladophora/silt/sand albedos with one variable
plt.close()
a = np.loadtxt('albedos.txt').T
b = np.zeros(a.shape)
pp = np.zeros((a.shape[0], 5))
for i in range(len(a)-1):
    p = np.polyfit(a[-1], a[i], 4)
    pp[i] = p[:]
    y = np.polyval(p, a[-1])
    plt.plot(a[i], y, '.')
    b[i] = y
    print p, ((y - a[i])**2).sum()

plt.show()
plt.plot(b[:-1], 'o-');plt.show()
plt.plot(pp[:-1], 'o-');plt.show()
