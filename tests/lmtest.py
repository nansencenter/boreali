#!/usr/bin/python
from boreali import lm, Boreali

import numpy as np
import matplotlib.pyplot as plt

#from boreali import Boreali
from scipy import interpolate

parameters=[0.01, 5.,
            0.01, 2.,
            0.01, .5,]

#wavelen = [412, 443, 490, 510, 555, 670]
wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709]
b = Boreali('michigan', wavelen)
# get matrix with HO-model
model = b.get_homodel()

# get sand albedo
albedo = b.get_albedo([0])[0]

print 'MODEL:', model
print 'ALBEDO:', albedo

# deep
r0 = lm.get_rrsw_deep(model, [0.3, 0.2, 0.1], 25, len(wavelen))[1]
print 'R0: ', r0
c0 = lm.get_c_deep(parameters, model, [r0], [25], 4)[1]
print 'C0: ', c0

# shallow
albedo = [0.043,0.053,0.096,0.134,0.204,0.128,0.109,0.111,0.117]
print 'ALBEDO:', albedo
r1 = lm.get_rrsw_shal(model, [0.3, 0.2, 0.1], 25, 5, albedo, len(wavelen))[1]
print 'R1: ', r1
c1 = lm.get_c_shal(parameters, model, [r1], [25], [5], [albedo], 4)[1]
print 'C1:', c1

# shallow with albedo parametrization
r2 = lm.get_rrsw_albe(model, [0.3, 0.2, 0.1, 0.4, 0.4], 25, 5, wavelen, len(wavelen))[1]
print 'R2: ', r2
c2 = lm.get_c_albe(parameters, model, [r2], [25], [5], wavelen, 6)[1]
print 'C2:', c2
    
#plt.plot(wavelen, r0, 'o-', wavelen, r1, 'o-', );plt.show()
