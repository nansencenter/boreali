#!/usr/bin/python
from boreali import lm, Boreali

import numpy as np
import matplotlib.pyplot as plt


wavelen = [412, 443, 490, 510, 555, 670]
albedo = [0, 0, 0, 0, 0, 0]
h = -10
theta = 0

parameters=[len(wavelen), 3, 1e-6, 
            0.01, 2,
            0.01, 1,
            0.01, 1,]

b = Boreali('ladoga', 'ladoga', wavelen)

model = b.get_homodel()

chls = [0.01, 0.05, 0.1, 0.5, 1, 5]
legendVals = []
for chl in chls:
    r = lm.get_rrsw(parameters, model, [chl, 0.01, 0.01], albedo, h, theta, len(wavelen))[1]
    plt.plot(wavelen, r, '.-')
    legendVals.append('%5.2f mg m-3' % chl)

plt.legend(legendVals)
plt.xlabel('wavelength, nm')
plt.ylabel('Rrsw, sr-1')
plt.title('Rrsw spectra for various chl values')
plt.savefig('rrsw_chl.png')
plt.close()
