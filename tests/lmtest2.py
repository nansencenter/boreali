#!/usr/bin/python
from boreali import lm, Boreali

import numpy as np
import matplotlib.pyplot as plt

import time

parameters=[0.01, 5,
            0.01, .5,
            0.01, .2,]

wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709]


b = Boreali('michigan', wavelen)
# get matrix with HO-model
model = b.get_homodel()
albedo = b.get_albedo([0])[0]

theta = 25

# plot Rrsw / C dependence for semi shallow water
h = 5

cCombinations = [
    [
        np.arange(0.01, 5, 0.8),
        'H=5 m, TSM=0.01 mg l-1, DOC=0.01 mgC l-1',
        'CHL: %4.2f',
        'h_on_c_chl.png',
    ],
    [

        np.arange(0.01, 4, 0.6),
        'H=5 m, CHL=0.01 mg l-1, DOC=0.01 mgC l-1',
        'TSM: %4.2f',
        'h_on_c_tsm.png',
    ],
    [
        np.arange(0.01, 4, 0.6),
        'H=5 m, CHL=0.01, TSM=0.01',
        'DOC: %4.2f',
        'h_on_c_doc.png',
    ],
    [
        np.arange(0.01, 4, 0.6),
        'H=5 m',
        'CHL=TSM=DOC: %4.2f',
        'h_on_c_all.png',
    ]
]

for ci in range(3):
    cComb = cCombinations[ci]
    cVec = cComb[0]
    titleString = cComb[1]
    lableString = cComb[2]
    ofilename = cComb[3]
    
    c = [0.01, 0.01, 0.01]
    plt.close()
    legendVals = []
    rrr = []
    for c0 in cVec:
        c[ci] = c0
        #c[1] = c0
        #c[2] = c0
        #r = lm.get_rrsw(model, c, albedo, h, theta, len(wavelen))[1]
        r = lm.get_rrsw_shal(model, c, theta, h, albedo, len(wavelen))[1]
        plt.plot(wavelen, r, 'o-')
        legendVals.append(lableString % c0)
        rrr.append(r)
        
    plt.legend(legendVals)
    plt.xlabel('wavelength, nm')
    plt.ylabel('Rrsw, sr-1')
    plt.title(titleString)
    plt.savefig(ofilename)
    
