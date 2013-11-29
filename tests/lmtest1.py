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

# plot Rrsw / depth dependence
cCombinations = [
    [
        'h_on_depth_wat.png',
        [0.01, 0.01, 0.01],
        'Pure water: CHL=0.01 ug l-1, TSM=0.01 mg l-1, DOC=0.01 mgC l-1',
    ],
    [
        'h_on_depth_chl.png',
        [1., 0.01, 0.01],
        'CHL dominated: CHL=1 ug l-1, TSM=0.01 mg l-1, DOC=0.01 mgC l-1',
    ],
    [
        'h_on_depth_tsm.png',
        [0.01, 0.5, 0.01],
        'TSM dominated: CHL=0.01 ug l-1, TSM=0.5 mg l-1, DOC=0.01 mgC l-1',
    ],
    [
        'h_on_depth_doc.png',
        [0.01, 0.01, 0.5],
        'DOC dominated: CHL=0.01 ug l-1, TSM=0.01 mg l-1, DOC=0.5 mgC l-1',
    ],
]

for cComb in cCombinations:
    plt.close()
    legendVals = []
    for h in range(1, 30, 4):
        r = lm.get_rrsw_shal(model, cComb[1], theta, h, albedo, len(wavelen))[1]
        plt.plot(wavelen, r, 'o-')
        legendVals.append('%02d m' % h)
    plt.legend(legendVals)
    plt.xlabel('wavelength, nm')
    plt.ylabel('Rrsw, sr-1')
    plt.title(cComb[2])
    plt.savefig(cComb[0])
    plt.close()
    
