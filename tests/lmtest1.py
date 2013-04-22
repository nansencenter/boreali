#!/usr/bin/python
from lmtest0 import *

# select silicon sand albedo
albedo = albedos[0](wavelen)

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
        'h_on_depth_tsw.png',
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
        r = lm.get_rrsw(parameters, model, cComb[1], albedo, h, theta, 6)[1]
        plt.plot(wavelen, r, 'o-')
        legendVals.append('%02d m' % h)
    plt.legend(legendVals)
    plt.xlabel('wavelength, nm')
    plt.ylabel('Rrsw, sr-1')
    plt.title(cComb[2])
    plt.savefig(cComb[0])
    plt.close()
    
