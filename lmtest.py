#!/usr/bin/python
import time

from boreali import lm, Boreali

import numpy as np
import matplotlib.pyplot as plt

#from boreali import Boreali
from scipy import interpolate

parameters=[0.01, 15,
            0.01, 5,
            0.01, 5,]


bands = 6
wavelen = [412, 443, 490, 510, 555, 670]
b = Boreali('michigan', 'michigan')
# get matrix with HO-model
model = []
for i in range(0, 8):
    model.append(b.homo[i](wavelen))

#sun at typical zenith angle
theta = 25

#tough test of LM-retrieval
#generate N random concentrations, solar zenith angles
pixels = 5000
cRange = np.array([15., 5., 10.])
ccc = np.multiply(cRange, np.random.rand(pixels, 3))
ttt = np.random.rand(pixels) * 10

# calculate N Rrsw spectra
rrr = []
for i, c in enumerate(ccc):
    r = lm.get_rrsw(model, c, ttt[i], 6)[1]
    # add random noise to Rrsw +/- %
    rNoise = 1 + (np.random.randn(1, r.shape[0]) - np.random.randn(1, r.shape[0])) * 0.01
    r *= rNoise[0]
    # append r to list
    rrr.append(r)

# retrieve N concetrations
t0 = time.time()
ccc2 = lm.get_c(parameters, model, rrr, ttt, 4*pixels)[1]
t1 = time.time()
print 'Time spent:', t1-t0
ccc2 = np.array(ccc2).reshape(pixels, 4)

# get relative difference and keep it
dccc = (ccc - ccc2[:, 0:3]) / ccc

# plot comparison of concentrations
cTitles = ['CHL-CHL2', 'TSM-TSM2', 'DOC-DOC2']
for ci in [0, 1, 2]:
    plt.close()
    # get relative difference for that depth for that concentration
    dc = dccc[:, ci]
    # add histogram
    plt.hist(dc, np.arange(-0.5, 0.5, 0.02), alpha=0.5);
    print len(dc[np.abs(dc) < 0.1]) / float(pixels)

    plt.title(cTitles[ci])
    plt.xlabel('concentration')
    plt.ylabel('pixels')
    plt.savefig('hist_%s_a.png' % cTitles[ci])
    
