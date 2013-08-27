#!/usr/bin/python
from lmtest0 import *
import time

# select silicon sand albedo
albedo = albedos[0](wavelen)

dcccH = []

for h0 in [5, 10, 20]:
    #tough test of LM-retrieval
    #generate N random concentrations, depths, albedos, solarzenith angles
    pixels = 5000
    cRange = np.array([1., .5, .5])
    ccc = np.multiply(cRange, np.random.rand(pixels, 3))
    #hhh = 5 + np.random.rand(pixels) * 5  # H in range 5 - 10 m
    hhh = h0 * np.ones(pixels)
    aaa = np.repeat(np.array([albedo]), pixels, axis=0)
    ttt = np.random.rand(pixels) * 10
    
    # calculate N Rrsw spectra
    rrr = []
    for i, c in enumerate(ccc):
        r = lm.get_rrsw(parameters, model, c, aaa[i], hhh[i], ttt[i], 6)[1]
        rrr.append(r)
    
    # add random noise to depth (+/- 1 m)
    hhh2 = hhh - 1 + np.random.rand(pixels) * 2
    dh = (hhh - hhh2) ** 2
    
    # add random noise to bottom type
    albedo2 = albedos[1](wavelen)
    aaa2 = np.repeat(np.array([albedo2]), pixels, axis=0)
    aaa2ratio = np.random.rand(pixels, 1) * 0.3
    aaa2 = aaa - np.multiply(aaa, aaa2ratio) + np.multiply(aaa2, aaa2ratio)
    
    # retrieve N concetrations
    t0 = time.time()
    ccc2 = lm.get_c(parameters, model, rrr, aaa2, hhh, ttt, 4*pixels)[1]
    t1 = time.time()
    print 'Time spent:', t1-t0
    ccc2 = np.array(ccc2).reshape(pixels, 4)
    
    # get relative difference and keep it
    dccc = (ccc - ccc2[:, 0:3]) / ccc
    dcccH.append(dccc)

# plot comparison of concentrations
cTitles = ['CHL-CHL2', 'TSM-TSM2', 'DOC-DOC2']
for ci in [0, 1, 2]:
    plt.close()
    for hi in [0, 1, 2]:
        # get relative difference for that depth for that concentration
        dch = dcccH[hi][:, ci]
        # add histogram
        plt.hist(dch, np.arange(-0.5, 0.5, 0.02), alpha=0.5);
        print len(dch[np.abs(dch) < 0.1]) / float(pixels)
    
    plt.title(cTitles[ci])
    plt.legend(['H=5', 'H=10', 'H=20'])
    plt.xlabel('concentration')
    plt.ylabel('pixels')
    plt.savefig('hist_%s_a.png' % cTitles[ci])
