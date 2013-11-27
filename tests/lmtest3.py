#!/usr/bin/python
from lmtest0 import *
import time

# select silicon sand albedo
albedo = albedos[0](wavelen)

dcccH = []

for h0 in [2, 5, 10, 20]:
    #simple test of LM-retrieval
    #generate N random concentrations, depths, albedos, solarzenith angles
    pixels = 5
    cRange = np.array([1., .5, .5])
    ccc = np.multiply(cRange, np.random.rand(pixels, 3))
    hhh = h0 * np.ones(pixels)
    aaa = np.repeat(np.array([albedo]), pixels, axis=0)
    ttt = np.random.rand(pixels) * 10
    
    # calculate N Rrsw spectra
    rrr = []
    for i, c in enumerate(ccc):
        r = lm.get_rrsw(model, c, aaa[i], hhh[i], ttt[i], len(wavelen))[1]
        # append r to list
        rrr.append(r)
    
    # retrieve N concetrations
    t0 = time.time()
    ccc2 = lm.get_c(parameters, model, rrr, aaa, hhh, ttt, 4*pixels)[1]
    t1 = time.time()
    print 'Time spent:', t1-t0
    ccc2 = np.array(ccc2).reshape(pixels, 4)
    
    # print results
    for ic, oc in zip(ccc, ccc2):
        print ic, oc
    
