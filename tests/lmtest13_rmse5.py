#!/usr/bin/python
from lmtest13_setup import *

b = Boreali('michigan', wavelen)
# get matrix with HO-model
model = b.get_homodel()

rmser0 = []
rmsea0 = []
for v0 in v0vector:

    rmser1 = []
    rmsea1 = []
    for v1 in v1vector:

        rmser2 = []
        rmsea2 = []
        for v2 in v2vector:
        
            # get sand albedo
            albedo = b.get_albedo([v2])[0]
    
            #tough test of LM-retrieval
            #generate N random concentrations, depths, albedos, solarzenith angles
            pixels = 500
            cRange = np.array([2., 0.2, .1, 1., 1.])
            ccc = np.multiply(cRange, np.random.rand(pixels, 5))
            hhh = np.ones(pixels) * v1
            aaa = np.repeat(np.array([albedo]), pixels, axis=0)
            ttt = np.random.rand(pixels) * 10
            
            # calculate N Rrsw spectra
            rrr = []
            for i, c in enumerate(ccc):
                #r = lm.get_rrsw_deep(model, c, ttt[i], len(wavelen))[1]
                r = lm.get_rrsw_albe(model, c, ttt[i], hhh[i], wavelen, len(wavelen))[1]
                rrr.append(r)
            rrr = np.array(rrr)
            
            rrr2 = rrr
            hhh2 = hhh
            aaa2 = aaa
            
            # get noise in [-1, 1] interval
            noise = np.random.normal(size=(pixels))
            #noise /= np.abs(noise).max()

            # add random noise to Rrsw (%)
            rrr2 = rrr + np.multiply(noise, rrr.T).T * v0 / 2 / 100
        
            # add random noise to depth (meters)
            #hhh2 += noise * v0 / 2
            
            # add random noise to bottom albedo (%)
            #aaa2 += np.multiply(noise, aaa.T).T * v0 / 2 / 100
            
            # retrieve N concetrations
            t0 = time.time()
            #ccc2 = lm.get_c_deep(parameters, model, rrr, ttt, 4*pixels)[1]
            ccc2 = lm.get_c_albe(parameters, model, rrr2, ttt, hhh2, wavelen, 6*pixels)[1]
            t1 = time.time()
            print 'Time spent:', t1-t0
            ccc2 = np.array(ccc2).reshape(pixels, 6)
    
            #plt.plot(ccc[:, 0], ccc2[:, 0], '.')
            #plt.plot(ccc[:, 1], ccc2[:, 1], '.')
            #plt.plot(ccc[:, 2], ccc2[:, 2], '.');plt.show()
    
            # get relative and absolute RMSE and keep it
            rmser = 100 * np.sqrt(np.median(((ccc - ccc2[:, 0:5]) / ccc) ** 2, axis=0))
            rmsea = np.sqrt(np.median((ccc - ccc2[:, 0:5]) ** 2, axis=0))
            print rmsea
            
            rmser2.append(rmser)
            rmsea2.append(rmsea)
        
        rmser1.append(rmser2)
        rmsea1.append(rmsea2)

    rmser0.append(rmser1)
    rmsea0.append(rmsea1)

rmse = [rmser0, rmsea0]
rmsenp = np.array(rmse)
# dimensions = [absrel, noise, depth, bottype, cpa]
np.save(v0Name + 'rmse5.npy', rmsenp)

import lmtest13_plot
