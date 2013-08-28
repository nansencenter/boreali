#!/usr/bin/python
from lmtest0 import *

# select silicon sand albedo
albedo = albedos[2](wavelen)

# plot typical simulated Rrsw for shallow waters
for h0 in [1, 5, 10, 20]:
    #generate N random concentrations, depths, albedos, solarzenith angles
    pixels = 5000
    cRange = np.array([1., .5, .2])
    ccc = np.multiply(cRange, np.random.rand(pixels, 3))
    hhh = h0 * np.ones(pixels)
    aaa = np.repeat(np.array([albedo]), pixels, axis=0)
    ttt = 37 * np.ones(pixels)
    
    # calculate N Rrsw spectra
    rrr = []
    for i, c in enumerate(ccc):
        r = lm.get_rrsw(parameters, model, c, aaa[i], hhh[i], ttt[i], 6)[1]
        rrr.append(r)

    rrrn = np.array(rrr)
    rmean = rrrn.mean(axis=0)
    rstd = rrrn.std(axis=0)
    rmax = rrrn.max(axis=0)
    rmin = rrrn.min(axis=0)

    plt.close()    
    plt.plot(wavelen, rmean, 'o-r')
    plt.fill_between(wavelen, rmean+rstd, rmean-rstd, facecolor='red', alpha=0.3)
    plt.fill_between(wavelen, rmax, rmin, facecolor='red', alpha=0.3)
    a = plt.gca()
    a.set_ylim([0, 0.05])
    
    plt.xlabel('wavelength, nm')
    plt.ylabel('Rrsw, sr-1')
    plt.title('Mean Rrsw for H = %d m, C in range 1, .5, .2' % h0)
    
    plt.savefig('rrsw_mean_h%02d_algae8.png' % h0)
    plt.close()
