#!/usr/bin/python
from lmtest0 import *

parameters=[6, 3, 1e-6, 
            0.01, 1.,
            0.01, .5,
            0.01, .5,]


wavelen = [413., 443., 490., 510., 560., 620., 665., 681., 709.]

# get matrix with HO-model
model = []
for i in range(0, 8):
    model.append(b.homo[i](wavelen))


dcccH = []

for h0 in [5, 10, 20]:
    #tough test of LM-retrieval
    #generate N random concentrations, depths, albedos, solarzenith angles
    pixels = 1000
    cRange = np.array([1., .5, .5, 1., 1.])
    ccc = np.multiply(cRange, np.random.rand(pixels, 5))
    #ccc[:,3] = 0.5
    #ccc[:,4] = 0.5
    
    #hhh = 5 + np.random.rand(pixels) * 5  # H in range 5 - 10 m
    hhh = h0 * np.ones(pixels)
    ttt = np.random.rand(pixels) * 10

    # calculate N Rrsw spectra
    rrr = []
    for i, c in enumerate(ccc):
        r = lm.get_rrsw_al(parameters, model, c[0:3], wavelen, hhh[i], ttt[i], c[3], c[4], len(wavelen))[1]
        rrr.append(r)
    
    # add random noise to depth (+/- 1 m)
    hhh2 = hhh - 1 + np.random.rand(pixels) * 2
    dh = (hhh - hhh2) ** 2
    
    # retrieve N concetrations
    ccc2 = lm.get_c_al(parameters, model, rrr, wavelen, hhh, ttt, 6*pixels)[1]
    ccc2 = np.array(ccc2).reshape(pixels, 6)
    
    # get relative difference and keep it
    dccc = (ccc - ccc2[:, 0:5]) / ccc
    dcccH.append(dccc)

# plot comparison of concentrations
cTitles = ['CHL-CHL2', 'TSM-TSM2', 'DOC-DOC2', 'AL1-AL12', 'AL2-AL22']
for ci in [0, 1, 2, 3, 4]:
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
    plt.savefig('hist_%s_c.png' % cTitles[ci])
