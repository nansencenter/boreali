#!/usr/bin/python
import os
import glob

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

from boreali import Boreali
import lm

from nansat import *
from obpg_l2_image import OBPGL2Image

from lmtest0 import *

iDir = '/files/michi/raw/'
oDir = '/files/michi/'
#iFiles = glob.glob(iDir + 'M*.L2_FRS_OC')
#iFiles = glob.glob(iDir + '*.nc')
iFiles = glob.glob(iDir + 'subset_0_of_MER_FRS_1PNUPA20090713_160805_000005122080_00398_38528_7699.N1_C2IOP.nc')
#iFiles = glob.glob(iDir + 'subset_1_of_MER_FRS_1PNUPA20090720_154806_000005122080_00498_38628_8130.N1_C2IOP.nc')

b = Boreali('michigan', 'michigan')


nBathy = Nansat('/files/michi/glbathy.tif')
#southern part
#d = Domain(4326, '-te -88 41.5 -86 44 -ts 600 1000')
#chicago
#d = Domain(4326, '-te -88 41.5 -87 42.3 -ts 1000 1000')
#north
d = Domain(4326, '-te -87 44.5 -85 45.5 -ts 640 400')

nBathy.reproject(d)
nBathy.write_figure(oDir + 'north_bathymetry.png', 1, clim=[-60, 0], legend=True, caption='Depth [ft]', fontSize=20)
bathy = -0.3 * nBathy[1]

wm = nBathy.watermask()[1]

cloudBits = [1, 4, 5, 6, 9, 10, 15, 20, 29]

wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709]
parameters=[len(wavelen), 3, 1e-6, 
            0.001, 2,
            0.001, .5,
            0.001, .2]


for aNo in [4, 5]:
    # choose albedo    
    albedo = albedos[aNo](wavelen)
    
    #choose CPA limits for visualization
    cMin = [0, 0, 0, 1e-5]
    cMax = [2, 0.3, 0.1, 2e-3]
    
    maxDepth = 30
    
    model = []
    for i in range(0, 8):
        model.append(b.homo[i](wavelen))
    
    for ifile in iFiles:
        ifilename = os.path.split(ifile)[1]
        print ifile
        n = OBPGL2Image(ifile)
    
        # for meris L2 from OBPG
        #n.reproject(d, tps=True, use_geolocationArray=False)
        #mask = n.l2c_mask(cloudBits, [2], 'Rrs_413')
    
        # for the c2r image
        n.reproject(d)
        mask = np.ones(wm.shape)*64
        mask[wm == 0] = 2
        mask[bathy > maxDepth] = 1
        gpi = mask == 64
        pixels = len(mask[gpi])
    
        # generate remote sensing reflectance
        rrr = []
        for w in wavelen:
            print w
            rrsw = n['Rrsw_%03d' % w]
            rrsw = rrsw[gpi]
            rrr.append(rrsw)
        rrr = np.array(rrr).T
    
        # generate sun zenith
        ttt = n['sun_zenith']
        ttt = ttt[gpi]
        
        # generate albedo
        aaa = np.repeat([albedo], pixels, axis=0)
        
        # generate depth
        hhh = np.copy(bathy)
        hhh = hhh[gpi]
        hhh[hhh > maxDepth] = -100 # mark deep waters
        
        ccc = lm.get_c(parameters, model, rrr, aaa, hhh, ttt, 4*pixels)[1]
        cccr = ccc.reshape(pixels, 4)
        for ci in range(0, 4):
            c = np.zeros(bathy.shape)
            c[gpi] = cccr[:, ci]
            np.save(oDir + ifilename + '%d_%s' % (ci, albedoNames[aNo]), c)
            f = Figure(nparray=c, cmin=[cMin[ci]], cmax=[cMax[ci]], mask_array=wm, mask_lut={0:[128,128,128]}, legend=True)
            f.process()
            f.save(oDir + ifilename + '%d_%s.png' % (ci, albedoNames[aNo]))
    
    
