#!/usr/bin/python
import glob

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

from boreali import Boreali
import lm

from nansat import *
from obpg_l2_image import OBPGL2Image

from lmtest0 import *

iDir = '/files/michi/'
#iFiles = glob.glob(iDir + 'M*.L2_FRS_OC')
iFiles = glob.glob(iDir + '*.nc')

b = Boreali('michigan', 'michigan')


nBathy = Nansat('/files/michi/glbathy.tif')
#southern part
#d = Domain(4326, '-te -88 41.5 -86 44 -ts 600 1000')
#chicago
d = Domain(4326, '-te -88 41.5 -87 42.3 -ts 1000 1000')
nBathy.reproject(d)
nBathy.write_figure(iDir + 'bathymetry.png', 1, clim=[-60, 0], legend=True, caption='Depth [ft]', fontSize=20)
bathy = -0.3 * nBathy[1]

wm = nBathy.watermask()[1]

cloudBits = [1, 4, 5, 6, 9, 10, 15, 20, 29]

wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709]
parameters=[len(wavelen), 3, 1e-6, 
            0.001, 2.,
            0.001, .5,
            0.001, .5,]

# choose albedo
aNo = 3
albedo = albedos[aNo](wavelen)

model = []
for i in range(0, 8):
    model.append(b.homo[i](wavelen))

michipins = '''
pin_1	379.50134	415.49777	-87.6205	41.9676	Pin 1	
pin_2	409.50012	397.50195	-87.5905	41.982	Pin 2	
pin_3	439.4989	381.4993	-87.5605	41.9948	Pin 3	
pin_6	538.4979	318.49957	-87.4615	42.0452	Pin 6	
pin_10	666.4963	227.49995	-87.3335	42.118	Pin 10	
'''

pinX, pinY = [],[]
for michipin in michipins.splitlines()[1:]:
    pinX.append(int(float(michipin.split('\t')[1])))
    pinY.append(int(float(michipin.split('\t')[2])))


for ifile in iFiles:
    print ifile
    n = OBPGL2Image(ifile)

    # for meris L2 from OBPG
    #n.reproject(d, tps=True, use_geolocationArray=False)
    #mask = n.l2c_mask(cloudBits, [2], 'Rrs_413')

    # for the c2r image
    n.reproject(d)
    mask = np.ones(wm.shape)*64
    mask[wm == 0] = 2
    gpi = np.ravel_multi_index([pinY, pinX], n.shape())
    pixels = len(mask.flat[gpi])

    # generate remote sensing reflectance
    rrr = []
    for w in wavelen:
        print w
        rrsw = n['Rrsw_%03d' % w]
        rrsw = rrsw.flat[gpi]
        rrr.append(rrsw)
    rrr = np.array(rrr).T

    # generate sun zenith
    ttt = n['sun_zenith']
    ttt = ttt.flat[gpi]
    
    # generate albedo
    aaa = np.repeat([albedo], pixels, axis=0)
    
    # generate depth
    maxDepth = 20
    hhh = np.copy(bathy)
    hhh = hhh.flat[gpi]
    hhh[hhh > maxDepth] = -100 # mark deep waters
"""
    ccc = lm.get_c(parameters, model, rrr, aaa, hhh, ttt, 4*pixels)[1]
    cccr = ccc.reshape(pixels, 4)

legendStr = []
for iii in range(pixels):
    legendStr.append('%02dm %6.3f %6.3f %6.3f' % (
                        np.round(bathy.flat[gpi[iii]]),
                        cccr[iii, 0],
                        cccr[iii, 1],
                        cccr[iii, 2]))
    legendStr.append('reconstructed')
                        
clrs = ['r', 'g', 'b', 'c', 'm', 'k']
for iii in range(pixels):
    r = lm.get_rrsw(parameters, model, cccr[iii, 0:3], aaa[iii], hhh[iii], float(ttt[iii]), 9)[1]
    plt.plot(wavelen, rrr[iii], 'o-'+clrs[iii], wavelen, r, '--'+clrs[iii])
plt.legend(legendStr);plt.show();
"""

# plot 5m spectrum with diff bottomtypes
# choose albedo
iii = 0
plt.close()
plt.plot(wavelen, rrr[iii], 'o-')
legStrings = ['measured', 'sand    ', 'algae5  ', 'silt    ', 'algae6  ', 'deep    ']

for aNo in range(4):
    albedo = albedos[aNo](wavelen)

    ccc = lm.get_c(parameters, model, [rrr[iii]], [albedo], [hhh[iii]], [ttt[iii]], 4)[1]
    print ccc
    r = lm.get_rrsw(parameters, model, ccc[0:3], albedo, hhh[iii], float(ttt[iii]), 9)[1]
    plt.plot(wavelen, r, 'o-')
    legStrings[aNo+1] += '%5.2f %5.2f %5.2f %5.2e' % (ccc[0], ccc[1], ccc[2], ccc[3])

ccc = lm.get_c(parameters, model, [rrr[iii]], [albedo], [hhh[iii]], [ttt[iii]], 4*pixels)[1]
r = lm.get_rrsw(parameters, model, ccc[0:3], albedo, -90, float(ttt[iii]), 9)[1]
legStrings[5] += '%5.2f %5.2f %5.2f %5.2e' % (ccc[0], ccc[1], ccc[2], ccc[3])
plt.plot(wavelen, r, 'o-')
    
plt.legend(legStrings)
plt.show();

