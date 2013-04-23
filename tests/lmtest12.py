#!/usr/bin/python
import os
import glob

import numpy as np
import matplotlib.pyplot as plt

from nansat import *
from obpg_l2_image import OBPGL2Image

from objectivezoning import ObjectiveZoning

""" Generate maps of zones for planning sampling campaign """

oz = ObjectiveZoning()

iDir = '/files/michi/raw/'
oDir = '/files/michi/'

ifile = iDir + 'subset_0_of_MER_FRS_1PNUPA20090713_160805_000005122080_00398_38528_7699.N1_C2IOP.nc'

wavelen = [413., 443., 490., 510., 560., 620., 665., 681., 709.]
#north
d = Domain(4326, '-te -86.3 44.8 -85.7 45.3 -ts 800 800')

#CPA
dCPA = Domain(4326, '-te -87 44.5 -85 45.5 -ts 640 400')
chl = np.load('/files/michi/subset_0_of_MER_FRS_1PNUPA20090713_160805_000005122080_00398_38528_7699.N1_C2IOP.nc0_sand.npy')
tsm = np.load('/files/michi/subset_0_of_MER_FRS_1PNUPA20090713_160805_000005122080_00398_38528_7699.N1_C2IOP.nc1_sand.npy')
doc = np.load('/files/michi/subset_0_of_MER_FRS_1PNUPA20090713_160805_000005122080_00398_38528_7699.N1_C2IOP.nc2_sand.npy')
cldmap = plt.imread('chladophoramap.png')
nCPA = Nansat(domain=dCPA)
nCPA.add_band(array=chl, parameters={'name': 'chl'})
nCPA.add_band(array=tsm, parameters={'name': 'tsm'})
nCPA.add_band(array=doc, parameters={'name': 'doc'})
nCPA.add_band(array=cldmap, parameters={'name': 'cldmap'})
nCPA.reproject(d)


nBathy = Nansat('/files/michi/glbathy.tif')                            
nBathy.reproject(d)
bathy = -0.3 * nBathy[1]
f = Figure(nparray=bathy, cmin=[0], cmax=[30])
f.process()
f.save('manitou.png')

wm = nBathy.watermask()[1]



n = OBPGL2Image(ifile)
    
# for the c2r image
n.reproject(d)

badMask = np.zeros(bathy.shape)
badMask[bathy > 40] = 1

bands = [bathy]
for w in wavelen:
    b = n['Rrsw_%03d' % w]
    badMask[b <= 0] = 1
    bands.append(b)

for cpaname in ['chl', 'tsm', 'doc', 'cldmap']:
    cpa = nCPA[cpaname]
    bands.append(cpa)
    badMask[cpa <= 0] = 1

bands = np.array(bands)
for band in bands:
    band[badMask == 1] = np.nan


bandsPCA = oz.pca(bands, 5)
zones = 6
zoneMap = oz.kmeans(bandsPCA, zones, 'test6f', norm=False)
zm1, zstd1, zn1 = oz.timeseries(bands, zoneMap)
plt.imshow(zoneMap);plt.colorbar(boundaries=range(1,zones+2), values=range(1,zones+1));plt.savefig('testzones_cbar.png')
plt.close()

plt.plot(zm1[0], zm1[10], 'o', zm1[0], zm1[11], '.',);plt.show()


