#!/usr/bin/python
import glob
import lm

import numpy as np
import matplotlib.pyplot as plt

from nansat import *
from obpg_l2_image import OBPGL2Image

from boreali import Boreali

iDir = '/files/michi/raw/'
oDir = '/files/michi/'

parameters=[6, 3, 1e-6, 
            0.01, 2,
            0.01, 1,
            0.01, 1,]
b = Boreali('michigan', 'michigan')
bands = 6
wavelen = [412, 443, 490, 510, 555, 670]
# get matrix with HO-model
model = []
for i in range(0, 8):
    model.append(b.homo[i](wavelen))

print model

print b.albedos[0](wavelen)
print b.albedos[1](wavelen)
print b.albedos[2](wavelen)

d = Domain(4326, '-te -87 44.5 -85 45.5 -ts 640 400')

iFiles = glob.glob(iDir + 'subset_0_of_MER_FRS_1PNUPA20090713_160805_000005122080_00398_38528_7699.N1_C2IOP.nc')
nBathy = Nansat('/files/michi/glbathy.tif')

dCPA = Domain(4326, '-te -87 44.5 -85 45.5 -ts 640 400')
bot = plt.imread('bottomtype.png')
nCPA = Nansat(domain=dCPA)
nCPA.add_band(array=cldmap, parameters={'name': 'cldmap'})

for ifile in iFiles:
    n = OBPGL2Image(ifile)
    mask = n.l2c_mask([1, 4, 5, 6, 9, 10, 15, 20, 29])
    n.add_band(array=mask, parameters={'name': 'mask'})
    
    nBathy.reproject(n)
    depth = nBathy[1]
    n.add_band(array=depth, parameters={'name': 'depth'})
    print n

"""
#create albedo interpolators
albedoString = '''
lambda	400	425	450	475	500	525	550	575	600	625	650	675	700	725	750
sand	0.11	0.13	0.15	0.165	0.183	0.196	0.209	0.222	0.235	0.248	0.261	0.274	0.287	0.3	0.31
sargassum	0.02	0.02	0.02	0.02	0.02	0.03	0.04	0.06	0.07	0.05	0.05	0.02	0.15	0.31	0.35
silt	0.05	0.065	0.08	0.095	0.11	0.125	0.14	0.155	0.17	0.175	0.18	0.185	0.19	0.195	0.2
boodlea	0.05	0.05	0.06	0.07	0.11	0.25	0.28	0.25	0.22	0.19	0.15	0.08	0.25	0.43	0.45
limestone	0.065	0.075	0.09	0.1	0.115	0.13	0.14	0.15	0.145	0.14	0.14	0.125	0.135	0.1	0.07
enteromorpha	0.05	0.05	0.05	0.05	0.07	0.135	0.165	0.155	0.145	0.13	0.1	0.075	0.22	0.3	0.31
cladophora	0.05	0.05	0.05	0.05	0.07	0.135	0.165	0.155	0.145	0.13	0.1	0.075	0.14	0.175	0.18
chara	0.015	0.015	0.015	0.015	0.015	0.04	0.05	0.04	0.035	0.03	0.025	0.02	0.05	0.14	0.16
'''

albedoLines = albedoString.splitlines()
albedoNames = []
x = np.array(map(float, albedoLines[1].split('\t')[1:]))
albedos = []
for awi in range(2, len(albedoLines)):
    awiList = albedoLines[awi].split('\t')
    albedoNames.append(awiList[0])
    y = np.array(map(float, awiList[1:]))
    albedos.append(interpolate.interp1d(x, y))

#sun at typical zenith angle
theta = 25
"""
