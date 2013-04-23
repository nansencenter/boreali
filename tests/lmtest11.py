#!/usr/bin/python
import glob

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from nansat import *
from boreali import *

""" Test Boreali """

iDir = '/files/michi/raw/'
oDir = '/files/michi/'

# set wavelengths to use
wavelen = [413.0, 443.0, 490.0, 510.0, 560.0, 620.0, 665.0, 681.0, 709.0]
opts=[9, 3, 1e-6, 
      0.01, 2,
      0.01, 1,
      0.01, 1,]
      
b = Boreali('michigan', 'michigan')

#dstDomain = Domain(4326, '-te -87 44.5 -85 45.5 -ts 640 400')
dstDomain = Domain(4326, '-te -87 44.5 -85 45.5 -ts 60 40')


# get projected bathymetry
iFiles = glob.glob(iDir + 'subset_0_of_MER_FRS_1PNUPA20090713_160805_000005122080_00398_38528_7699.N1_C2IOP.nc')
nBathy = Nansat('/files/michi/glbathy.tif')
nBathy.reproject(dstDomain)
depth = -0.3 * nBathy[1]

# get projected bottom type
dCPA = Domain(4326, '-te -87 44.5 -85 45.5 -ts 640 400')
cldmap = sp.misc.imread('bottomtype.png')
nCPA = Nansat(domain=dCPA)
nCPA.add_band(array=cldmap, parameters={'name': 'cldmap'})
nCPA.reproject(dstDomain)
bottom = nCPA[1]

for ifile in iFiles:
    # get projected input image
    n = Nansat(ifile)
    n.reproject(dstDomain)

    # get mask
    Rrsw_560 = n['Rrsw_560']
    mask = np.zeros(Rrsw_560.shape, 'uint8') + 1
    mask[Rrsw_560 > 0] = 64
    mask[nBathy[1] == -32768] = 2
    
    # get sun zenith
    theta = n['sun_zenith']
    
    # run processing and get CPA
    cpa = b.process(n, opts, wavelen=wavelen, mask=mask, depth=depth, bottom=bottom, theta=theta)
    cpa = b.process(n, opts, wavelen=wavelen, mask=mask, depth=depth, bottom=bottom)
    cpa = b.process(n, opts, wavelen=wavelen, mask=mask, depth=depth)
    cpa = b.process(n, opts, wavelen=wavelen, mask=mask)
    cpa = b.process(n, opts, wavelen=wavelen)
    cpa = b.process(n, opts)
