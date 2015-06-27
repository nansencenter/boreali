#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
""" Boreali turorial
Import modules including Boreali and lm
Boreali is for processing satellite images
lm is Levenberg-Marquardt algorithm (c++). Used processing single pixels
within boreali.

Set all parameters and create a Boreali object

Use lm for calculation of Rrsw from given concentrations

Use lm for retrieval of vector of ceoncentrations from given Rrsw

Use Boreali for processing a test image and return 3D matrix with maps
of colour producing agents (CPA)

Create new nansat and save CPA

Create PNG images with maps of CPA
"""


import numpy as np
import matplotlib.pyplot as plt

from boreali import lm, Boreali
from nansat import *

# set forward model parameters: bands, albedo, depth, solar zenith angle
wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709] # MERIS bands
albedoType = 0 # sand bottom reflection
h = 2 # boreali will use bottom correction
theta = 0 # sun in zenith

# set LM parameters
cpaLimits=[0.01, 2,
            0.01, 1,
            0.01, 1, 10]

config = {
    'water': {'a':'a_smb', 'bb':'bb_smb'},
    'chl': {'a':'a_chl_MICHIGAN', 'bb':'bb_chl_MICHIGAN'},
    'tsm': {'a': 'a_tsm_MICHIGAN', 'bb': 'bb_tsm_MICHIGAN'},
    'doc': {'a': 'a_doc_MICHIGAN'}
}

# create Boreali object
b = Boreali('michigan', wavelen)
# get HO-model from boreali
model = b.get_homodel()

albedo = b.get_albedo([albedoType])[0]


# calculate Rrsw for 1 values of CHL in deep waters
chl = 0.05
# call lm.get_rrsw to calculate Rrs from input concentrations
rrsw = lm.get_rrsw_deep(model, [chl, 0.2, 0.01], theta, len(wavelen))[1]
# retrieve concentration vector from the last Rrsw
c = lm.get_c_deep(cpaLimits, model, [rrsw], [theta], 4)[1]
print 'chl=%5.2f, tsm=%5.2f, doc=%5.2f, rmse=%5.2f' % tuple(c)


# calculate and plot Rrsw for 6 values of CHL in shallow waters
chls = [0.01, 0.05, 0.1, 0.5, 1, 5]
legendVals = []
for chl in chls:
    # call lm.get_rrsw to calculate Rrs from input concentrations
    rrs = lm.get_rrsw_shal(model, [chl, 0.2, 0.01], theta, h, albedo, len(wavelen))[1]
    # plot R vs. wavelength
    plt.plot(wavelen, rrs, '.-')
    # add chl concentration to legends
    legendVals.append('%5.2f mg m-3' % chl)

    # retrieve concentration vector from the last Rrsw
    c = lm.get_c_shal(cpaLimits, model, [rrs], [theta], [h], [albedo], 4)[1]
    print 'chl=%5.2f, tsm=%5.2f, doc=%5.2f, rmse=%5.2f' % tuple(c)

# add legend, lables, title and save the plot to PNG file
plt.legend(legendVals)
plt.xlabel('wavelength, nm')
plt.ylabel('Rrsw, sr-1')
plt.title('Rrsw spectra for various chl values')
plt.savefig('test_rrsw.png')
plt.close()
raise

# open test image (subimage of lake michigan)
n = Nansat('test.tif')

# make a subset from input image
d = Domain(4326, '-lle -88 42.5 -87.33 42.75 -ts 200 100')
n.reproject(d)

# use Boreali object to process the input image, keep the result in cpa
cpa, mask = b.process(n, cpaLimits, threads=4)

# create new Nansat object to store the results of retrieval
nCPA = Nansat(domain=n)
# add bands with color producing agents to the nansat
nCPA.add_band(array=cpa[0], parameters={'name': 'chl', 'long_name': 'Chlorophyl-a', 'units': 'mg m-3'})
nCPA.add_band(array=cpa[1], parameters={'name': 'tsm', 'long_name': 'Total suspended matter', 'units': 'g m-3'})
nCPA.add_band(array=cpa[2], parameters={'name': 'doc', 'long_name': 'Dissolved organic carbon', 'units': 'gC m-3'})
nCPA.add_band(array=cpa[3], parameters={'name': 'mse', 'long_name': 'Root Mean Square Error', 'units': 'sr-1'})
nCPA.add_band(array=mask, parameters={'name': 'mask', 'long_name': 'L2 Boreali mask', 'units': '1'})

# save as GeoTIF file
nCPA.export('test_cpa.nc')

# generate PNG files with CPA spatial distribution and with RGB
figParams = {'legend': True, 'LEGEND_HEIGHT': 0.5, 'NAME_LOCATION_Y': 0, 'mask_array': mask, 'mask_lut': {1: [255, 255, 255], 2:[128,128,128], 4:[200,200,255]}}
nCPA.write_figure('test_chl.png', 'chl', clim=[0, 1.], **figParams)
nCPA.write_figure('test_tsm.png', 'tsm', clim=[0, 1.], **figParams)
nCPA.write_figure('test_doc.png', 'doc', clim=[0, .2], **figParams)
nCPA.write_figure('test_mse.png', 'mse', clim=[1e-5, 1e-2], logarithm=True, **figParams)
n.write_figure('test_rgb.png', [9, 5, 1], clim=[[0, 0, 0], [0.006, 0.04, 0.024]], mask_array=mask, mask_lut={2:[128,128,128]})
