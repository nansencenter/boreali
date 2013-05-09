#!/usr/bin/python
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
albedo = [0, 0, 0, 0, 0, 0, 0, 0, 0] # no bottom reflection
h = -10 # boreali will not use bottom correction if h < 0
theta = 0 # sun in zenith

# set LM parameters
parameters=[len(wavelen), 3, 1e-6, 
            0.01, 2,
            0.01, 1,
            0.01, 1,]

# create Boreali object
b = Boreali('michigan', wavelen)
# get HO-model from boreali
model = b.get_homodel()

# calculate and plot Rrsw for 6 values of CHL 
chls = [0.01, 0.05, 0.1, 0.5, 1, 5]
legendVals = []
for chl in chls:
    # call lm.get_rrsw to calculate Rrsw from input concentrations
    r = lm.get_rrsw(parameters, model, [chl, 0.01, 0.01], albedo, h, theta, len(wavelen))[1]
    # plot R vs. wavelength
    plt.plot(wavelen, r, '.-')
    # add chl concentration to legends
    legendVals.append('%5.2f mg m-3' % chl)

# add legend, lables, title and save the plot to PNG file
plt.legend(legendVals)
plt.xlabel('wavelength, nm')
plt.ylabel('Rrsw, sr-1')
plt.title('Rrsw spectra for various chl values')
plt.savefig('test_rrsw.png')
plt.close()

# retrieve concentration vector from the last Rrsw
c = lm.get_c(parameters, model, [r], [albedo], [h], [theta], 4)[1]
print 'chl=%5.2f, tsm=%5.2f, doc=%5.2f, rmse=%5.2f' % tuple(c)

# open test image (subimage of lake michigan)
n = Nansat('test.tif')
# use Boreali object to process the input image, keep the result in cpa
cpa = b.process(n, parameters)

# create new Nansat object to store the results of retrieval
nCPA = Nansat(domain=n)
# add bands with color producing agents to the nansat
nCPA.add_band(array=cpa[0], parameters={'name': 'chl', 'long_name': 'Chlorophyl-a', 'units': 'mg m-3'})
nCPA.add_band(array=cpa[1], parameters={'name': 'tsm', 'long_name': 'Total suspended matter', 'units': 'g m-3'})
nCPA.add_band(array=cpa[2], parameters={'name': 'doc', 'long_name': 'Dissolved organic carbon', 'units': 'gC m-3'})
nCPA.add_band(array=cpa[3], parameters={'name': 'mse', 'long_name': 'Root Mean Square Error', 'units': 'sr-1'})
# save as GeoTIF file
nCPA.export('test_cpa.tif', driver='GTiff')

# get watermask for plotting
wm = n.watermask()[1]

# generate PNG files with CPA spatial distribution and with RGB
nCPA.write_figure('test_chl.png', 'chl', clim=[0, 0.5], legend=True, LEGEND_HEIGHT=0.5, NAME_LOCATION_Y=0, mask_array=wm, mask_lut={2:[128,128,128]})
nCPA.write_figure('test_tsm.png', 'tsm', clim=[0, 1], legend=True, LEGEND_HEIGHT=0.5, NAME_LOCATION_Y=0, mask_array=wm, mask_lut={2:[128,128,128]})
nCPA.write_figure('test_doc.png', 'doc', clim=[0, .5], legend=True, LEGEND_HEIGHT=0.5, NAME_LOCATION_Y=0, mask_array=wm, mask_lut={2:[128,128,128]})
nCPA.write_figure('test_mse.png', 'mse', clim=[1e-5, 1e-2], logarithm=True, LEGEND_HEIGHT=0.5, NAME_LOCATION_Y=0, legend=True, mask_array=wm, mask_lut={2:[128,128,128]})
n.write_figure('test_rgb.png', [9, 5, 1], clim=[[0, 0, 0], [0.006, 0.04, 0.024]], mask_array=wm, mask_lut={2:[128,128,128]})
