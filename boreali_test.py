#!/usr/bin/python
from nansat import *

from boreali import Boreali
from obpg_l2_image import OBPGL2Image


srsString = "+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs";
extentString = "-lle 120 27 127 37 -ts 500 1000"

d = Domain(srs=srsString, ext=extentString)
d.write_map('ecs_domain.png')

iFile = 'M2011217014149.L2_RR_OC.x.hdf'
iPath = './'
iFilename = iPath + iFile

img = OBPGL2Image(iFilename)
l2c_mask = img.l2c_mask()
img.write_figure(iFile + '_raw.png', ['Rrsw_709', 'Rrsw_560', 'Rrsw_413'],
                 clim=[[0, 0, 0], [0.007, 0.05, 0.03]],
                 mask_array=l2c_mask, mask_lut={2:[128,128,128], 1:[255,255,255]})

img.reproject(d)

l2c_mask = img.l2c_mask()
img.write_figure(iFile + '_rgb.png', ['Rrsw_709', 'Rrsw_560', 'Rrsw_413'],
                 clim=[[0, 0, 0], [0.007, 0.05, 0.03]],
                 mask_array=l2c_mask, mask_lut={2:[128,128,128], 1:[255,255,255]})
print img

# set forward model parameters: bands, albedo, depth, solar zenith angle
wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709] # MERIS bands
# create Boreali object
b = Boreali('michigan', wavelen)

# set LM parameters
parameters=[0.01, 15,
            0.01, 100,
            0.01, 5,]
# use Boreali object to process the input image, keep the result in cpa
cpa, mask = b.process(img, parameters, threads=2)

mask[l2c_mask == 2] = 2
mask[l2c_mask == 1] = 1

nCPA = Nansat(domain=d)
# add bands with color producing agents to the nansat
nCPA.add_band(array=cpa[0], parameters={'name': 'chl', 'long_name': 'Chlorophyl-a', 'units': 'mg m-3'})
nCPA.add_band(array=cpa[1], parameters={'name': 'tsm', 'long_name': 'Total suspended matter', 'units': 'g m-3'})
nCPA.add_band(array=cpa[2], parameters={'name': 'doc', 'long_name': 'Dissolved organic carbon', 'units': 'gC m-3'})
nCPA.add_band(array=cpa[3], parameters={'name': 'mse', 'long_name': 'Root Mean Square Error', 'units': 'sr-1'})
nCPA.add_band(array=mask, parameters={'name': 'mask', 'long_name': 'L2 Boreali mask', 'units': '1'})

# save as GeoTIF file
nCPA.export('test_cpa.nc')

figParams = {'legend': True, 'LEGEND_HEIGHT': 0.1, 'NAME_LOCATION_Y': 0, 'mask_array': mask, 'mask_lut': {1: [255, 255, 255], 2:[128,128,128], 4:[200,200,255]}}
nCPA.write_figure(iFile +  '_chl.png', 'chl', clim=[0, 5], **figParams)
nCPA.write_figure(iFile +  '_tsm.png', 'tsm', clim=[0, 15], **figParams)
nCPA.write_figure(iFile +  '_doc.png', 'doc', clim=[0, 1], **figParams)
nCPA.write_figure(iFile +  '_mse.png', 'mse', clim=[1e-5, 1e-1], logarithm=True, **figParams)
