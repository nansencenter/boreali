#!/usr/bin/python
from nansat import *

from boreali import Boreali
from obpg_l2_image import OBPGL2Image


srsString = "+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs";
extentString = "-lle 120 27 127 37 -ts 500 1000"

d = Domain(srs=srsString, ext=extentString)
d.write_map('ecs_domain.png')

iFile = 'M2011217014149.L2_RR_OC.x.hdf'
iPath = '/files/yellowsea/'
iFilename = iPath + iFile

img = OBPGL2Image(iFilename)
img.reproject(d)

mask = img.watermask()[1]
img.write_figure(iFile + '_rgb.png', ['Rrsw_709', 'Rrsw_560', 'Rrsw_413'],
                 clim=[-0.005, 0.03],
                 mask_array=mask, mask_lut={2:[128,128,128]})
print img

# set forward model parameters: bands, albedo, depth, solar zenith angle
wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709] # MERIS bands
# create Boreali object
b = Boreali('michigan', wavelen)

# set LM parameters
parameters=[len(wavelen), 10, 0, 
            0.01, 15,
            0.01, 100,
            0.01, 5,]
# use Boreali object to process the input image, keep the result in cpa
cpa, mask = b.process(img, parameters, threads=2)

plt.imshow(cpa[0, :, :], vmin=0, vmax=5, interpolation='nearest');plt.colorbar();plt.show()
plt.imshow(cpa[1, :, :], vmin=0, vmax=3, interpolation='nearest');plt.colorbar();plt.show()
plt.imshow(cpa[2, :, :], vmin=0, vmax=3, interpolation='nearest');plt.colorbar();plt.show()
plt.imshow(np.log10(cpa[3, :, :]), vmin=-8, vmax=-1, interpolation='nearest');plt.colorbar();plt.show()
