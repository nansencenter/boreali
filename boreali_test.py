#!/usr/bin/python
from boreali import Boreali
from domain import Domain
from modis_l2_image import ModisL2Image
import numpy as np
import matplotlib.pyplot as plt
b = Boreali()

srsString = "+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs";
extentString = "-lle -10 50 10 65 -ts 50 50"

d = Domain(srs=srsString, ext=extentString)
d.write_map('/data/domain.png')

img = ModisL2Image('/data/gdal_test/A2012166115500.L2_LAC.NorthNorwegianSeas.hdf')
img.reproject(d)
img.write_figure('/data/rrs_412.png', 'Rrs_412', clim='hist', legend=True)
print img

cImg = b.process(img, [412, 443, 488, 531, 555, 667], start=10)

plt.imshow(cImg[0, :, :], vmin=0, vmax=5, interpolation='nearest');plt.colorbar();plt.show()
plt.imshow(cImg[1, :, :], vmin=0, vmax=3, interpolation='nearest');plt.colorbar();plt.show()
plt.imshow(cImg[2, :, :], vmin=0, vmax=3, interpolation='nearest');plt.colorbar();plt.show()
plt.imshow(np.log10(cImg[3, :, :]), vmin=-8, vmax=-5, interpolation='nearest');plt.colorbar();plt.show()
