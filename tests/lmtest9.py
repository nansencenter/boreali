#!/usr/bin/python
from nansat import *
from obpg_l2_image import OBPGL2Image

from boreali import Boreali
from numpy import ma

from lmtest0 import *

#"""
# retrieve CPA from in situ Rrsw
rrr = np.loadtxt('insitu_rrsw.txt')
aaa = np.loadtxt('insitu_albedo.txt')
hhh = [
4.8,
8.4,
13.5,
2.1,
2.4,
0.6,
4.8,
6.6,
6,
6.3,
4.8,
]
"""
hhh = [
-1,
-1,
-1,
-1,
-1,
-1,
-1,
-1,
-1,
-1,
-1
]
"""
ttt = [
0.7501425125,
0.692895713,
0.6571164634,
0.5904448859,
0.5806710421,
0.5750859885,
0.5775294495,
0.5921902152,
0.5921902152,
0.6740461571,
0.7283258969,
]

wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709]
bands = len(wavelen)
parameters=[len(wavelen), 10, 1e-6, 
            0.001, 12.,
            0.001, 1.,
            0.001, 1.]

b = Boreali('michigan', 'michigan')
model = []
for i in range(0, 8):
    model.append(b.homo[i](wavelen))

pixels = len(hhh)


ccc = lm.get_c(parameters, model, rrr, aaa, hhh, ttt, 4*pixels)[1]

cccr = ccc.reshape(pixels, 4)
print cccr
for iii in range(pixels):
    r = lm.get_rrsw(parameters, model, cccr[iii, 0:3], aaa[iii], hhh[iii], float(ttt[iii]), bands)[1]
    print np.array_str(r, 200)


#c = np.array([12.43792, 0.74, 1.34])
#r = lm.get_rrsw(parameters, model, c, aaa[iii], hhh[iii], float(ttt[iii]), bands)[1]
#print np.array_str(r, 200)
"""

lat = [
42.62658,
42.62559,
42.62386,
42.62663,
42.62676,
42.6272,
42.62671,
42.62737,
42.64528,
42.6634166667,
42.67293,
]

lon = [
-86.23351,
-86.23817,
-86.2464,
-86.22969,
-86.2299,
-86.22839,
-86.23252,
-86.23708,
-86.23238,
-86.2222666667,
-86.21919,
]
nBathy = Nansat('/files/michi/glbathy.tif')
n = OBPGL2Image('/files/michi/A2012242185000.L2_LAC_OC')
#d = Domain(4326, '-lle -86.5 42.5 -86 42.8 -ts 700 800')
d = Domain('+proj=stere +lat_0=42 +lon_0=-86 +ellps=WGS84 +datum=WGS84 +units=m +no_defs', '-lle -86.5 42.5 -86 42.8 -tr 100 100') #86 x 30 km (*10)

n.reproject(d)
wm = n.watermask()[1]
#n.write_figure('rrsw_555.png', 16, clim=[0, 0.01], mask_array=wm, mask_lut={0:[128,128,128]})

nBathy.reproject(d)
#nBathy.write_figure('bathy.png', 1, clim=[-300, 0], mask_array=wm, mask_lut={0:[128,128,128]})

rrsw = n[16] #rrsw_555
#rrsw = n[4] #rrsw_412
bathy = -0.3 * nBathy[1]
bathy[bathy > 9000] = 0

rrswM = ma.masked_where((wm==0) + (rrsw == 0), rrsw)

nMap = Nansatmap(n, figsize=(10,10), lat_0=42, lon_0=-86)
#nMap.pcolormesh(rrswM, vmin=-0.01, vmax=0.01) # 412
nMap.pcolormesh(rrswM, vmin=0, vmax=0.01) # 555

nMap.contour(bathy, fontsize=8, colors='k', levels=[10, 20, 30, 50])
nMap.add_colorbar(fontsize=10)
nMap.drawgrid()
x,y = nMap(lon, lat)
nMap.plot(x, y, 'ow')
plt.suptitle('MODIS Rrsw at 555 nm ')
nMap.save('rrsw_555_bathy.png', landmask=False)
#"""
