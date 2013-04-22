#!/usr/bin/python
from lmtest0 import *

# plot typical real Rrsw for shallow waters

from nansat import *

from obpg_l2_image import OBPGL2Image


import glob
iDir = '/files/michi/'
iFiles = glob.glob(iDir + 'M*.L2_FRS_OC')
iFiles = glob.glob(iDir + '*.nc')

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


h0s = [0, 5, 10]
h1s = [5, 10, 20]
wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709]

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
    #n.write_figure(ifile + '_Rrsw_413_masked.png', 'Rrsw_413', clim=[-0.001, 0.02], legend=True, fontSize=20, mask_array=mask, mask_lut={0:[100]*3, 1:[255]*3, 2:[128]*3})
    n.write_figure(ifile + '_Rrsw_560_masked.png', 'Rrsw_560', clim=[-0.001, 0.02], legend=True, fontSize=20, mask_array=mask, mask_lut={0:[100]*3, 1:[255]*3, 2:[128]*3})
    #n.write_figure(ifile + '_Rrsw_709.png', 'Rrsw_709', clim='hist', legend=True, fontSize=20)
    #n.write_figure(ifile + '_RGB.png', ['Rrsw_709', 'Rrsw_560', 'Rrsw_413'], clim='hist')
    raise

    for h0, h1 in zip(h0s, h1s):
        print h0
        rr = []
        for w in wavelen:
            print w
            rrsw = n['Rrsw_%03d' % w]
            subrrsw = rrsw[(mask == 64) * (bathy > h0) * (bathy <= h1)]
            rr.append([subrrsw.mean(), subrrsw.std(), subrrsw.min(), subrrsw.max()])
    
        rrn = np.array(rr).T
    
        plt.close()
        plt.plot(wavelen, rrn[0], 'o-b')
        plt.fill_between(wavelen, rrn[0]+rrn[1], rrn[0]-rrn[1], facecolor='b', alpha=0.3)
        plt.fill_between(wavelen, rrn[2], rrn[3], facecolor='b', alpha=0.3)
        a = plt.gca()
        a.set_ylim([0, 0.05])
        
        plt.xlabel('wavelength, nm')
        plt.ylabel('Rrsw, sr-1')
        plt.title('Mean Rrsw for %d < H < %d m' % (h0, h1))
        
        plt.savefig(ifile + 'rrsw_mean_h%02db.png' % h0)
        plt.close()
    
    
