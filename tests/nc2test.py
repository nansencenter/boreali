from nansat import *
bands = [
'Rrsw_413',
'Rrsw_443',
'Rrsw_490',
'Rrsw_510',
'Rrsw_560',
'Rrsw_620',
'Rrsw_665',
'Rrsw_681',
'Rrsw_709',
'sun_zenith',
'chl_conc',
'tsm',
'a_ys_443',
]

n = Nansat('/files/michi/raw/subset_0_of_MER_FRS_1PNUPA20090716_161348_000005122080_00441_38571_7871.N1_C2IOP.nc')
#southern part
d = Domain(4326, '-te -88 42.5 -87 42.75 -ts 300 100')

n.reproject(d)
#Rrsw_560 = n['Rrsw_560']
#plt.imshow(Rrsw_560);plt.colorbar();plt.show()

#'''
n2 = Nansat(domain=d)
for b in bands:
    bmeta = n.get_metadata(bandID = b)
    bdata = n[b]
    print bmeta
    n2.add_band(array=bdata, parameters=bmeta)

Rrsw_413 = n['Rrsw_413']
Rrsw_413[Rrsw_413 > 0] = 64
n2.add_band(array=Rrsw_413, parameters={'name': 'mask'})
n2.export('test.tif', driver='GTiff')

#'''
