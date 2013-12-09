#!/usr/bin/python
from boreali import lm, Boreali

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

from matplotlib import cm
cmap = cm.get_cmap('jet')


# RETRIEVE C FROM INSITU SPECTRA
# load insitu data
from scipy import interpolate

idata = np.loadtxt('rrsw0_spectra.txt')
iwavelen = idata[:, 0]
ied0 = idata[:, 1::6]
ied1 = idata[:, 2::6]
ilu0 = idata[:, 3::6]
ilu1 = idata[:, 4::6]


#wavelen = range(412, 709, 20)
#wavelen.append(709)
#wavelen = [412, 443, 490, 510, 555, 670]
wavelen = [413, 443, 490, 510, 560, 620, 665, 709]


ed0 = interpolate.interp1d(iwavelen, ied0.T)(wavelen)
ed1 = interpolate.interp1d(iwavelen, ied1.T)(wavelen)
lu0 = interpolate.interp1d(iwavelen, ilu0.T)(wavelen)
lu1 = interpolate.interp1d(iwavelen, ilu1.T)(wavelen)


rr0 = lu0 / ed0
rr1 = lu1 / ed1

#plt.plot(wavelen, rr0[0], 'o-')
#plt.plot(wavelen, rr1, 'o-')
#plt.show()
#raise

b = Boreali('michigan', wavelen)
# get matrix with HO-model
model = b.get_homodel()

hh = [
4.8,
9.8,
15,
20.4,
56.8,
12.6,
7.4,
5.1,
5.8,
]

tt = [
49.96000,
45.40000,
43.64000,
45.27000,
47.80000,
51.70000,
54.53000,
58.23000,
63.19000,
]

bb = [
0,
0,
0,
0,
0,
0,
0,
0,
4,
]

ccaqua = [
1.08,
0.78,
1.40,
1.11,
1.67,
2.56,
1.40,
1.06,
1.21,
]

cci = [
[1.76, 0.01, 0.01],
[1.45, 0.01, 0.01],
[1.48, 0.01, 0.01],
[1.38, 0.01, 0.01],
[1.39, 0.01, 0.01],
[1.28, 0.01, 0.01],
[1.34, 0.01, 0.01],
[1.58, 0.01, 0.01],
[1.77, 0.01, 0.01],
]
cci = np.array(cci)

aa = [b.get_albedo([bi])[0] for bi in bb]

aa = np.array(aa)


########### retrieve C shallow
parameters=[0.01, 5.0,
            0.01, 2.0,
            0.01, 0.5,
            10]

c0 = lm.get_c_shal(parameters, model, rr0, tt, hh, aa, 4*len(tt))[1]
c1 = lm.get_c_shal(parameters, model, rr1, tt, hh, aa, 4*len(tt))[1]


cc0 = np.array(c0).reshape(9, 4)[:, :-1]
cc1 = np.array(c1).reshape(9, 4)[:, :-1]

param1n = 2
param2n = 6

stations = [0,1,5,6,8]

#### OPTIMIZE TSM a and bb

a0 = model[param1n][:]
bb0 = model[param2n][:]

a0s = []
bb0s = []

for step in range(5):

    a0s.append(a0)
    bb0s.append(bb0)

    model[param1n][:] = a0[:]
    model[param2n][:] = bb0[:]

    mina = a0.min()
    maxa = a0.max()
    minbb = bb0.min()
    maxbb = bb0.max()
    
    avec  = np.r_[0.8*mina :1.2*maxa: 10j]
    bbvec = np.r_[0.8*minbb:1.2*maxbb:10j]

    abest = []
    bbbest = []
    for bn in range(len(wavelen)):
        bErr = np.ones((len(avec), len(bbvec)))
        for ci in stations:
            cErr = np.zeros((len(avec), len(bbvec)))
            for ai in range(len(avec)):
                for bbi in range(len(bbvec)):
                    model = b.get_homodel()
                    
                    model[param1n][bn] = avec[ai]
                    model[param2n][bn] = bbvec[bbi]
                    
                    params = parameters[:]
                    params[0] = cci[ci, 0]*0.95
                    params[1] = cci[ci, 0]*1.05
                    c1 = lm.get_c_shal(params, model, [rr1[ci]], [tt[ci]], [hh[ci]-1], [aa[ci]], 4)[1]
                    cErr[ai, bbi] = c1[3]
                    #cErr[ai, bbi] = (c1[0] - cci[ci, 0]) ** 2
    
            cErr = gaussian_filter(cErr, 1)
            plt.imshow(cErr, interpolation='nearest', extent =[bbvec[0], bbvec[-1], avec[0], avec[-1]])
            plt.colorbar()
            ax = plt.gca()
            ax.set_aspect((bbvec[-1]-bbvec[0])/(avec[-1]-avec[0]))
            plt.savefig('cerr_%02d_%02d_%02d.png' % (step, bn, ci))
            plt.close()
            
            bErr += cErr
    
        plt.imshow(bErr, interpolation='nearest', extent =[bbvec[0], bbvec[-1], avec[0], avec[-1]])
        plt.colorbar()
        ax = plt.gca()
        ax.set_aspect((bbvec[-1]-bbvec[0])/(avec[-1]-avec[0]))
        plt.savefig('berr_%02d_%02d.png' % (step, bn))
        plt.close()
        
        iBest = np.nonzero(bErr == bErr.min())
        abest.append(avec[iBest[0]])
        bbbest.append(bbvec[iBest[1]])
    
    a0 = np.array(abest).flatten()
    bb0 = np.array(bbbest).flatten()

a0s = np.array(a0s)
bb0s = np.array(bb0s)
plt.plot(wavelen, a0s.T);plt.legend(range(5));plt.show()
plt.plot(wavelen, bb0s.T);plt.legend(range(5));plt.show()

#%cpaste
model = b.get_homodel()
model[param1n][:] = abest[:]
model[param2n][:] = bbbest[:]

legval = []
#for ci in [0,5,6]:
for ci in stations:
    params = parameters[:]
    params[0] = cci[ci, 0]*0.95
    params[1] = cci[ci, 0]*1.05
    
    c1 = lm.get_c_shal(params, model, [rr1[ci]], [tt[ci]], [hh[ci]-1], [aa[ci]], 4)[1]
    print c1
    r1 = lm.get_rrsw_shal(model, c1[:-1], tt[ci], hh[ci]-1, aa[ci], len(wavelen))[1]
    plt.plot(wavelen, rr1[ci], '.-', c=cmap(ci/8.))
    plt.plot(wavelen, r1, '--', c=cmap(ci/8.))
    legval.append('%d input' % (ci + 1))
    legval.append('shallow: CHL=%5.2f TSM=%5.2f DOC=%5.2f' % (c1[0], c1[1], c1[2]))

x1,x2,y1,y2 = plt.axis()
plt.axis((x1, x2 * 1.1, y1, y2))
plt.xlabel('wavelength, nm')
plt.ylabel('Rrsw(-1), W m-2, nm')
plt.legend(legval)
plt.show()

