#!/usr/bin/python
from boreali import lm, Boreali

import numpy as np
import matplotlib.pyplot as plt

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
wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709]

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

b = Boreali('michi2', wavelen)
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
            100]

c0 = lm.get_c_shal(parameters, model, rr0, tt, hh, aa, 4*len(tt))[1]
c1 = lm.get_c_shal(parameters, model, rr1, tt, hh, aa, 4*len(tt))[1]


cc0 = np.array(c0).reshape(9, 4)[:, :-1]
cc1 = np.array(c1).reshape(9, 4)[:, :-1]

"""
######## plot Ed(0) vs Ed_est(0)
legval = []
for ci in [0,1,2,3,4,5,6,7,8]:
    kd0 = b.get_kd(cc0[ci])
    kd1 = b.get_kd(cc1[ci])
    ed1c = b.get_attenuated_ed(ed0[ci], kd0, 1.)
    ed0c = ed1[ci] / np.exp(- kd1 * 1.0)
    plt.plot(wavelen, ed0[ci], '.-', c=cmap(ci/8.))
    plt.plot(wavelen, ed1[ci], '-', c=cmap(ci/8.))
    plt.plot(wavelen, ed1c, '--', c=cmap(ci/8.))
    plt.plot(wavelen, ed0c, '-.', c=cmap(ci/8.))
    
    legval.append('%d Ed(0)' % (ci + 1))
    legval.append('%d Ed(-1)' % (ci + 1))
    legval.append('%d Ed_est(-1)' % (ci + 1))
    legval.append('%d Ed_est(0)' % (ci + 1))
    plt.title('Station %d' % (ci + 1))
    plt.xlabel('wavelength, nm')
    plt.ylabel('Ed, W m-2, nm')
    plt.legend(legval)
    plt.savefig('ed_vs_ed_%02d.png' % (ci + 1))
    plt.close()
"""




"""
######### plot Rrsw retr vs Rrsw insitu
legval = []
for ci in [0,1,2,3,4,5,6,7,8]:
    r1 = lm.get_rrsw_shal(model, cc1[ci], tt[ci], hh[ci], aa[ci], len(wavelen))[1]
    plt.plot(wavelen, rr1[ci], '.-', c=cmap(ci/8.))
    plt.plot(wavelen, r1, '--', c=cmap(ci/8.))
    legval.append('%d input' % (ci + 1))
    legval.append('from: CHL=%5.2f TSM=%5.2f DOC=%5.2f' % (cc0[ci, 0], cc0[ci, 1], cc0[ci, 2]))

x1,x2,y1,y2 = plt.axis()
plt.axis((x1, x2 * 1.1, y1, y2))
plt.xlabel('wavelength, nm')
plt.ylabel('Rrsw(-1), W m-2, nm')
plt.legend(legval)
plt.show()
"""




"""
######### plot C vs C
left=np.r_[0:9]*5
plt.bar(left, cc[:, 0], 1, color=[1,0,0])

left=np.r_[0:9]*5 + 1
plt.bar(left, ccaqua, 1, color=[1,1,0])

left=np.r_[0:9]*5 + 2
plt.bar(left, cc0[:, 0], 1, color=[0,1,0])

left=np.r_[0:9]*5 + 3
plt.bar(left, cc1[:, 0], 1, color=[0,0,1])

plt.legend(['lab', 'aquaflor', 'from -0', 'from -1'])
plt.show()
"""




#### plot Rrsw_insitu vs Rrsw_rec(C_insitu)
legval = []

#for ci in [0,5,6]:
for ci in range(9):
    params = parameters[:]
    #params[0] = cci[ci, 0]*0.95
    #params[1] = cci[ci, 0]*1.05
    c1 = lm.get_c_shal(params, model, [rr1[ci]], [tt[ci]], [hh[ci]-1], [aa[ci]], 4)[1]
    #c1 = lm.get_c_deep(params, model, [rr1[ci]], [tt[ci]], 4)[1]
    r1 = lm.get_rrsw_shal(model, c1[:-1], tt[ci], hh[ci]-1, aa[ci], len(wavelen))[1]
    rd = lm.get_rrsw_deep(model, c1[:-1], tt[ci], len(wavelen))[1]
    plt.plot(wavelen, rr1[ci], '.-', c=cmap(ci/8.))
    plt.plot(wavelen, r1, '--', c=cmap(ci/8.))
    #plt.plot(wavelen, rd, '-.', c=cmap(ci/8.))
    legval.append('%d input' % (ci + 1))
    legval.append('shallow: CHL=%5.2f TSM=%5.2f DOC=%5.2f' % (c1[0], c1[1], c1[2]))
    #legval.append('deep: CHL=%5.2f TSM=%5.2f DOC=%5.2f' % (c1[0], c1[1], c1[2]))

x1,x2,y1,y2 = plt.axis()
plt.axis((x1, x2 * 1.1, y1, y2))
plt.xlabel('wavelength, nm')
plt.ylabel('Rrsw(-1), W m-2, nm')
plt.legend(legval)
plt.show()

