#!/usr/bin/python
from boreali import lm, Boreali

import numpy as np
import matplotlib.pyplot as plt

#from boreali import Boreali
from scipy import interpolate

parameters=[0.01, 15,
            0.01, 5,
            0.01, 5,]

bands = 6
#wavelen = [412, 443, 490, 510, 555, 670]
wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709]
b = Boreali('michigan', 'michigan')
# get matrix with HO-model
model = []
for i in range(0, 8):
    model.append(b.homo[i](wavelen))

#create albedo interpolators
albedoString = '''
lambda	400	425	450	475	500	525	550	575	600	625	650	675	700	725	750
sand	0.11	0.13	0.15	0.165	0.183	0.196	0.209	0.222	0.235	0.248	0.261	0.274	0.287	0.3	0.31
sargassum	0.02	0.02	0.02	0.02	0.02	0.03	0.04	0.06	0.07	0.05	0.05	0.02	0.15	0.31	0.35
silt	0.05	0.065	0.08	0.095	0.11	0.125	0.14	0.155	0.17	0.175	0.18	0.185	0.19	0.195	0.2
boodlea	0.05	0.05	0.06	0.07	0.11	0.25	0.28	0.25	0.22	0.19	0.15	0.08	0.25	0.43	0.45
limestone	0.065	0.075	0.09	0.1	0.115	0.13	0.14	0.15	0.145	0.14	0.14	0.125	0.135	0.1	0.07
enteromorpha	0.05	0.05	0.05	0.05	0.07	0.135	0.165	0.155	0.145	0.13	0.1	0.075	0.22	0.3	0.31
cladophora	0.05	0.05	0.05	0.05	0.07	0.135	0.165	0.155	0.145	0.13	0.1	0.075	0.14	0.175	0.18
chara	0.015	0.015	0.015	0.015	0.015	0.04	0.05	0.04	0.035	0.03	0.025	0.02	0.05	0.14	0.16
'''

albedoLines = albedoString.splitlines()
albedoNames = []
x = np.array(map(float, albedoLines[1].split('\t')[1:]))
albedos = []
for awi in range(2, len(albedoLines)):
    awiList = albedoLines[awi].split('\t')
    albedoNames.append(awiList[0])
    y = np.array(map(float, awiList[1:]))
    albedos.append(interpolate.interp1d(x, y))

#sun at typical zenith angle
theta = 25
