#!/usr/bin/python
from boreali import lm, Boreali

import numpy as np
import matplotlib.pyplot as plt

import time

parameters=[0.01, 5,
            0.01, .5,
            0.01, .2,]

wavelen = [412, 443, 490, 510, 555, 670]
#wavelen = [413, 443, 490, 510, 560, 620, 665, 681, 709]

#v0vector = [0, 1, 2, 5, 10] # r-noise [%]
#xlabel = 'noise level in reflectance, %'
#v0Name = 'rnoise_'


v0vector = [0, 0.2, 0.5, 1., 2.] # h-noise [meters]
xlabel = 'noise level in depth, m'
v0Name = 'hnoise_'


#v0vector = [0, 1, 2, 5, 10] # a-noise [%]
#xlabel = 'noise level in albedo, %'
#v0Name = 'anoise_'

v1vector = [5, 10] # depth

v2vector = [0, 2, 4, 6, 7] # bottom types sand, silt, limestone, cladophora, chara


raTitles = ['relative', 'absolute']
varTitles = ['CHL', 'TSM', 'aCDOM']
v1Titles = ['h=5 ', 'h=10 ']
v2Titles = ['sand', 'silt', 'lime', 'clad', 'char']

v1Lines = ['.-', '-.']
v2Lines = ['r', 'g', 'b', 'c', 'm']
