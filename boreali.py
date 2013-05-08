import os
from time import time 
import inspect

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

import lm

class Boreali():
    '''Perform processing of optical image with the BOREALI boreali algorithm'''
    zones = {'ladoga': {'ho_model_file': 'models.txt',
                        'min': [0, 0, 0],
                        'max': [10, 5, 15],
                        'variables': ['lmchl', 'lmtsm', 'lmdoc', 'lmmse'],
                        },
            'northsea': {'ho_model_file': 'models.txt',
                        'min': [0, 0, 0],
                        'max': [5, 2, 2],
                        'variables': ['lmchl', 'lmtsm', 'lmdoc', 'lmmse'],
                        },
            'michigan': {'ho_model_file': 'models.txt',
                        'min': [0, 0, 0],
                        'max': [5, 2, 2],
                        'variables': ['lmchl', 'lmtsm', 'lmdoc', 'lmmse'],
                        },
            }
    def __init__(self, zone='ladoga', model='ladoga', wavelen=[412, 443, 490, 510, 555, 670]):
        '''Set up parameters'''
        self.z = self.zones[zone]
        self.homo = self.read_ho_models(self.z['ho_model_file'], model)
        self.set_albedos()
        self.wavelen = wavelen
    
    def set_albedos(self):
        '''Set spectra of albedos'''
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
        x = np.array(map(float, albedoLines[1].strip().split('\t')[1:]))
        albedos = []
        for awi in range(2, len(albedoLines)):
            awiList = albedoLines[awi].strip().split('\t')
            albedoNames.append(awiList[0])
            y = np.array(map(float, awiList[1:]))
            if len(y) == len(x):
                albedos.append(interpolate.interp1d(x, y))
        
        self.albedos = albedos
        self.albedoNames = albedoNames

    def get_homodel(self):
        '''Get matrix with hydro-optical model'''
        abCoef = []
        for i in range(0, 8):
            abCoef.append(self.homo[i](self.wavelen))
        
        return abCoef

    def get_albedo(self, bottom):
        '''Get array with albedos for each valid pixel'''
        #import pdb; pdb.set_trace()
        albedo = np.zeros((len(bottom), len(self.wavelen)))
        uniqBottom = np.unique(bottom)
        for ubt in uniqBottom:
            albedo[bottom == ubt] = self.albedos[ubt](self.wavelen) 
            
        return albedo
        
        
    def read_ho_models(self, hoFileName='', iModelName='ladoga'):
        '''Read file with HO-models and prepare interpolators'''
        
        # READ HO_MODELS FROM FILE
        selfDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))        
        f = open(os.path.join(selfDir, hoFileName))
        noOfModels = int(f.readline())
        hoModels = {}
        for i in range(0, noOfModels):
            fLine = f.readline().strip().split(' ')
            modelName = fLine[0]
            modelLines = int(fLine[1])
            hoModels[modelName] = np.zeros((modelLines, 7))
            for j in range(0, modelLines):
                fLine = f.readline().strip().split(' ')
                for k, x in enumerate(fLine):
                    hoModels[modelName][j, k] = float(x)

        # PREPARE INTERPOLATORS
        hoModelInterp = []
        for hoModelName in ['water', iModelName]:
            if hoModelName == 'water':
                abn = 2
            else:
                abn = 6
            #abn = {'water': 2}.get(hoModelName, 6)
            hoModel = hoModels[hoModelName]
            x = hoModel[:, 0]
            for yn in range(1, abn+1):
                hoModelInterp.append(interpolate.interp1d(x, hoModel[:, yn]))
            
        #shuffle crossections:
        #aW, bbW, aCHL, bbCHL, aTSM, bbTSM, aDOC, bbDOC => 
        #aW, aCHL, aTSM, aDOC, bbW, bbCHL, bbTSM, bbDOC
        hoModelInterp2 = hoModelInterp[:]
        hoModelInterp[0] = hoModelInterp2[0]
        hoModelInterp[1] = hoModelInterp2[2]
        hoModelInterp[2] = hoModelInterp2[4]
        hoModelInterp[3] = hoModelInterp2[6]
        hoModelInterp[4] = hoModelInterp2[1]
        hoModelInterp[5] = hoModelInterp2[3]
        hoModelInterp[6] = hoModelInterp2[5]
        hoModelInterp[7] = hoModelInterp2[7]

        return hoModelInterp
    
    def get_wavelengths(self, img, opts, wavelen=None):
        '''Return list of wavelengths from the available bands in img'''
        imgBands = img.bands()
        wavelengths = []
        
        bands = 0
        for imgBand in imgBands:
            bandName = imgBands[imgBand].get('name', '')
            if 'Rrsw_' in bandName or 'rrsw_' in bandName:
                wl = float(imgBands[imgBand]['wavelength'])
                if (wavelen is None or
                        (wavelen is not None and wl in wavelen)):
                    wavelengths.append(wl)
                    bands += 1
                if bands >= opts[0]:
                    break
        return wavelengths
    
    def get_rrsw(self, img, mask=None):
        '''Get array with RRSW for valid pixels and mask of valid pixels'''
        rrsw = []
        for wl in self.wavelen:
            bandName = 'Rrsw_%d' % wl
            bandArray = img[bandName]
            rrsw.append(bandArray[mask == 64])
        rrsw = np.array(rrsw).T
        
        return rrsw

    def process(self, img, opts,
                    mask=None,
                    depth=None,
                    bottom=None,
                    theta=None):
        '''Process Nansat object <img> with lm.cpp'''
        # get wavelengths
        self.wavelen = self.get_wavelengths(img, opts, self.wavelen)
        print 'wavelen', self.wavelen

        # assume all pixel are valid if mask is not given
        if mask is None:
            mask = np.zeros(img.shape(), 'uint8') + 64
        plt.imshow(mask);plt.title('mask');plt.colorbar();plt.show()

        # get Rrsw for valid pixels
        rrsw = self.get_rrsw(img, mask)
        pixels = rrsw.shape[0]

        # get bathymetry in valid pixels (or assume deep waters)
        if depth is None:
            depth = np.zeros(img.shape()) - 10
        plt.imshow(depth);plt.title('depth');plt.colorbar();plt.show()
        depth = depth[mask == 64]

        # get albedo for valid pixels (or assume sand bottom)
        if bottom is None:
            bottom = np.zeros(img.shape(), 'uint8')
        plt.imshow(bottom);plt.title('bottom');plt.colorbar();plt.show()
        albedo = self.get_albedo(bottom[mask == 64])

        # get solar zenith of valid pixels (or assume sun in zenith)
        if theta is None:
            theta = np.zeros(img.shape())
        plt.imshow(theta);plt.title('theta');plt.colorbar();plt.show()
        theta = theta[mask==64]

        abCoef = self.get_homodel(self.wavelen)

        # process RRSW spectra with LM
        t0 = time()
        print 'rrsw.shape:', rrsw.shape

        # retrieve N concetrations
        c = lm.get_c(opts, abCoef, rrsw, albedo, depth, theta, 4*pixels)[1]
        c = np.array(c).reshape(pixels, 4)
        print 'spent: ', time() - t0

        cImg = np.zeros((4, img.shape()[0], img.shape()[1]))
        for i in range(0, 4):
            cArray = np.zeros(img.shape())
            cArray[mask == 64] = c[:, i]
            cImg[i, :, :] = cArray
        
        return cImg
