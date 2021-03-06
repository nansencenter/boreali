import os
from time import time
import inspect
from multiprocessing import Process, Array
import string
import re

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

import lm

class Boreali():
    '''Perform processing of optical image with the BOREALI algorithm'''
    B_WAT = 0.5
    B_CHL = 0.0011
    B_TSM = 0.08
    B_DOC = 1

    KD0 = -0.218
    KD1 =  0.473

    mu0 = 1
    def __init__(self, model='ladoga', wavelen=[412, 443, 490, 510, 555, 670]):
        '''Generate Boreali and Set up retrieval options

        Parameters:
        -----------
        model : string, name of the model in the models.txt
        wavelen : list of int, wavelength of bands

        Modifies:
        ---------
        Reads hydropotical model from given name and or given wavelengths
        Read albedos for given wavelengths
        '''
        self.homo = self.read_ho_models('models.txt', model)
        self.set_albedos()
        self.wavelen = wavelen
        self.model = self.get_homodel()

    def read_ho_models(self, hoFileName='', iModelName='ladoga'):
        '''Read file with HO-models and prepare interpolators

        Parameters:
        -----------
            hoFileName : string
                name of the text file with hydrooptical model
            iModelName : string
                name of the ho-model in the file
        Returns:
        --------
            interpolators for each a/bb component. Eah interpolator can
            be applied to retrieve a/bb for any wavelength
            '''

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

    def set_albedos(self):
        '''Set spectra of albedos

        Spectra of albedos are kept in this method in a tring. Albedo
        values are read from string and a list of interpolators
        which relate wavelengths and albed is created for each bottom

        Modifies:
        ---------
        self.albedos : list of interpolation functions for each albedo
        self.albedoName : names of albedos
        '''

        albedoString = '''
        lambda  400 425 450 475 500 525 550 575 600 625 650 675 700 725 750
        00sand  0.11    0.13    0.15    0.165   0.183   0.196   0.209   0.222   0.235   0.248   0.261   0.274   0.287   0.3 0.31
        01sargassum 0.02    0.02    0.02    0.02    0.02    0.03    0.04    0.06    0.07    0.05    0.05    0.02    0.15    0.31    0.35
        02silt  0.05    0.065   0.08    0.095   0.11    0.125   0.14    0.155   0.17    0.175   0.18    0.185   0.19    0.195   0.2
        03boodlea   0.05    0.05    0.06    0.07    0.11    0.25    0.28    0.25    0.22    0.19    0.15    0.08    0.25    0.43    0.45
        04limestone 0.065   0.075   0.09    0.1 0.115   0.13    0.14    0.15    0.145   0.14    0.14    0.125   0.135   0.1 0.07
        05enteromorpha  0.05    0.05    0.05    0.05    0.07    0.135   0.165   0.155   0.145   0.13    0.1 0.075   0.22    0.3 0.31
        06cladophora    0.05    0.05    0.05    0.05    0.07    0.135   0.165   0.155   0.145   0.13    0.1 0.075   0.14    0.175   0.18
        07chara 0.015   0.015   0.015   0.015   0.015   0.04    0.05    0.04    0.035   0.03    0.025   0.02    0.05    0.14    0.16
        08sand2 0.20    0.25    0.32    0.33    0.35    0.4 0.63    0.64    0.66    0.67    0.69    0.71    0.72    0.73    0.74
        09sand3 0.22    0.26    0.30    0.33    0.27    0.40    0.41    0.444   0.47    0.49    0.52    0.54    0.57    0.6 0.62
        10charasand 0.0625  0.0725  0.0825  0.09    0.099   0.118   0.1295  0.131   0.135   0.139   0.143   0.147   0.1685  0.22    0.235
        11cladforaosand 0.09    0.1033333333    0.1166666667    0.1266666667    0.1453333333    0.1756666667    0.1943333333    0.1996666667    0.205   0.2086666667    0.2073333333    0.2076666667    0.238   0.2583333333    0.2666666667
        12sandmichi 0.06    0.0627907   0.06529778  0.06707556  0.07053953  0.07583023  0.08794583  0.09242222  0.09445 0.09647778  0.09749394  0.09446364  0.09590278  0.09795139  0.1
        13cladomichi    0.025   0.0277907   0.03020444  0.03164889  0.03424186  0.03807907  0.04699583  0.05014444  0.05145 0.05275556  0.05317879  0.05007273  0.05131944  0.05315972  0.055

        '''
        #sand   0.11    0.13    0.15    0.165   0.183   0.196   0.209   0.222   0.235   0.248   0.261   0.274   0.287   0.3 0.31
        #sand        0.2    0.25    0.32    0.33    0.35    0.4 0.63    0.64    0.66    0.67    0.69    0.71    0.72    0.73    0.74

        # split the string into line
        albedoLists = []
        for al in albedoString.splitlines():
            if len(al.strip()) > 0:
                alSingleSpace = re.sub(' +', ' ', al.strip())
                alList = alSingleSpace.replace('\t', ' ').split(' ')
                albedoLists.append(alList)

        x = np.array(map(float, albedoLists[0][1:]))
        albedoNames = []
        albedos = []
        for alist in albedoLists[1:]:
            albedoNames.append(alist[0])
            y = np.array(map(float, alist[1:]))
            if len(y) == len(x):
                albedos.append(interpolate.interp1d(x, y))

        self.albedos = albedos
        self.albedoNames = albedoNames


    def get_homodel(self):
        '''Get matrix with hydro-optical model for curent wavelength

        Access the list of HO-model interpolators and create a matrix
        with a and bb coefficients for the curent wavelengths

        Returns:
        --------
        Numpy 2D array : Matrix with values of
        aW, aCHL, aTSM, aDOC, bbW, bbCHL, bbTSM, bbDOC
        for each wavelength.
        '''

        abCoef = []
        for i in range(0, 8):
            abCoef.append(self.homo[i](self.wavelen))

        return abCoef

    def get_albedo(self, bottom):
        '''Get matrix with albedos for each input pixel

        Each albedo spectrum is calculated from albedo interpolators
        using the current wavelengths and input index of bottom type

        Parameters:
        -----------
        bottom : numpy 1D integer matrix
            array with indeces of bottom types

        Returns:
        ---------
        Numpy 2d matrix: matrix with albedo spectra for each input index
        '''

        #import pdb; pdb.set_trace()
        albedo = np.zeros((len(bottom), len(self.wavelen)))
        uniqBottom = np.unique(bottom)
        for ubt in uniqBottom:
            albedo[bottom == ubt] = self.albedos[ubt](self.wavelen)

        return albedo

    def get_wavelengths(self, img, wavelen=None):
        '''Return list of wavelengths from the available bands in img

        Open input Nansat image and loop through band names. Find for
        which wavelengths bands exist.

        Parameters:
        -----------
        img : Nansat
            input Nansat object with bands with names 'Rrsw_NNN'
        wavelen : list of integers
            list of wavelengths (in nm) to find in the input file.
            If empty - all Rrsw_ bands are read

        Returns:
        --------
        list of int, wavelenghs available in the file

        '''
        imgBands = img.bands()
        wavelengths = []

        bands = 0
        for imgBand in imgBands:
            bandName = imgBands[imgBand].get('name', '')
            if 'Rrsw_' in bandName or 'rrsw_' in bandName:
                if 'wavelength' not in imgBands[imgBand]:
                    continue
                wl = float(imgBands[imgBand]['wavelength'])
                if (wavelen is None or
                        (wavelen is not None and wl in wavelen)):
                    wavelengths.append(wl)
                    bands += 1
        return wavelengths

    def get_rrsw(self, img, flexible, osw):
        '''Get spectra of Rrsw from input image

        Get bands with Rrsw from input image and reshape into 2D array
        with Rrsw spectra (bands x pixels)

        Parameters:
        -----------
        img : Nansat
            input Nansat image
        flexible : boolean
            if True also bands with name Rrsw_NNNZZZ are read where
            NNN is wavelenth and ZZZ are any symbols
        osw : boolean
            is it for optically shallow waters?

        Returns:
        --------
        2D Numpy array with Rrsw spectra (bands x pixels)
        '''
        if osw:
            rName = 'Rrs_'
        else:
            rName = 'Rrsw_'
        print 'R Name is ', rName
        wlNames = ['%s%d' % (rName, wl) for wl in self.wavelen]
        # get names of bands with Rrsw
        bands = img.bands()
        rrsw = []
        for bn in bands:
            bandName = bands[bn]['name']
            for wlName in wlNames:
                if bandName == wlName or (flexible and wlName in bandName):
                    bandArray = img[bandName]
                    rrsw.append(bandArray)

        rrsw = np.array(rrsw)
        rrsw = rrsw.reshape(rrsw.shape[0], rrsw.shape[1] * rrsw.shape[2])

        return rrsw

    def filter(self, rrsw, mask, negative=True):
        '''Remove invalid Rrsw spectra using mask and custom filters

        Search for negative values in Rrsw and if found mask this pixel
        with value=4 in the mask

        Parameters:
        -----------
            rrsw : 2D Numpy array with Rrsw spectra
            mask : 2D Numpy array with L2 Nansat mask
            negative : boolean, apply negative Rrsw filter?

        Returns:
        --------
            rrsw : 2D Numpy array with valid Rrsw spectra only
            mask : 2D Numpy array with updated L2 mask (Boreali mask)
        '''
        # apply 'negative values' filter to water pixels
        resMask = np.array(mask)
        if negative:
            negativePixels = rrsw.min(axis=0) < 0
            resMask.flat[negativePixels * (resMask.flat == 64)] = 4

        return rrsw[:, resMask.flat == 64], resMask

    def get_c(self, opts, abCoef, osw,
                    rrsw, albedo, depth, theta,
                    values):
        '''Retrieve CPA concentrations from Rrsw spectra using LM-optimization

        Wraper around lm.get_c_shal() to be used in multi-processing
        Parameters:
        -----------
        opts : list with 9 values
            number of bands; number of starting vectors; min. RMSE;
            min/max of each CPA to retreive
        abCoef :  Numpy 2D matrix
            Values of absorption and bacscattering og water and CPA for
            each wavelength
        rrsw : numpy 2D array
            Rrsw spectra
        albedo : numpy 2D array
            bottom albedo spectra
        depth : Numpy 1D array
            Bathymetry of each pixel. If negative - shallow water
            retrieval is not applied
        theta : Numpy 1D array
            Sun zenith angle in each pixel.
        values : list
            container of output results
        '''
        valuesLen = len(values)
        if osw:
            print 'launch shallow'
            c = lm.get_c_shal(opts, abCoef, rrsw, theta, depth, albedo, valuesLen)[1]
        else:
            print 'launch deep'
            c = lm.get_c_deep(opts, abCoef, rrsw, theta, valuesLen)[1]
        values[:] = c

    def process(self, img, opts,
                    mask=None,
                    depth=None,
                    bottom=None,
                    theta=None,
                    flexible=True,
                    threads=1):
        '''Retrieve CPA concentrations from input Nansat image

        Read Rrsw values from bands of input Nansat object,
        filter only valid pixels, apply Levenberg-Marquardt optimization
        for retrieval of CPA concentrations from Rrsw spectra

        Parameters:
        -----------
        img : Nansat
            input image with Rrsw bands. The bands should have name
            'Rrsw_NNN' where NNN is an integer value of wavelength in nm
            NNN should match the wavelengths from the self.wavelen
            created when initializing Boreali object
        opts : list with 6 values
            min/max of each CPA to retreive
        mask : Numpy 2D array (same shape as img)
            Nansat L2 mask with values 0, 1, 2, 64 meaning:
            out_of_swath, cloud, land, valid_pixel
            If None, all pixels are considered valid
        depth : Numpy 2D array (same shape as img)
            Bathymetry of the ROI. If none - shallow water retrieval is not applied
        bottom : Numpy 2D array (same shape as img)
            Type of bottom cover. Matrix with integer values of index of
            bottom type. See Boreali.get_albedo()
            If None, sand bottom is assumed
        theta : Numpy 2D array (same shape as img)
            Sun zenith angle in each pixel. If None sun is in nadir
        flexible : Boolean
            Allow the name of the band to be 'Rrsw_NNNZZZZ' where ZZZ
            is any number. Required for ability to process NC-files
            previously saved by Nansat
        threads : int
            Number of subprocesses

        Returns:
        --------
            cpa : Numpy 3D array (4 x (same shape as img))
                Matrix with CPA and RMSE. CPA are the agents from the
                hydro-opticalmodel. RMSE is the retrieval error (diff
                between the measured and reconstructed spectra)
            mask : Numpy 2D array (same shape as img)
                Boreali L2 Mask. In addition to Nansat mask value 4
                designates negative Rrsw values
        '''
        # get wavelengths
        self.wavelen = self.get_wavelengths(img, self.wavelen)
        print 'wavelen', self.wavelen

        # assume all pixel are valid if mask is not given
        if mask is None:
            try:
                mask = img['mask']
            except:
                mask = None
        if mask is None:
            mask = np.zeros(img.shape(), 'uint8') + 64
        #plt.imshow(mask);plt.title('mask');plt.colorbar();plt.show()

        # get bathymetry in valid pixels (or assume deep waters)
        osw = True
        if depth is None:
            osw = False
            depth = np.zeros(img.shape()) + 1000

        # get Rrsw (or Rrs) for all pixels [bands x pixels]
        rrsw = self.get_rrsw(img, flexible, osw)

        # get Rrsw for valid pixels [bands x pixels] and new mask
        rrsw, mask = self.filter(rrsw, mask)

        # transpose Rrsw and get number of pixels
        rrsw = rrsw.T
        pixels = rrsw.shape[0]

        # keep only good pixels
        depth = depth[mask == 64]

        # get albedo for valid pixels (or assume sand bottom)
        if bottom is None:
            bottom = np.zeros(img.shape(), 'uint8')
        #plt.imshow(bottom);plt.title('bottom');plt.colorbar();plt.show()
        albedo = self.get_albedo(bottom[mask == 64])

        # get solar zenith of valid pixels (or assume sun in zenith)
        if theta is None:
            theta = np.zeros(img.shape())
        #plt.imshow(theta);plt.title('theta');plt.colorbar();plt.show()
        theta = theta[mask==64]

        #abCoef = self.get_homodel()
        abCoef = self.model

        # prepare for multiprocessing: create last and first indices for
        # dividing input matrix with rrsw into sub-matrices
        thIndex = []
        thi0 = 0
        for thn in range(threads):
            thIndex.append([thi0, thi0+pixels/threads])
            thi0 += pixels/threads
        thIndex[-1][-1] = pixels

        t0 = time()
        print 'Pixels x bands to process: ', rrsw.shape

        # create and launch subprocesss (collect in procs)
        # keep results in cVals
        procs = []
        cVals = []
        for thn in range(threads):
            # select sub-matrix for this subprocess
            procPixels = thIndex[thn][1] - thIndex[thn][0]
            # allocate memory to store results of rerieval from submatrix
            cVals.append(Array('d', range(procPixels*4)))
            # subset rrsw, albedo, depth and theta
            r = rrsw[thIndex[thn][0]:thIndex[thn][1], :]
            a = albedo[thIndex[thn][0]:thIndex[thn][1], :]
            h = depth[thIndex[thn][0]:thIndex[thn][1]]
            t = theta[thIndex[thn][0]:thIndex[thn][1]]
            # create new subprocess args are:
            # name of func to retrieve CPA and
            # arguments for tht func
            procs.append(Process(target=self.get_c,
                                        args=(opts, abCoef, osw,
                                                r, a, h, t,
                                                cVals[thn])))
            # launch sub-process
            procs[thn].start()

        # wait for subprocess to finish
        for thn in range(threads):
            procs[thn].join()

        # collect retrieval results into one list
        c = []
        for thn in range(threads):
            c += cVals[thn]

        # convert list with retrievals into 2d matrix with column for each CPA
        c = np.array(c).reshape(pixels, 4)
        print 'spent: ', time() - t0

        # reshape 2D-CPA matix into 3D with 2D images of CPA
        cImg = np.zeros((4, img.shape()[0], img.shape()[1])) + np.nan
        for i in range(0, 4):
            cArray = np.zeros(img.shape()) + np.nan
            cArray[mask == 64] = c[:, i]
            cImg[i, :, :] = cArray

        return np.vstack([cImg, mask[None]])

    def get_kd(self, c):
        '''Get spectral Kd'''
        #aW, aCHL, aTSM, aDOC, bbW, bbCHL, bbTSM, bbDOC
        awat = self.model[0]
        achl = self.model[1]
        atsm = self.model[2]
        adoc = self.model[3]
        bbwat = self.model[4]
        bbchl = self.model[5]
        bbtsm = self.model[6]
        bbdoc = self.model[7]
        bwat = bbwat / self.B_WAT
        bchl = bbchl / self.B_CHL
        btsm = bbtsm / self.B_TSM
        bdoc = 0

        a = awat + achl * c[0] + atsm * c[1] + adoc * c[2]
        bb = bbwat + bbchl * c[0] + bbtsm * c[1] + bbdoc * c[2]
        b = bwat + bchl * c[0] + btsm * c[1] + bdoc * c[2]
        kd = np.sqrt(a ** 2 + a * b * (self.KD0 + self.KD1 * self.mu0)) / self.mu0

        return kd

    def get_attenuated_ed(self, ed0, kd, h):
        '''Calculate attenuated Ed'''

        ed = ed0 * np.exp(- kd * h)
        return ed
