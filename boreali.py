import os
from time import time 
import inspect
from multiprocessing import Process, Array

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

import lm

class Boreali():
    '''Perform processing of optical image with the BOREALI algorithm'''
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
        self.wavelen = wavelen
    
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

    def get_wavelengths(self, img, opts, wavelen=None):
        '''Return list of wavelengths from the available bands in img
        
        Open input Nansat image and loop through band names. Find for
        which wavelengths bands exist.
        
        Parameters:
        -----------
        img : Nansat
            input Nansat object with bands with names 'Rrsw_NNN'
        opts : list with options
            opts[0] should contain number of bands to read from file
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
                wl = float(imgBands[imgBand]['wavelength'])
                if (wavelen is None or
                        (wavelen is not None and wl in wavelen)):
                    wavelengths.append(wl)
                    bands += 1
                if bands >= opts[0]:
                    break
        return wavelengths
    
    def get_rrsw(self, img, flexible):
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

        Returns:
        --------
        2D Numpy array with Rrsw spectra (bands x pixels)
        '''
        wlNames = ['Rrsw_%d' % wl for wl in self.wavelen]
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
        if negative:
            negativePixels = rrsw.min(axis=0) < 0
            mask.flat[negativePixels * (mask.flat == 64)] = 4
        
        return rrsw[:, mask.flat == 64], mask
    
    def get_c(self, opts, abCoef,
                    rrsw, theta,
                    values):
        '''Retrieve CPA concentrations from Rrsw spectra using LM-optimization
        
        Wraper around lm.get_c() to be used in multi-processing
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
        theta : Numpy 1D array
            Sun zenith angle in each pixel.
        values : list
            container of output results
        '''
        valuesLen = len(values)
        c = lm.get_c(opts, abCoef, rrsw, theta, valuesLen)[1]
        values[:] = c
        

    def process(self, img, opts,
                    mask=None,
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
        opts : list with 9 values
            number of bands; number of starting vectors; min. RMSE;
            min/max of each CPA to retreive
        mask : Numpy 2D array (same shape as img)
            Nansat L2 mask with values 0, 1, 2, 64 meaning:
            out_of_swath, cloud, land, valid_pixel
            If None, all pixels are considered valid
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
        self.wavelen = self.get_wavelengths(img, opts, self.wavelen)
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

        # get Rrsw for all pixels [bands x pixels]
        rrsw = self.get_rrsw(img, flexible)
        
        # get Rrsw for valid pixels [bands x pixels] and new mask
        rrsw, mask = self.filter(rrsw, mask)
        
        # transpose Rrsw and get number of pixels
        rrsw = rrsw.T
        pixels = rrsw.shape[0]

        # get solar zenith of valid pixels (or assume sun in zenith)
        if theta is None:
            theta = np.zeros(img.shape())
        #plt.imshow(theta);plt.title('theta');plt.colorbar();plt.show()
        theta = theta[mask==64]

        abCoef = self.get_homodel()

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
            # subset rrsw and theta
            r = rrsw[thIndex[thn][0]:thIndex[thn][1], :]
            t = theta[thIndex[thn][0]:thIndex[thn][1]]
            # create new subprocess args are:
            # name of func to retrieve CPA and
            # arguments for tht func
            procs.append(Process(target=self.get_c,
                                        args=(opts, abCoef,
                                                r, t,
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
        
        return cImg, mask
