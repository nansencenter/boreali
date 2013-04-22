import os, os.path
from datetime import datetime, timedelta
import numpy as np
from xml.etree.ElementTree import XML, ElementTree, tostring

import matplotlib.pyplot as plt

import gdal

from nansat import Domain, Nansat, Figure

#from boreali import Boreali


class OBPGL2Image(Nansat):

    def l2c_mask(self, cloudBits, landBits=[2], tmpProdName='latitude', invalidTmp=-999):
        '''Create l2c_flags:
        Flag coding:
        1 - cloud   (value = 1)
        2 - land    (value = 2)
        8 - water   (value = 64)

        L2 BITS:
        f01_name=ATMFAIL            #
        f02_name=LAND
        f03_name=PRODWARN           
        f04_name=HIGLINT            #
        f05_name=HILT               #
        f06_name=HISATZEN           #
        f07_name=COASTZ
        f08_name=SPARE
        f09_name=STRAYLIGHT         .#
        f10_name=CLDICE             #
        f11_name=COCCOLITH          .
        f12_name=TURBIDW            .
        f13_name=HISOLZEN           #
        f14_name=SPARE
        f15_name=LOWLW              #
        f16_name=CHLFAIL
        f17_name=NAVWARN
        f18_name=ABSAER
        f19_name=SPARE
        f20_name=MAXAERITER         #
        f21_name=MODGLINT           .#
        f22_name=CHLWARN
        f23_name=ATMWARN            #
        f24_name=SPARE
        f25_name=SEAICE
        f26_name=NAVFAIL
        f27_name=FILTER
        f28_name=SSTWARN            .
        f29_name=SSTFAIL			#
        f30_name=HIPOL              #
        f31_name=PRODFAIL           .
        f32_name=SPARE
        '''

        # get l2_flags (try to change WorkingDataType if VRT is warped)
        try:
            self._modify_warpedVRT2('GDALWarpOptions/WorkingDataType', 'Int32')
        except:
            self.logger.error('Unable to modify WorkingDataType in VRT!')
        
        l2_flags = self['flags']
        
        try:
            self._modify_warpedVRT2('GDALWarpOptions/WorkingDataType', 'Float32')
        except:
            self.logger.error('Unable to modify WorkingDataType in VRT!')

        tmpVar = self[tmpProdName]
                
        # == FOR DEBUG ==
        #import matplotlib.pyplot as plt
        #plt.imshow(l2_flags);plt.colorbar();plt.show()
        #for bit in range(0, 32):
        #    maskTmp = np.bitwise_and(l2_flags, np.power(np.uint32(2), bit))
        #    fName = '%s_%02d.png' % (self.name, bit)
        #    print fName
        #    plt.imsave(fName, maskTmp)
        
        # create l2c_falg matrix with water bit by default
        l2c_mask = 64 * np.ones(l2_flags.shape, np.uint8)
        
        # process cloud and land masks
        maskLand = np.zeros(l2c_mask.shape)
        maskCloud = np.zeros(l2c_mask.shape)
        # check every bit in a set
        for bit in landBits:
            maskTmp = np.bitwise_and(l2_flags, np.power(np.uint32(2), bit-1))
            maskLand[maskTmp > 0] = 1
        # check every bit in a set
        for bit in cloudBits:
            maskTmp = np.bitwise_and(l2_flags, np.power(np.uint32(2), bit-1))
            maskCloud[maskTmp > 0] = 1
        
        # set cloud bit
        l2c_mask[maskCloud > 0] = 1
        # set land bit
        l2c_mask[maskLand > 0] = 2
        # erase data out of swath
        l2c_mask[tmpVar == 0] = 0
        # erase bad values where latitutde is invalid
        l2c_mask[tmpVar == invalidTmp] = 0
        
        return l2c_mask

    def _modify_warpedVRT2(self, key, value):
        ''' Modify workingDataType in the warped VRT

        Parameters
        ----------
            workingDataType: string
                desired WorkingDataType 
        Modifies
        --------
            the VRT file which keepes warped vrt is modified
        '''
        # Get XML content from VSI-file
        vsiFileContent = self.vrt.read_xml()

        # Get element from the XML content and modify it
        element = XML(vsiFileContent)
        tree = ElementTree(element)
        elem = tree.find(key)
        elem.text = value

        # Overwrite element
        element = tree.getroot()

        # Write the modified elemements into VSI-file
        self.vrt.write_xml(tostring(element))
