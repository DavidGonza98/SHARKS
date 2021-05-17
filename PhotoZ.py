# -*- coding: utf-8 -*-
"""
Created on Sun May 16 12:18:10 2021

@author: lpama
"""

class PhotoZ():
    """Class to characterize photo-z"""
    
    def __init__(self, codeName, parentList, tag):
        """
         Instantiate a PhotoZ.

        Parameters
        ----------
        codeName : `str`
            The name of the photoZ code, "eazy".
        parentList : `list`
             The list of files to analise: ['gals_dessharks_1.txt', 'gals_dessharks_2.txt', 'gals_dessharks_3.txt', 'gals_dessharks_4.txt', 'gals_dessharks_5.txt'].
        tag : `str`
            A tag associated with the parentList, "sharks".

        Returns
        -------
        None.

        """
        
        
        self.nombresalida = 'mytest'
        self.__Prepare()
        self.__mergeResults()
        
        
    def __Prepare(self):
        """
        Create the list of filenames we need.

        Returns
        -------
        None.

        """
        self.zout_results = []
        self.zero_points = []
        self.translate = []
        self.params = []
        self.datafits = []
        for f in self.parentList:
            self.zout_results.append(f + self.tag + '.zout.fits')
            self.zero_points.append(f + self.tag + '.zphot.zeropoint')
            self.translate.append(f + self.tag + '.zphot.translate')
            self.params.append(f + self.tag + '.zphot.param')
            self.datafits.append(f + self.tag + '.data.fits')
            
            
    def unify_data(self):
        import astropy.io.fits as pf
        import numpy as np

        files = ['gals_des_1.txt.eazypy.data.fits','gals_des_2.txt.eazypy.data.fits','gals_des_3.txt.eazypy.data.fits','gals_des_4.txt.eazypy.data.fits','gals_des_5.txt.eazypy.data.fits']
        hdul = pf.open(files[0])
        primary_hdu = hdul[0]

        zgrid_hdu = hdul[2]

        data_zbest = []
        for f in files:
            hda = pf.open(f)[1]
            for i in hda.data:
                data_zbest.append(i)
        data_zbest = np.array(data_zbest)
        zbest_hdu = pf.ImageHDU(data_zbest, name='ZBEST')

        data_chi = []
        for f in files:
            hda = pf.open(f)[3].data
            for i,vv in enumerate(hda):
                data_chi.append(vv)

        data_chi = np.array(data_chi)
        chi_hdu = pf.ImageHDU(data_chi, name='CHI2')

        data_coeffs = []
        for f in files:
            hda = pf.open(f)[4]
            for i in hda.data:
                data_coeffs.append(i)

        data_coeffs = np.array(data_coeffs)
        coeffs_hdu = pf.ImageHDU(data_zbest, name='COEFFS')

        hd_out = pf.HDUList([primary_hdu, zbest_hdu, zgrid_hdu, chi_hdu, coeffs_hdu])

        hd_out.writeto('test.data.fits')


        
    def __mergeResults(self):
        """
        This function merge results into a single file (see unify_data.py).

        Create a single file zout, parentcat, translate and zeropoints (only change name to first), change params (change name and set MAIN_OUTPUT_FILE and CATALOG_FILE), datafits.


        Returns
        -------
        None.

        """
        from astropy.table import Table, vstack
        from astropy.io import writeto
        lista_temp = []
        for f in  self.zout_results:
            
            t = Table.read(f, format= 'fits')
            lista_temp.append(t)
            
        self.mergeZout = vstack(lista_temp)
        self.mergeZout.writeto('Zout.fits')
        
    def getStats(self):
        import numpy as np  
        
        self.zspec_original=self.mergeZout['zspec']
        self.zspec=[]
        if np.greater_equal(self.zspec_original, 0.0):
            self.zspec.append(self.zspec_original)
            
        self.zphot= self.mergeZout['zphot']
        
        
        self.zbias= np.mean(self.delta_z_1pz())
        
        self.NMAD=self.arr/(1+self.zspec)
        self.outliers=self.arr
        
        self.t=self.sigma68/(1+self.NMAD)
        
        
    def delta_z(self):
        self.arr=(self.zspec-self.zphot)
        return self.arr
    
    def delta_z_1pz(self):
        return self.arr/(1+self.zspec)
    
    def sigma_68(arr, axis=None):
        """Input: an (multi-dimensional) array
           Optional input: the axis along which to calculate the metric
           Outputs: the 68% spread of data about the median value of the array
        """
        import numpy as np 
        upper, lower = np.percentile(arr, [84.075, 15.825], axis=axis)
        return (upper - lower) / 2.0

        
        
        
        
        
        
        

