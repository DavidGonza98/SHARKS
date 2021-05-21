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
        self.codeName=codeName
        self.parentList=parentList
        self.tag=tag
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
            self.zout_results.append(f +'.'+ self.tag + '.zout.fits')
            self.zero_points.append(f + '.'+self.tag + '.zphot.zeropoint')
            self.translate.append(f +'.'+ self.tag + '.zphot.translate')
            self.params.append(f +'.'+ self.tag + '.zphot.param')
            self.datafits.append(f +'.'+ self.tag + '.data.fits')
            
            
    def unify_data(self):
        import astropy.io.fits as pf
        import numpy as np

        files = self.datafits
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

        hd_out.writeto(self.nombresalida+'.data.fits')


        
    def __mergeResults(self):
        """
        This function merge results into a single file (see unify_data.py).

        Create a single file zout, parentcat, translate and zeropoints (only change name to first), change params (change name and set MAIN_OUTPUT_FILE and CATALOG_FILE), datafits.


        Returns
        -------
        None.

        """
        from astropy.table import Table, vstack
        import astropy.io.fits
        lista_temp = []
        for f in  self.zout_results:
            
            t = Table.read(f, format= 'fits')
            lista_temp.append(t)
            
        self.mergeZout = vstack(lista_temp)
        self.mergeZout.write(self.nombresalida + 'Zout.fits', overwrite=True)
        
    def getStats(self):
        import numpy as np 
        from astropy.table import Table, vstack
        import astropy.io.fits

 
        self.zspec_original= self.mergeZout['z_spec']
        mask = (self.zspec_original>0)
        self.zspec = self.zspec_original[mask]
        self.zphot= self.mergeZout['z_phot'][mask]
        '''
        self.zspec=[]
        if np.greater_equal(self.zspec_original, 0.0) :
            if self.zspec_original.all()!=-1:
                self.zspec.append(self.zspec_original.all())
        self.zspec=np.array(self.zspec)
        print(self.zspec)
        self.zphot= self.mergeZout['z_phot']
        '''
        self.arr=(self.zspec-self.zphot)
        self.zbias= np.mean(self.delta_z_1pz())
        
        self.NMAD = np.mean(np.abs(self.arr/(1+self.zspec))) 
        self.outliers=self.outlier_rate(self.arr)
        
        arr = self.arr/(1+self.zspec)
        
        self.t=self.sigma_68(arr)
        
        print('zbias',self.zbias)
        print('Outliers', self.outliers)
        print ('sigma_68', self.t)
        print('NMAD', self.NMAD)
        
    def getStats_Bin(self, lims):
        
        
        import numpy as np
        import matplotlib.pyplot as plt
        
        #lims = [0.,0.1,0.2,0.3,0.4,0.5]
        self.zspec_original= self.mergeZout['z_spec']
        self.zphot= self.mergeZout['z_phot']

        bin_central = []
        sigma68 = []
        outliers = []
        zbias = []
        NMAD = []
        for i in range(len(lims)-1):
            mask = (self.zspec_original>lims[i])&(self.zspec_original<lims[i+1])

            zspec_tem = self.zspec_original[mask]
            zphot_tem = self.zphot[mask]
            bin_central.append((lims[i]+lims[i+1])/2.)
            arr = (zspec_tem-zphot_tem)/(1+zspec_tem)
            zbias.append(np.mean(zspec_tem-zphot_tem))
            outliers.append(self.outlier_rate(zspec_tem-zphot_tem))
            sigma68.append(self.sigma_68(arr))
            NMAD.append(np.mean(np.abs((zspec_tem-zphot_tem)/(1+zspec_tem))))
            #obtienes zbias, sigma68 y outliers
        
        # 1 plot para metrica, para zbias, sigma68 y outliers, donde eje x = bin_central, y = metrica
        plt.plot(bin_central, zbias, label='zbias')
        plt.plot(bin_central, sigma68, label='sigma68')
        plt.plot(bin_central, outliers, label='outliers')
        plt.plot(bin_central, NMAD, label='NMAD')
        plt.xlabel('zspec') 
        plt.ylabel('zphot')
        plt.title('Metrical measures')
        plt.legend()
        plt.savefig('Metrical_measures.png')
            
    def delta_z(self):
        
        return self.arr
    
    def delta_z_1pz(self):
        return self.arr/(1+self.zspec)
    
    def sigma_68(self, arr, axis=None):
        """Input: an (multi-dimensional) array
           Optional input: the axis along which to calculate the metric
           Outputs: the 68% spread of data about the median value of the array
        """
        import numpy as np 
        upper, lower = np.percentile(arr, [84.075, 15.825], axis=axis)
        return (upper - lower) / 2.0
    
    def outlier_rate(self, arr, outR=None):
        """assumes frac outliers >0.15"""
        import numpy as np
    
        if outR is None:
            outR = 0.15
        return np.sum(np.abs(arr) > outR)*1.0/len(arr)
        
        
'''        
results_eazy = PhotoZ('eazy', ['gals_dessharks_1.txt', 'gals_dessharks_2.txt', 'gals_dessharks_3.txt', 'gals_dessharks_4.txt', 'gals_dessharks_5.txt'], 'sharks')

results_eazy.getStats()  
results_eazy.getStats_Bin([0.,0.1,0.2,0.3,0.4,0.5])    
#results_eazy.delta_z_1pz()        
'''        
        

