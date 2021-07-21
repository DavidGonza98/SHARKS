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
        parentList : `list` or `str`
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
        self.nombresalida = self.tag
        
        if self.codeName=="eazy":
            self.__Prepare()
            self.__mergeResults()
            
        elif self.codeName=="lephare":
            self.__readLephare()
        
        elif self.codeName=="lephare_master_cat":
            self.__readMaster_Cat()
            
            
        else:
            print('Introduzca un codeName correcto (eazy, lephare o lephare_master_cat)')
        
        
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
            
            
    def __readLephare(self):
        import numpy as np
        
        z_phot, z_spec = np.loadtxt(self.parentList, usecols=(1,41), comments='#', unpack=True)
        self.mergeZout= {'z_phot':z_phot, 'z_spec':z_spec}
        
    def __readMaster_Cat(self):
        from astropy.table import Table
        if self.tag=="Master_Cat":
            self.datos = Table.read(self.parentList, format= 'fits')
            self.z_spec= self.datos["ZSPEC"]
            #mask= (self.z_spec!=1)
            self.datos=self.datos#[mask]
            self.zspec=self.z_spec#[mask]
            self.zphot_doy= self.datos["ZPHOT_DOY"]#[mask]
            self.zphot_shd=self.datos["ZPHOT_SHD"]#[mask]
            self.zphot_cww_doy=self.datos["ZPHOT_CWW_DOY"]#[mask]
            self.zphot_cww_shd=self.datos["ZPHOT_CWW_SHD"]#[mask]
        
        elif self.tag=="Master_Cat_bpz":
            self.datos = Table.read(self.parentList, format= 'fits')
            self.z_spec= self.datos["ZSPEC"]
            #mask= (self.datos["TYPE_SHD"]=="G")&(self.datos["TYPE_DOY"]=="G")
            self.datos=self.datos#[mask]
            self.zspec=self.z_spec#[mask]
            self.z_best_doy= self.datos["Z_BEST_DOY"]#[mask]
            self.z_best_shd=self.datos["Z_BEST_SHD"]#[mask]
            self.mag_z=self.datos["MAG_Z"]
            self.mag_i=self.datos["MAG_I"]
            self.mag_Y=self.datos["MAG_Y"]
            self.mag_Ks=self.datos["MAG_KS"]
    
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
        
    def getStats_Bin(self, lims, split=False):
        
        
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns
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
        if split==False:
            plt.clf()
            plt.figure(1)
            plt.plot(bin_central, zbias, label='zbias')
            plt.plot(bin_central, sigma68, label='sigma68')
            plt.plot(bin_central, outliers, label='outliers')
            plt.plot(bin_central, NMAD, label='NMAD')
            plt.xlabel('zspec') 
            plt.ylabel('zphot')
            plt.title('Metrical measures')
            plt.legend()
            plt.savefig('%s_Metrical_measures.png' % self.tag)
            plt.figure(2)
            sns.displot(zphot_tem, kde=True)
            plt.title('Zphot_'+self.tag)
            plt.xlabel('Zphot')
            plt.savefig('Zphot_'+self.tag+'.png')
            
        else:
            plt.figure(1)
            plt.plot(bin_central, zbias, label='zbias_'+self.tag)
            plt.xlabel('zspec') 
            plt.ylabel('zphot')
            plt.title('Metrical measures')
            plt.legend()
            plt.savefig('Zbias_measures.png')
            plt.figure(2)
            plt.plot(bin_central, sigma68, label='sigma68_'+self.tag)
            plt.xlabel('zspec') 
            plt.ylabel('zphot')
            plt.title('Metrical measures')
            plt.legend()
            plt.savefig('Sigma68_measures.png')
            plt.figure(3)
            plt.plot(bin_central, outliers, label='outliers_'+self.tag)
            plt.xlabel('zspec') 
            plt.ylabel('zphot')
            plt.title('Metrical measures')
            plt.legend()
            plt.savefig('Outliers_measures.png')
            plt.figure(4)
            plt.plot(bin_central, NMAD, label='NMAD_'+self.tag)
            plt.xlabel('zspec') 
            plt.ylabel('zphot')
            plt.title('Metrical measures')
            plt.legend()
            plt.savefig('NMAD_measures.png')
        
   
        
        
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
    
    def plot(self):
        import matplotlib.pyplot as plt
        mask=(self.z_spec!=-1)&(self.z_spec!=1)
        z_spec=self.z_spec[mask]
        zphot_doy=self.zphot_doy[mask]
        zphot_shd=self.zphot_shd[mask]
        zphot_cww_doy=self.zphot_cww_doy[mask]
        zphot_cww_shd=self.zphot_cww_shd[mask]

        plt.clf()
        fig, axs = plt.subplots(2, 2)
        axs[0, 0].plot(z_spec, zphot_doy,'.', label='Zphot_doy')
        axs[0, 0].plot(z_spec, z_spec, 'r--')
        axs[0, 0].set_title('DOY')
        axs[0, 0].legend()
        axs[0, 1].plot(z_spec, zphot_shd,'.', label='Zphot_shd')
        axs[0, 1].plot(z_spec, z_spec, 'r--')
        axs[0, 1].set_title('SHD')
        axs[0, 1].legend()
        axs[1, 0].plot(z_spec, zphot_cww_doy,'.', label='Zphot_cww_doy')
        axs[1, 0].plot(z_spec, z_spec, 'r--')
        axs[1, 0].set_title('CWW_DOY')
        axs[1, 0].legend()
        axs[1, 1].plot(z_spec, zphot_cww_shd,'.', label='Zphot_cww_shd')
        axs[1, 1].plot(z_spec, z_spec, 'r--')
        axs[1, 1].set_title('CWW_SHD')
        axs[1, 1].legend()
        
        for ax in axs.flat:
            ax.set(xlabel='zspec', ylabel='zphot')
            ax.set(xlim=[0., 1], ylim=[0, 1])

        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()
        plt.savefig('Subplot_Zphot.png') 
        plt.clf()
    def get_stats_Master(self, lims, without_bin=False, phot_type=False):
        
        import matplotlib.pyplot as plt
        import numpy as np
        import seaborn as sns
        
        bin_central = []
        sigma68 = []
        outliers = []
        zbias = []
        NMAD = []
        
        self.phot_type= self.datos['PHOT_TYPE_SHD']
        mask1 = (self.phot_type=='G')&(self.zspec>0)&(self.zspec!=-1)&(self.zspec!=1)

        lista=[self.zphot_doy, self.zphot_shd, self.zphot_cww_doy, self.zphot_cww_shd]
        zspec=self.zspec[mask1]
        zphot_doy=self.zphot_doy[mask1]
        zphot_shd=self.zphot_shd[mask1]
        zphot_cww_doy=self.zphot_cww_doy[mask1]
        zphot_cww_shd=self.zphot_cww_shd[mask1]
        
        plt.figure(1)    
        plt.clf()
        fig, axs = plt.subplots(2, 2)
        axs[0, 0].plot(zspec, zphot_doy,'.', label='Zphot_doy')
        axs[0, 0].plot(zspec, zspec, 'r--')
        axs[0, 0].set_title('DOY')
        axs[0, 0].legend()
        axs[0, 1].plot(zspec, zphot_shd,'.', label='Zphot_shd')
        axs[0, 1].plot(zspec, zspec, 'r--')
        axs[0, 1].set_title('SHD')
        axs[0, 1].legend()
        axs[1, 0].plot(zspec, zphot_cww_doy,'.', label='Zphot_cww_doy')
        axs[1, 0].plot(zspec, zspec, 'r--')
        axs[1, 0].set_title('CWW_DOY')
        axs[1, 0].legend()
        axs[1, 1].plot(zspec, zphot_cww_shd,'.', label='Zphot_cww_shd')
        axs[1, 1].plot(zspec, zspec, 'r--')
        axs[1, 1].set_title('CWW_SHD')
        axs[1, 1].legend()
        
        for ax in axs.flat:
            ax.set(xlabel='zspec', ylabel='zphot')
            ax.set(xlim=[0., 1], ylim=[0, 1])

        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()
        plt.savefig('Subplot_Zphot_phot_type.png')
        plt.clf()
        if without_bin:
            self.tag='not_bin'
            for i in lista:
                
                arr = (self.zspec-i)/(1+self.zspec)
                zbias.append(np.mean(self.zspec-i))
                outliers.append(self.outlier_rate(self.zspec-i))
                sigma68.append(self.sigma_68(arr))
                NMAD.append(np.mean(np.abs((self.zspec-i)/(1+self.zspec))))
                
            print (zbias)
        else:
            for i in range(len(lims)-1):
                
                if phot_type:
                    mask=(self.zspec>lims[i])&(self.zspec<lims[i+1])&(self.phot_type=='G')&(self.zspec>0)&(self.zspec!=-1)&(self.zspec!=1)
                    self.tag='Master_Cat_phot_type'
                    
                else:
                    mask = (self.zspec>lims[i])&(self.zspec<lims[i+1])
            
                zspec_tem = self.zspec[mask]
                bin_central.append((lims[i]+lims[i+1])/2.)
                
                    
                for i in lista:
                    i=i[mask]
                    arr = (zspec_tem-i)/(1+zspec_tem)
                    zbias.append(np.mean(zspec_tem-i))
                    outliers.append(self.outlier_rate(zspec_tem-i))
                    sigma68.append(self.sigma_68(arr))
                    NMAD.append(np.mean(np.abs((zspec_tem-i)/(1+zspec_tem))))
            
            print (sigma68)
            zbias_doy = []
            zbias_shd = []
            zbias_cww_doy = []
            zbias_cww_shd = []
            outliers_doy = []
            outliers_shd = []
            outliers_cww_doy = []
            outliers_cww_shd = []
            sigma68_doy = []
            sigma68_shd = []
            sigma68_cww_doy = []
            sigma68_cww_shd = []
            NMAD_doy = []
            NMAD_shd = []
            NMAD_cww_doy = []
            NMAD_cww_shd = []
        
            for i in range(0, len(zbias), 4):
                zbias_doy.append(zbias[i])
                zbias_shd.append(zbias[i+1])
                zbias_cww_doy.append(zbias[i+2])
                zbias_cww_shd.append(zbias[i+3])
                outliers_doy.append(outliers[i])
                outliers_shd.append(outliers[i+1])
                outliers_cww_doy.append(outliers[i+2])
                outliers_cww_shd.append(outliers[i+3])
                sigma68_doy.append(sigma68[i])
                sigma68_shd.append(sigma68[i+1])
                sigma68_cww_doy.append(sigma68[i+2])
                sigma68_cww_shd.append(sigma68[i+3])
                NMAD_doy.append(NMAD[i])
                NMAD_shd.append(NMAD[i+1])
                NMAD_cww_doy.append(NMAD[i+2])
                NMAD_cww_shd.append(NMAD[i+3])
            print (sigma68_doy)
        
        
        
            plt.figure(2)
            plt.clf()
            lista_zbias=[zbias_doy, zbias_shd, zbias_cww_doy, zbias_cww_shd]
            for i in lista_zbias:
                plt.plot(bin_central, i)
                plt.xlabel('zspec') 
                plt.ylabel('zphot')
                plt.title('Zbias measures')
                plt.legend(('Zphot_DOY', 'Zphot_SHD', 'Zphot_CWW_DOY', 'Zphot_CWW_SHD'))
                plt.savefig('Zbias_measures_'+self.tag+'.png')
        
            plt.figure(3)
            lista_outliers=[outliers_doy, outliers_shd, outliers_cww_doy, outliers_cww_shd]
            for i in lista_outliers:
                plt.plot(bin_central, i)
                plt.xlabel('zspec') 
                plt.ylabel('zphot')
                plt.title('Outliers measures')
                plt.legend(('Zphot_DOY', 'Zphot_SHD', 'Zphot_CWW_DOY', 'Zphot_CWW_SHD'))
                plt.savefig('Outliers_measures_'+self.tag+'.png')
            plt.figure(4)
            lista_sigma68=[sigma68_doy, sigma68_shd, sigma68_cww_doy, sigma68_cww_shd]
            for i in lista_sigma68:
                plt.plot(bin_central, i)
                plt.xlabel('zspec') 
                plt.ylabel('zphot')
                plt.title('Sigma68 measures')
                plt.legend(('Zphot_DOY', 'Zphot_SHD', 'Zphot_CWW_DOY', 'Zphot_CWW_SHD'))
                plt.savefig('Sigma68_measures_'+self.tag+'.png')
            plt.figure(5)
            lista_NMAD=[NMAD_doy, NMAD_shd, NMAD_cww_doy, NMAD_cww_shd]
            for i in lista_NMAD:
                plt.plot(bin_central, i)
                plt.xlabel('zspec') 
                plt.ylabel('zphot')
                plt.title('NMAD measures')
                plt.legend(('Zphot_DOY', 'Zphot_SHD', 'Zphot_CWW_DOY', 'Zphot_CWW_SHD'))
                plt.savefig('NMAD_measures_'+self.tag+'.png')
                
    def Histograma(self):
        
        import matplotlib.pyplot as plt
        import numpy as np
        import seaborn as sns
        
        #self.phot_type= self.datos['PHOT_TYPE_SHD']
        
        
        #mask1 = (self.phot_type=='G')&(self.zspec>0)
        '''
        zphot_cww_doy=self.zphot_cww_doy#[mask1]
        zphot_cww_shd=self.zphot_cww_shd#[mask1]
        zphot_doy=self.zphot_doy#[mask1]
        zphot_shd=self.zphot_shd#[mask1]
        '''
        lims= np.arange(0, 1, 0.02)
        
        
        plt.figure(1)
        plt.hist(x=self.zspec, bins=lims, bottom=True)
        plt.title('Distribución de Redshift Espectroscópico')
        plt.xlabel('ZSPEC')
        #plt.savefig('Distribucion_redshift.png')
        lim= np.arange(0, 3, 0.05)
        
        plt.figure(2)
        plt.hist(x=self.z_best_doy, bins=lim, label='DES only')#, alpha=0.7)
        plt.hist(x=self.z_best_shd, bins=lim, label='SHARKS-DES', alpha=0.7)
        plt.legend()
        plt.title('Comparación Zphot SHARKS-DES y solo DES')
        plt.xlabel('Z_BEST')
        #plt.savefig('Comparacion_zphot_doy_shd.png')
        #'''
        plt.figure(3)
        plt.figure(figsize=(10, 10), dpi=80)
        plt.plot(self.z_best_doy, self.z_best_shd, '.')
        plt.title('Comparación del zphot de SHARKS-DES frente solo a DES')
        plt.xlabel('Z_BEST_DOY')
        plt.ylabel('Z_BEST_SHD')
        #plt.savefig('Comparacion_zbest_doy_shd.png')
        #'''
        self.outliers_doy=self.datos["OUTLIERS_DOY"]
        self.outliers_shd=self.datos["OUTLIERS_SHD"]
        
        mask_doy = (self.outliers_doy==True)
        mask_shd = (self.outliers_shd==True)
        
        self.z_best_doy_mask1= self.z_best_doy[mask_doy]
        self.z_best_shd_mask1= self.z_best_shd[mask_doy]
        self.z_best_doy_mask2= self.z_best_doy[mask_shd]
        self.z_best_shd_mask2= self.z_best_shd[mask_shd]
        
        plt.figure(4)
        plt.plot(self.z_best_doy, self.z_best_shd, 'r.')
        plt.plot(self.z_best_doy_mask1, self.z_best_shd_mask1, 'b.')
        plt.plot(self.z_best_doy_mask2, self.z_best_shd_mask2, 'g.')
        plt.title('Comparación del zphot con diferencia de outliers por colores')
        plt.xlabel('Z_BEST_DOY')
        plt.ylabel('Z_BEST_SHD')
        #plt.savefig('Comparacion_zbest_doy_colores.png')
        
        
        self.mag_i_mask1=self.mag_i[mask_doy]
        self.mag_i_mask2=self.mag_i[mask_shd]
        
        self.mag_z_mask1=self.mag_z[mask_doy]
        self.mag_z_mask2=self.mag_z[mask_shd]
        
        self.mag_Ks_mask1=self.mag_Ks[mask_doy]
        self.mag_Ks_mask2=self.mag_Ks[mask_shd]
        limit= np.arange(20, 24, 0.1)
                
        plt.figure(5)       
        plt.hist(self.mag_i_mask1, bins=limit, density=True, label='DES only', color='blue')
        plt.hist(self.mag_i_mask2, bins=limit, density=True, label='SHARKS-DES', alpha=0.7, color='green')
        plt.xlabel('MAG_I')
        plt.legend()
        plt.savefig('MAG_I_doy_shd.png')
        
        plt.figure(6)
        plt.hist(self.mag_z_mask1, bins=limit, density=True, label='DES only', color='blue')
        plt.hist(self.mag_z_mask2, bins=limit, density=True, label='SHARKS-DES', alpha=0.7, color='green')
        plt.xlabel('MAG_Z')
        plt.legend()
        #plt.savefig('MAG_Z_doy_shd.png')
        plt.figure(7)
        plt.hist(self.mag_Ks_mask1, bins=limit, density=True, label='DES only', color='blue')
        plt.hist(self.mag_Ks_mask2, bins=limit, density=True, label='SHARKS-DES', alpha=0.7, color='green')
        plt.xlabel('MAG_Ks')
        plt.legend()
        #plt.savefig('MAG_KS_doy_shd.png')
        
        
        '''
        plt.figure(1)
        plt.hist(x=zphot_cww_doy, bins=lims, color='#F2AB6D')
        plt.title('CWW_DOY Template')
        plt.savefig('CWW_DOY_Template.png')
        plt.figure(2)
        plt.hist(x=zphot_cww_shd, bins=lims, color='#F2AB6D')
        plt.title('CWW_SHD Template')
        plt.savefig('CWW_SHD_Template.png')
        plt.figure(3)
        plt.hist(x=zphot_doy, bins=lims, color='#F2AB6D')
        plt.title('Pogiantti_DOY Template')
        plt.savefig('Pogiantti_DOY_Template.png')
        plt.figure(4)
        plt.hist(x=zphot_shd, bins=lims, color='#F2AB6D')
        plt.title('Pogiantti_SHD Template')
        plt.savefig('Pogiantti_SHD_Template.png')
        
        #plt.xlabel('ZPHOT')
        #plt.savefig(title+'.png')
        '''
    def magnitudes(self):
        import matplotlib.pyplot as plt
        import numpy as np
        from scipy.stats import kde
        
        #self.mag_ks=self.datos['MAG_Ks']
        plt.figure(1)
        nbins=200
        x=self.mag_i-self.mag_z
        y=self.mag_z-self.mag_Ks
        k = kde.gaussian_kde([x,y])
        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
 
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
        
        plt.plot( x, x-0.1)
        plt.xlim(-0.5, 1.2)
        plt.ylim(-1.75, 2.25)
        plt.title('Separación galaxias y estrellas')
        plt.ylabel('MAG_Z - MAG_Ks')
        plt.xlabel('MAG_I - MAG_Z')
        #plt.ylim(0, 40)
        #plt.xlim(0, 2)
        plt.show()
        plt.savefig('Separación_galaxias_y_estrellas.png')
        '''
        plt.figure(2)
        plt.plot( self.zphot_shd, self.mag_ks, '.')
        plt.title('Ks vs SHD')
        plt.ylabel('Ks')
        plt.xlabel('Zphot_SHD')
        plt.ylim(0, 40)
        plt.xlim(0, 2)
        plt.savefig('Ks_vs_SHD.png')
        plt.figure(3)
        plt.plot( self.zphot_cww_doy, self.mag_ks, '.')
        plt.title('Ks vs CWW_DOY')
        plt.ylabel('Ks')
        plt.xlabel('Zphot_CWW_DOY')
        plt.ylim(0, 40)
        plt.xlim(0, 2)
        plt.savefig('Ks_vs_CWW_DOY.png')
        plt.figure(4)
        plt.plot( self.zphot_cww_shd, self.mag_ks, '.')
        plt.title('Ks vs CWW_SHD')
        plt.ylabel('Ks')
        plt.xlabel('Zphot_CWW_SHD')
        plt.ylim(0, 40)
        plt.xlim(0, 2)
        plt.savefig('Ks_vs_CWW_SHD.png')
        '''
           
        
        
'''        
results_eazy = PhotoZ('eazy', ['gals_dessharks_1.txt', 'gals_dessharks_2.txt', 'gals_dessharks_3.txt', 'gals_dessharks_4.txt', 'gals_dessharks_5.txt'], 'sharks')

results_eazy.getStats()  
results_eazy.getStats_Bin([0.,0.1,0.2,0.3,0.4,0.5])    
#results_eazy.delta_z_1pz()        
'''        
        
results = PhotoZ("lephare_master_cat", 'master_cat_5s_bpz.fits', 'Master_Cat_bpz')

#results.get_stats_Master([0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7], False, True)
#results.plot()
results.Histograma()
#results.magnitudes()
