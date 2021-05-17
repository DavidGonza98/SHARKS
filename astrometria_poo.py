"""
Created on Wed Mar 31 16:39:46 2021

@author: lpama
"""

#from astropy.table import Table
#import numpy as np
#import smatch
#import math as mt
#import seaborn as sns
#import matplotlib.pyplot as plt
#from scipy import optimize, stats
#import pandas as pd

class astrometria():
    """
    Class to define the astrometry of the catalog.
    """
    
    def __init__(self, ObjectCatalog):
        """
        Instantiate the astrometry of the catalog.

        Parameters
        ----------
        ObjectCatalog : `str`
            Name of the catalog we have defined before on catalogo_sharks.

        Returns
        -------
        None.

        """
        
        self.datos = ObjectCatalog.datos
        self.ra1 = ObjectCatalog.RA
        self.ra2 = ObjectCatalog.RA2
        self.dec1 = ObjectCatalog.DEC
        self.dec2 = ObjectCatalog.DEC2
        self.nombre = ObjectCatalog.nombre
        self.nombreMatch = ObjectCatalog.nombreMatch
        ObjectCatalog.estado()
        self.DefineVariables()
        
           
    def DefineVariables(self):
        """
        Define variables to differentiate RA and DEC from each catalog

        Returns
        -------
        None.

        """
        
        self.ra2_matched=self.datos[self.ra2]
        self.dec2_matched=self.datos[self.dec2]
        self.ra1_matched=self.datos[self.ra1]
        self.dec1_matched=self.datos[self.dec1]
    
       
    
    def Distancia_Angular(self):
        """
        Angular distance

        Returns
        -------
        `np.ndarray`
            Angular distance.

        """        
        import numpy as np
        import math as mt
        self.cosgamma = np.zeros(len(self.datos))

        
        for i in range(len(self.datos)):
            self.gamma = mt.cos(90-self.dec1_matched[i])*mt.cos(90-self.dec2_matched[i])+mt.sin(90-self.dec1_matched[i])*mt.sin(90-self.dec2_matched[i])*mt.cos(self.ra1_matched[i]-self.ra2_matched[i])
            self.cosgamma[i] = (mt.acos(self.gamma)*3600)
        return self.cosgamma
    
    def Histograma(self, tipo):
        """
        Representate a histogram of tipo.

        Parameters
        ----------
        tipo : `str`
            String describing what histogram we want.

        Returns
        -------
        None.

        """
        import matplotlib.pyplot as plt
        import seaborn as sns
    
        plt.figure(1)
        if tipo == 'Distancia angular':
            
            sns.displot(self.cosgamma, kde=True)
            plt.title('Histograma de distancia angular')
            title='Histograma_Distancia_Angular' + self.nombre+ '_y_'  + self.nombreMatch
            plt.xlabel('Segundos de Arco')
            plt.savefig(title+'.png')
            
        
        elif tipo == 'Jointplot':
            sns.jointplot(x=self.KMAG, y=self.MAG_AUTO, kind="reg", truncate=False, xlim=(5, 25), ylim=(5, 25), color="m", height=7)
        
    def plot(self, tipo):
        """
        Method plotting the variable that is specified in tipo.

        Parameters
        ----------
        tipo : `str`
            String defining of what we want to plot.

        Returns
        -------
        None.

        """
        import matplotlib.pyplot as plt
        import numpy as np
        plt.clf()
        dra=self.ra1_matched-self.ra2_matched
        ddec=self.dec1_matched-self.dec2_matched
        #plt.plot(dra*3600, ddec*3600, '.')
        plt.plot(0,0, 'r+')
        
        if tipo=='X' or tipo=='Y':
            sc=plt.scatter(dra*3600, ddec*3600, c= self.datos[tipo]- np.max(self.datos[tipo])/2)
            cbar=plt.colorbar(sc)
            cbar.set_label(tipo, rotation=270)
            plt.xlabel('dRA [arcsec]')
            plt.ylabel('dDEC [arcsec]')
            title= self.nombre+ '-'  + self.nombreMatch+' pos offsets of '+ tipo
            plt.title(title)
            plt.savefig(title+'.png')     
        else:
            sc=plt.scatter(dra*3600, ddec*3600, c= self.datos[tipo])
            cbar=plt.colorbar(sc)
            cbar.set_label(tipo, rotation=270)
            plt.xlabel('dRA [arcsec]')
            plt.ylabel('dDEC [arcsec]')
            title= self.nombre+ '-'  + self.nombreMatch+' pos offsets of ' +tipo
            plt.title(title)
            plt.savefig(title+'.png')    
        
        
        
        
        
        
        

         
    

















