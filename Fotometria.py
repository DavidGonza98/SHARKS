# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 20:52:44 2021

@author: lpama
"""
from astropy.table import Table
import numpy as np
import smatch
import math as mt
import seaborn as sns
import pylab as plt
from scipy import optimize, stats
import pandas as pd
     


class fotometria():
    """
    Class to define the fotometry of the catalog.
    """
    
    def __init__(self,ObjectCatalog):
        """
        Instantiate the fotometry of the catalog.

        Parameters
        ----------
        ObjectCatalog : `str`
            Name of the catalog we have defined before on catalogo_sharks.

        Returns
        -------
        None.

        """
        
        self.datos = ObjectCatalog.datos
                             
        
   
    def ajuste_lineal(self, magnitud1, magnitud2, ObjectCatalog, ObjectCatalog1, quitarOutliers=False):
        """
        Method that performs a linear adjustment of the parameters magnitud1 and magnitud2

        Parameters
        ----------
        magnitud1 : `str`
            X magnitude for the linear adjustment.
        magnitud2 : `str`
            Y magnitude for the linear adjustment.
        ObjectCatalog : `str`
            Catalog from which magnitud1 comes from.
        ObjectCatalog1 : `str`
            DESCRIPTION.
        quitarOutliers : `bool`, optional
            Catalog from which magnitud1 comes from.. The default is False.

        Returns
        -------
        None.

        """
        
        plt.figure(2)      
        self.p=plt.polyfit(self.datos[magnitud1], self.datos[magnitud2], 1)
        Y=self.p[0]*self.datos[magnitud1] + self.p[1] 
        plt.figure(1)
        plt.plot(self.datos[magnitud1], Y, label = "Ajuste y = %.3f + %.3f x" % (self.p[1],self.p[0]))
        plt.plot(self.datos[magnitud1], self.datos[magnitud2], '.', label='Datos')
        plt.legend()
        plt.title('Ajuste sin errores y con outliers')
        plt.xlabel(magnitud1)
        plt.ylabel(magnitud2)
        title='Sinerrores_y_conoutliers_de_archivoconjutode_' + ObjectCatalog.nombre+ '_y_'  + ObjectCatalog1.nombre
        #plt.savefig(title+'.png')
        print ('El valor de la pendiente del ajuste lineal es', self.p[0], 'y su ordenada en el origen', self.p[1]  )
        
        
            
    def errorfunc(self, magnitud1, magnitud2, errormagnitud1, ObjectCatalog, ObjectCatalog1, quitarOutliers=False, pendienteuno=False):
        """
        Method that calculates the error function. Here we have needed from the library optimize the function leastsq to obtain m and b from the linear adjustmet including errors y=m*x+b.

        Parameters
        ----------
        magnitud1 : `str`
            X magnitude for the linear adjustment.
        magnitud2 : `str`
            Y magnitude for the linear adjustment.
        errormagnitud1 : `str`
            Error of magnitud1.
        ObjectCatalog : `str`
            Catalog from which magnitud1 comes from.
        ObjectCatalog1 : `str`
            Catalog from which magnitud2 comes from.
        quitarOutliers : `bool`, optional
            Here we decide if we want out of our linear adjustment the outliers. The default is False.
        pendienteuno : `bool`, optional
            Here we decide if we want to force the slope of the linear adjustment to be equal to 1. The default is False.

        Returns
        -------
        None.

        """
        
        fitfunc = lambda p, x: p[0] + p[1] * x
        errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
        
        if self.p.any(None):
            print ('Todo esta bien')

        else:
            self.ajuste_lineal(magnitud1, magnitud2)

        #pinit = [self.p[1],self.p[0]]
        pinit = [1,1]
        out = optimize.leastsq(errfunc, pinit, args=(self.datos[magnitud1], self.datos[magnitud2],  self.datos[errormagnitud1]), full_output=1)
        #print(pinit)
        self.pfinal = out[0]
        print(self.pfinal)
        ajuste_con_error2= self.pfinal[1]*self.datos[magnitud1] + self.pfinal[0]
        plt.figure(3)
        plt.plot(self.datos[magnitud1], ajuste_con_error2, label = "Ajuste y = %.3f + %.3f x" % (self.pfinal[0],self.pfinal[1]))
        plt.plot(self.datos[magnitud1], self.datos[magnitud2], '.', label='Datos')
        plt.xlabel(magnitud1)
        plt.ylabel(magnitud2)
        plt.legend()
        plt.title('Ajuste con error y con outliers')
        title1='Conerrores_y_conoutliers_de_archivoconjutode_' + ObjectCatalog.nombre+ '_y_'  + ObjectCatalog1.nombre
        plt.savefig(title1+'.png')
        if quitarOutliers:
            z = (self.datos[magnitud2] - ajuste_con_error2)/self.datos[errormagnitud1]

            pf = pd.DataFrame(zip(self.datos[magnitud2], self.datos[magnitud1], self.datos[errormagnitud1]))
            pf_sin_outliers = pf[(np.abs(z) < 2.5)]
            plt.figure(4)
            m, b = plt.polyfit(pf_sin_outliers[1], pf_sin_outliers[0], 1)
            plt.plot(pf_sin_outliers[1], m*pf_sin_outliers[1]+b, 'r-', label = "Ajuste y = %.3f + %.3f x" % (b,m))
            plt.plot(pf_sin_outliers[1], pf_sin_outliers[0], 'b.', label='Datos')
            plt.legend()
            plt.xlabel(magnitud1)
            plt.ylabel(magnitud2)
            plt.title('Ajuste sin errores y sin outliers')
            title2='Sinerrores_y_sinoutliers_de_archivoconjutode_' + ObjectCatalog.nombre+ '_y_'  + ObjectCatalog1.nombre
            plt.savefig(title2+'.png')
            
            pinit = [1,1]
            out_sin_outliers = optimize.leastsq(errfunc, pinit, args=(pf_sin_outliers[1], pf_sin_outliers[0],  pf_sin_outliers[2]), full_output=1)
            self.pfinal_sin_outliers= out_sin_outliers[0]
            ajuste_sin_outliers= self.pfinal_sin_outliers[1]*pf_sin_outliers[1] + self.pfinal_sin_outliers[0]
            
            plt.figure(5)
            plt.plot(pf_sin_outliers[1], ajuste_sin_outliers, 'm', label = "Ajuste y = %.3f + %.3f x" % (self.pfinal_sin_outliers[0],self.pfinal_sin_outliers[1]))
            plt.plot(pf_sin_outliers[1], pf_sin_outliers[0], 'y.', label='Datos')
            plt.legend()
            plt.xlabel(magnitud1)
            plt.ylabel(magnitud2)
            plt.title('Ajuste con error y sin outliers')
            title3='Conerrores_y_sinoutliers_de_archivoconjutode_' + ObjectCatalog.nombre+ '_y_'  + ObjectCatalog1.nombre
            plt.savefig(title3+'.png')
        
        if pendienteuno:
            fitfunc = lambda p, x: p[0] + x
            pinit = [1,1]
            out = optimize.leastsq(errfunc, pinit, args=(self.datos[magnitud1], self.datos[magnitud2],  self.datos[errormagnitud1]), full_output=1)
            self.pfinal = out[0]
            ajuste_con_error2= self.datos[magnitud1] + self.pfinal[0]
            z = (self.datos[magnitud2] - ajuste_con_error2)/self.datos[errormagnitud1]

            pf = pd.DataFrame(zip(self.datos[magnitud2], self.datos[magnitud1], self.datos[errormagnitud1]))
            pf_sin_outliers = pf[(np.abs(z) < 2.5)]
            
            m, b = plt.polyfit(pf_sin_outliers[1], pf_sin_outliers[0], 1)
            plt.figure(6)
            plt.plot(pf_sin_outliers[1], pf_sin_outliers[1]+b, 'r-', label='Ajuste para m=1')
                
            
            pinit = [1,1]
            out_sin_outliers = optimize.leastsq(errfunc, pinit, args=(pf_sin_outliers[1], pf_sin_outliers[0],  pf_sin_outliers[2]), full_output=1)
            self.pfinal_sin_outliers= out_sin_outliers[0]
                
            ajuste_sin_outliers= pf_sin_outliers[1] + self.pfinal_sin_outliers[0]
            
                
            plt.plot(pf_sin_outliers[1], ajuste_sin_outliers, 'm', label = "Ajuste y = %.3f + x" % (self.pfinal_sin_outliers[0]) )
            plt.plot(pf_sin_outliers[1], pf_sin_outliers[0], 'y.', label='Datos')
            plt.legend()
            plt.title('Ajuste con pendiente m=1, con errores y sin outliers')
            plt.xlabel(magnitud1)
            plt.ylabel(magnitud2)
            title4='Conerrores_y_sinoutliers_de_archivoconjutode_' + ObjectCatalog.nombre+ '_y_'  + ObjectCatalog1.nombre
            plt.savefig(title4+'.png')
        

      
