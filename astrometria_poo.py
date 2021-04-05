# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 16:39:46 2021

@author: lpama
"""

from astropy.table import Table
import numpy as np
import smatch
import math as mt
import seaborn as sns
import pylab as plt

class astrometria():
    
    def __init__(self, nombre1, archivo1, nombre2, archivo2):
        self.nombre1=nombre1
        self.archivo1=archivo1
        self.nombre2=nombre2
        self.archivo2=archivo2
        self.nside=4096 
        self.maxmatch=1 
        self.radius= 1/3600.
        
        
    def LeerArchivo(self):
        self.datos1 = Table.read(self.archivo1, format= 'fits')
        self.datos2 = Table.read(self.archivo2, format= 'fits')

        self.Lleno=True
        
        if (self.Lleno):
            self.area1= np.shape(self.datos1)
            self.longitud1= len(self.datos1)
            self.area2= np.shape(self.datos2)
            self.longitud2= len(self.datos2)
            print('El catálogo ',self.nombre1,' se ha cargado correctamente')
            print('El área del catálogo es ', self.area1, 'y el número de objetos será ', self.longitud1 )
            print('El catálogo ',self.nombre2,' se ha cargado correctamente')
            print('El área del catálogo es ', self.area2, 'y el número de objetos será ', self.longitud2 )
            
    def Extraer_columna(self, DameColumna):
        
                
        if DameColumna in self.datos1.keys():
            print('La columna', DameColumna, 'tiene una longitud de', len(self.datos1[DameColumna]), 'un maximo de', np.max(self.datos1[DameColumna]),'y un minimo de', np.min(self.datos1[DameColumna]))
            return self.datos1[DameColumna]
        
        if DameColumna in self.datos2.keys():
            print('La columna', DameColumna, 'tiene una longitud de', len(self.datos2[DameColumna]), 'un maximo de', np.max(self.datos2[DameColumna]),'y un minimo de', np.min(self.datos2[DameColumna]))
            return self.datos2[DameColumna]

        else:
            print ('El nombre de la columna ingresada no se encuentra ningun catalogo')
            
    def Match(self, ra1, dec1, ra2, dec2):
        if (self.Lleno):
            self.ra1=self.datos1[ra1]
            self.dec1=self.datos1[dec1]
            self.ra2=self.datos2[ra2]
            self.dec2=self.datos2[dec2]
            print('La longitud de', ra1, 'es', len(self.ra1), 'y la de', dec1, 'es', len(self.dec1) )
            print('La longitud de', ra2, 'es', len(self.ra2), 'y la de', dec2, 'es', len(self.dec2) )
            self.matches = smatch.match(self.ra1, self.dec1, self.radius, self.ra2, self.dec2, nside=self.nside, maxmatch=self.maxmatch)
            
            self.ra1_matched= self.ra1 [self.matches['i1']]
            self.dec1_matched= self.dec1 [self.matches['i1']]
            self.ra2_matched= self.ra2 [self.matches['i2']]
            self.dec2_matched= self.dec2 [self.matches['i2']]
            print('La longitud de', ra1,',', dec1,',', ra2,',', dec2, 'ya habiendo realizado el matched es', len(self.ra1_matched))
            
            
        else:
            print('No se ha podido hacer el match')
            
            
    def Matches(self, MAG_AUTO, MAGERR_AUTO, KMAG, E_KMAG):
        
        if self.Lleno:
            self.MAG_AUTO= self.datos1[MAG_AUTO] [self.matches['i1']]
            self.MAGERR_AUTO= self.datos1[MAGERR_AUTO] [self.matches['i1']]
            self.KMAG= self.datos2[KMAG] [self.matches['i2']]
            self.E_KMAG= self.datos2[E_KMAG] [self.matches['i2']]
            
            print('La longitud de', MAG_AUTO, 'es', len(self.datos1[MAG_AUTO]), 'y la de', MAGERR_AUTO, 'es', len(self.datos1[MAGERR_AUTO]) )
            print('La longitud de', KMAG, 'es', len(self.datos2[KMAG]), 'y la de', E_KMAG, 'es', len(self.datos2[E_KMAG]) )
            
            print('La longitud de', MAG_AUTO,',', MAGERR_AUTO,',', KMAG,',', E_KMAG, 'ya habiendo realizado el matched es', len(self.E_KMAG))
            
           
        
            
    def Distancia_Angular(self):
        self.cosgamma = np.zeros(len(self.matches))
        
        for i in range(len(self.matches)):
            self.gamma = mt.cos(90-self.dec1_matched[i])*mt.cos(90-self.dec2_matched[i])+mt.sin(90-self.dec1_matched[i])*mt.sin(90-self.dec2_matched[i])*mt.cos(self.ra1_matched[i]-self.ra2_matched[i])
            self.cosgamma[i] = (mt.acos(self.gamma)*3600)
            return self.cosgamma
    
    def Histograma(self, tipo):
        
        
        if tipo == 'Distancia angular':
            return sns.displot(self.cosgamma, kde=True)
        
        elif tipo == 'Jointplot':
            return sns.jointplot(x=self.KMAG, y=self.MAG_AUTO, kind="reg", truncate=False, xlim=(5, 25), ylim=(5, 25), color="m", height=7)
        
            
    def __makeCondition(self, condicion):
        
        if condicion[1] == 'greater':
            return np.greater(self.Extraer_columna(condicion[0]), condicion[2])
        elif condicion[1] == 'greater_equal':
            return np.greater_equal(self.Extraer_columna(condicion[0]), condicion[2])
        elif condicion[1] == 'less':
            return np.less(self.Extraer_columna(condicion[0]), condicion[2])
        elif condicion[1] == 'less_equal':
            return np.less_equal(self.Extraer_columna(condicion[0]), condicion[2])
        elif condicion[1] == 'equal':
            return np.equal(self.Extraer_columna(condicion[0]), condicion[2])
        elif condicion[1] == 'not_equal':
            return np.not_equal(self.Extraer_columna(condicion[0]), condicion[2])
        else:
            return 0


    def mascara(self, listaCondiciones):
        '''condiciones sera una lista de listas, algo como: [['MAG_AUTO','greater',12.3], ['MAGERR_AUTO','greater',0.]]
        Los valores aceptados de condicional son: greater, greater_equal, less, less_equal, equal, not_equal
        '''
        if self.Lleno:
            mask = (self.__makeCondition(listaCondiciones[0]))
            if len(listaCondiciones)>1:
                for m in range(len(listaCondiciones)-1):
                    mask = mask*(self.__makeCondition(listaCondiciones[m+1]))

        
            self.datos1_mask = self.datos1[mask]
            self.datos2_mask = self.datos2[mask]

        else:
            print('File not read')
    
    def ajuste_lineal(self, x, y):
        
        #print(self.datos1[x])
        
        return plt.polyfit(self.datos1_mask[x], self.datos2_mask[y], 1)
        
            
            
            
            
        
sharks=astrometria('Sharks','Sharks_sgp_e_2_cat_small.fits', '2mass', '2mass.fit')   

sharks.LeerArchivo()

#sharks.Extraer_columna('MAG_AUTO')
sharks.Extraer_columna('Kmag')

sharks.Match('ALPHA_J2000', 'DELTA_J2000', 'RAJ2000', 'DEJ2000')
       
sharks.Matches('MAG_AUTO', 'MAGERR_AUTO', 'Kmag', 'e_Kmag')       

sharks.Distancia_Angular()    




sharks.Histograma('Distancia angular')   

#sharks.mascara([['MAG_AUTO','greater_equal',12.3],['MAGERR_AUTO','greater',0.]])

#sharks.ajuste_lineal('MAG_AUTO', 'Kmag')  
        
        
        
        
        
        
        
        
        
        
        
        

         
    

















