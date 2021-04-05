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

        self.Lleno=False
        
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
            


    def Match(self, ra1, dec1, ra2, dec2):
        if (self.Lleno):
            ra1_matched=self.datos1[ra1]
            dec1_matched=self.datos1[dec1]
            ra2_matched=self.datos2[ra2]
            dec2_matched=self.datos2[dec2]
            print('La longitud de', ra1, 'es', len(ra1_matched), 'y la de', dec1, 'es', len(dec1_matched) )
            print('La longitud de', ra2, 'es', len(ra2_matched), 'y la de', dec2, 'es', len(dec2_matched) )
            self.matches = smatch.match(ra1_matched, dec1_matched, self.radius, ra2_matched, dec2_matched, nside=self.nside, maxmatch=self.maxmatch)
            
            self.assoc1 = self.datos1[self.matches['i1']]
            self.assoc2 = self.datos2[self.matches['i2']]


#            self.ra1_matched= self.ra1 [self.matches['i1']]
#            self.dec1_matched= self.dec1 [self.matches['i1']]
#            self.ra2_matched= self.ra2 [self.matches['i2']]
#            self.dec2_matched= self.dec2 [self.matches['i2']]
            print('La longitud del catalogo ya habiendo realizado el matched es', len(self.assoc1))
            
            
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
            
           
        
    def MainCatalog (self):
        nombres=[]
        arrays=[]
        
        for name in self.assoc1.columns:
            nombres.append(name)
            arrays.append(self.assoc1[name])
        
        for name in self.assoc2.columns:
            nombres.append(name)
            arrays.append(self.assoc2[name])
        
        self.unify= Table(data=arrays, names=nombres)    
        return self.unify

    def DefineMain(self, datos):
        if datos == 'sharks':
            self.MainCatalog = self.datos1
        else:
            self.MainCatalog = self.datos2

    def Extraer_columna(self, DameColumna, matched):
        
        if matched == 'Si':
            if DameColumna in self.unify.keys():
                print('La columna', DameColumna, 'tiene una longitud de', len(self.unify[DameColumna]), 'un maximo de', np.max(self.unify[DameColumna]),'y un minimo de', np.min(self.unify[DameColumna]))
                return self.unify[DameColumna]
        
        elif matched == 'No':
            if DameColumna in self.datos1.keys():
                print('La columna', DameColumna, 'tiene una longitud de', len(self.datos1[DameColumna]), 'un maximo de', np.max(self.datos1[DameColumna]),'y un minimo de', np.min(self.datos1[DameColumna]))
                return self.datos1[DameColumna]
        
            if DameColumna in self.datos2.keys():
                print('La columna', DameColumna, 'tiene una longitud de', len(self.datos2[DameColumna]), 'un maximo de', np.max(self.datos2[DameColumna]),'y un minimo de', np.min(self.datos2[DameColumna]))
                return self.datos2[DameColumna]

        else:
            print ('El nombre de la columna ingresada no se encuentra ningun catalogo')
            



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
            return np.greater(self.Extraer_columna(condicion[0], 'Si'), condicion[2])
        elif condicion[1] == 'greater_equal':
            return np.greater_equal(self.Extraer_columna(condicion[0], 'Si'), condicion[2])
        elif condicion[1] == 'less':
            return np.less(self.Extraer_columna(condicion[0], 'Si'), condicion[2])
        elif condicion[1] == 'less_equal':
            return np.less_equal(self.Extraer_columna(condicion[0], 'Si'), condicion[2])
        elif condicion[1] == 'equal':
            return np.equal(self.Extraer_columna(condicion[0], 'Si'), condicion[2])
        elif condicion[1] == 'not_equal':
            return np.not_equal(self.Extraer_columna(condicion[0], 'Si'), condicion[2])
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

        
            self.unify_mask = self.unify[mask]
            
            print('La longitud del catalogo habiendo hecho el matched y la mascara sera: ',len(self.unify_mask))
        else:
            print('File not read')
    
    def ajuste_lineal(self, x, y):
        
        #print(self.datos1[x])
        
        p=plt.polyfit(self.unify_mask[x], self.unify_mask[x], 1)
        Y=p[0]*self.unify_mask[x] + p[1] 
        plt.plot(self.unify_mask[x], Y)
        plt.plot(self.unify_mask[x], self.unify_mask[x], '.')
        
        
            
            
            
        
sharks=astrometria('Sharks','Sharks_sgp_e_2_cat_small.fits', '2mass', '2mass.fit')   

sharks.LeerArchivo()

sharks.Match('ALPHA_J2000', 'DELTA_J2000', 'RAJ2000', 'DEJ2000')
sharks.MainCatalog()

sharks.Extraer_columna('Kmag', 'Si')
      

#sharks.DefineMain('Sharks') #o podrias poner '2mass'



'''
    #self.MainCatalog = union assoc1 y assoc2
1) crear un dictionario nuevo vacio.
2) Leer las columnas de sharks y añadirlas al nuevo
3) leer las columnas de 2mass y añadirlas al nuevo

'''


#sharks.Matches('MAG_AUTO', 'MAGERR_AUTO', 'Kmag', 'e_Kmag')       

#sharks.Distancia_Angular()    




#plot = sharks.Histograma('Distancia angular')   
#plt.savefig('dist_angular.png')

sharks.mascara([['MAG_AUTO','greater_equal',12.3],['MAGERR_AUTO','greater',0.]])

sharks.ajuste_lineal('MAG_AUTO', 'Kmag')  
        
        
        
        
        
        
        
        
        
        
        
        

         
    

















