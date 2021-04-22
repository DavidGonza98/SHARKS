# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 12:19:08 2021

@author: lpama
"""

from astropy.table import Table
import numpy as np
from astropy.io import ascii
import smatch
import math as mt

class Catalogo():
    
    def __init__(self, nombre, archivo, RA, DEC):
        self.nombre= nombre
        self.archivo= archivo
        
        #self.lineas=False
        self.Lleno=False
        self.RA=RA
        self.DEC=DEC
        self.nside=4096 
        self.maxmatch=1 
        self.radius= 1/3600.
        self.match=False
        self.mask=False
        self.healpix = False
    
    def LeerArchivo(self):
        self.datos = Table.read(self.archivo, format= 'fits')

        self.Lleno=True
        
        if self.Lleno:
            self.area= np.shape(self.datos)
            self.longitud= len(self.datos)
            print('El catálogo ',self.nombre,' se ha cargado correctamente')
            print('El área del catálogo es ', self.area, 'y el número de objetos será ', self.longitud )
            
        else:
            print('No ha sido posible cargar los datos del catálogo ', self.nombre)
        
    def Extraer_columna(self, DameColumna):
        
        if DameColumna in self.datos.keys():
            print('La columna', DameColumna, 'tiene una longitud de', len(self.datos[DameColumna]), 'un maximo de', np.max(self.datos[DameColumna]),'y un minimo de', np.min(self.datos[DameColumna]))
                         
            return self.datos[DameColumna]

        else:
            print ('El nombre de la columna ingresada no se encuentra en el catalogo', self.nombre)


    def Match(self, ObjectCatalog):
        self.match=True
        if self.Lleno and self.nombre!=ObjectCatalog.nombre:
            self.ra1_matched=self.datos[self.RA]
            self.dec1_matched=self.datos[self.DEC]
            
            self.datos2 = ObjectCatalog.datos
            self.RA2 = ObjectCatalog.RA
            self.DEC2 = ObjectCatalog.DEC
            self.ra2_matched=self.datos2[self.RA2]
            self.dec2_matched=self.datos2[self.DEC2]
            #print('La longitud de', self.RA, 'es', len(ra1_matched), 'y la de', self.DEC, 'es', len(dec1_matched) )
            #print('La longitud de', ObjectCatalog[2], 'es', len(ra2_matched), 'y la de', ObjectCatalog[3], 'es', len(dec2_matched), 'donde estos datos pertenecen al catalogo', ObjectCatalog[0])
            self.matches = smatch.match(self.ra1_matched, self.dec1_matched, self.radius, self.ra2_matched, self.dec2_matched, nside=self.nside, maxmatch=self.maxmatch)
            
            self.assoc1 = self.datos[self.matches['i1']]
            self.assoc2 = self.datos2[self.matches['i2']]


            print('La longitud del catalogo ya habiendo realizado el matched es', len(self.assoc1))
            self.match=True
            self.nombreMatch = ObjectCatalog.nombre 
            self.MainCatalog()
            
        else:
            print('No se ha podido hacer el match')
            
    
       
    
    def AreaEspacial(self):
        
        lado_menor1=np.max(self.ra1_matched) - np.min(self.ra1_matched)
        lado_mayor1=mt.sin(np.max(self.dec1_matched)) - mt.sin(np.min(self.dec1_matched))
        self.area1= (180/mt.pi)*lado_menor1*lado_mayor1 
        
        lado_menor2=np.max(self.ra2_matched) - np.min(self.ra2_matched)
        lado_mayor2=mt.sin(np.max(self.dec2_matched)) - mt.sin(np.min(self.dec2_matched))
        self.area2= (180/mt.pi)*lado_menor2*lado_mayor2 
        
        print('El valor del area de la esfera del catalogo '+self.nombre+' es',self.area1,'grados cuadrados y el area del catalogo '+self.nombreMatch+' es',self.area2, 'grados cuadrados')

    def MainCatalog (self):
        nombres=[]
        arrays=[]
        
        if self.match:
            
            for name in self.assoc1.columns:
                nombres.append(name)
                arrays.append(self.assoc1[name])
            
            for name in self.assoc2.columns:
                arrays.append(self.assoc2[name])
                if name in nombres:
                    name= name+'_'+self.nombreMatch
                nombres.append(name)
            
            self.datos= Table(data=arrays, names=nombres)

            if self.RA2 in self.assoc1.columns:
                self.RA2 += '_'+self.nombreMatch
            if self.DEC2 in self.assoc1.columns:
                self.DEC2 += '_'+self.nombreMatch
            #print(self.datos.columns)
                          
        else:
            print('No se ha podido unificar')
        
        #return self.unify
    '''
    def DefineMain(self, datos):
        if datos == 'sharks':
            self.unify = self.datos
        else:
            self.unify = self.datos2
    '''
    
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

            self.mask=True
            
            self.datos = self.datos[mask]

        else:
            print('File not read')
        
   
    def saveSample(self, listasColumnas):
        nombres=[]
        arrays=[]
        for name in listasColumnas:
            nombres.append(name)
            arrays.append(self.datos[name])
       
        self.t= Table(data=arrays, names=nombres)
   
    def createSample(self,format='fits', nameSample='mysample'):
        #format puede ser csv, ascii o fits
        
            
        if format=='fits':
            self.t.write(nameSample+'.fits')

        else:            
            ascii.write(self.t, format=format, output=nameSample+'.'+format, overwrite=True)
        
    def giveFluxes(self, MAG, tipo):
        
        if tipo =='cgs':
            
            flux= 10**((-48.6 - (self.Extraer_columna(MAG[0])))/2.5)

            fluxerr=self.Extraer_columna(MAG[1])*flux/1.086
        
            print('El valor del flujo es', flux,'y del error',fluxerr)
            self.datos['FLUX_'+MAG[0]] = flux
        
            self.datos['FLUXERR_'+MAG[0]] = fluxerr 
            
        elif tipo== 'jansky':
            
            self.flux = 10**((self.Extraer_columna(MAG[0])+8.9)/-2.5)
            
            fluxerr=self.Extraer_columna(MAG[1])*flux/1.086
            
            print('El valor del flujo es', flux,'y del error',fluxerr)
            self.datos['FLUX_'+MAG[0]] = flux
        
            self.datos['FLUXERR_'+MAG[0]] = fluxerr
        
        else:
            print('El tipo de flujo introducido es incorrecto')
            
    def estado(self):
        #self.linea=True
        self.match=True
        self.Lleno=True
        self.mask=True
        if self.Lleno:
            
            print('El catalogo '+self.nombre+' se ha cargado correctamente')
            
            if  self.match:
                print('Se ha realizado correctamente el match entre los catalogos '+self.nombre+' y '+self.nombreMatch)
                
            else:
                print ('No se ha podido realizar el match entre los catalogos '+self.nombre+' y '+self.nombreMatch)
            
            if self.mask:
                print('Tras haber realizado correctamente la mascara tenemos ',len(self.datos),' objetos, siendo el area del catalogo ',self.area)
                
            else:
                print('No se ha podido realizar la mascara correctamente')
            
        else:
            print('No ha sido posible leer el catalogo '+self.nombre+' correctamente')

        
        
           
