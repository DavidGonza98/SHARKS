# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 12:19:08 2021

@author: lpama
"""

from astropy.table import Table
import numpy as np

class Catalogo():
    
    def __init__(self, nombre, archivo):
        self.nombre= nombre
        self.archivo= archivo
        
        self.lineas=False
        self.Lleno=False
    
    def LeerArchivo(self):
        self.datos = Table.read(self.archivo, format= 'fits')

        self.Lleno=True
        
        if (self.Lleno):
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


    def estado(self, RA, DEC):
        self.linea=True
        if self.linea:
            print('La linea se ha cargado correctamente')
            
            print('El valor maximo de la columna',RA, 'es',np.max(self.datos[RA]),' y el minimo', np.min(self.datos[RA]))
            print('El valor maximo de la columna',DEC, 'es',np.max(self.datos[DEC]), ' y el minimo', np.min(self.datos[DEC]))
            
            
        else:
            print('No ha sido posible cargar las lineas')

    def mascara(self, Columna1, Columna2):
        
        if self.datos[Columna1]>0 and self.datos[Columna2]>12.3:
            
            self.mascara=self.datos
            return len(self.mascara)
        
            
           
sharks=Catalogo('Sharks', 'Sharks_sgp_e_2_cat_small.fits')

sharks.LeerArchivo()
ra = sharks.Extraer_columna('ALPHA_J2000')
sharks.estado('ALPHA_J2000', 'DELTA_J2000')
#print(sharks.mascara('MAGERR_AUTO','MAG_AUTO'))
#Probando




        
        
        
        
        
        
        