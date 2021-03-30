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
        self.datos = Table.read(self.archivo, format= 'fits')
        
        self.Lleno=False
    
    def LeerArchivo(self):

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
            print('La columna', DameColumna, 'tiene una longitud de', len(self.datos[DameColumna]))
            return self.datos[DameColumna]

        else:
            print ('El nombre de la columna ingresada no se encuentra en el catalogo', self.nombre)

    def mascara(self, Columna1, Columna2):
        
        if self.datos[Columna1]>0 and self.datos[Columna2]>12.3:
            
            self.mascara=self.datos
            return len(self.mascara)
        
            
           
sharks=Catalogo('Sharks', 'sharks_sgpe.fits')


print(sharks.Extraer_columna('RA'))
print(sharks.mascara('RA','MAG_AUTO'))




        
        
        
        
        
        
        