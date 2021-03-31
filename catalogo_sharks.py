# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 12:19:08 2021

@author: lpama
"""

from astropy.table import Table
import numpy as np

class Catalogo():
    
    def __init__(self, nombre, archivo, RA, DEC):
        self.nombre= nombre
        self.archivo= archivo
        
        self.lineas=False
        self.Lleno=False
        self.RA=RA
        self.DEC=DEC
    
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



    def estado(self):
        self.linea=True
        if self.linea:
            print('La linea se ha cargado correctamente')
            
            self.Extraer_columna(self.RA)
            self.Extraer_columna(self.DEC)
           
            
            
        else:
            print('No ha sido posible cargar las lineas')

    #def mascara(self, Columna1, Columna2):

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

        
            self.datos = self.datos[mask]

        else:
            print('File not read')
        
            
           
sharks=Catalogo('Sharks', 'Sharks_sgp_e_2_cat_small.fits', 'ALPHA_J2000', 'DELTA_J2000')


sharks.LeerArchivo()
ra = sharks.Extraer_columna('ALPHA_J2000')
sharks.estado()


print(len(sharks.datos))

sharks.mascara([['MAG_AUTO','greater_equal',12.3],['MAGERR_AUTO','greater',0.]])
print('Despues de mascarar', len(sharks.datos))




        
        
        
        
        
        
        
