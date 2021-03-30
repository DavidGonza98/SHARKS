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
        
        self.Lleno=False
    
    def LeerArchivo(self):
        self.datos= Table.read(self.archivo)
        self.Lleno=True
        
        if (self.Lleno):
            self.area= np.shape(self.datos)
            self.longitud= len(self.datos)
            print('El catálogo ',self.nombre,' se ha cargado correctamente')
            print('El área del catálogo es ', self.area, 'y el número de objetos será ', self.longitud )
            
        else:
            print('No ha sido posible cargar los datos del catálogo ', self.nombre)
        
    def Extraer_columna(self, DameColumna):
        self.DameColumna= self.datos(DameColumna)
        
    def mascara(self, mayor_o_menor, valor):
        
        if (mayor_o_menor=='mayor'):
            self.newdatos= self.datos > valor
            self.datos= self.newdatos
            
        elif (mayor_o_menor=='menor'):
            self.newdatos= self.datos < valor
            self.datos= self.newdatos
        
        else:
            print('La variable ingresada no es correcta')
            
        
        
        
        
        
        
        