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
from scipy import optimize, stats
import pandas as pd

class astrometria():
    
    def __init__(self, ObjectCatalog):
        print('hola')
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
        
        self.ra2_matched=self.datos[self.ra2]
        self.dec2_matched=self.datos[self.dec2]
        self.ra1_matched=self.datos[self.ra1]
        self.dec1_matched=self.datos[self.dec1]
    
       
    
    def Distancia_Angular(self):
        #sharks.Match(self.datos)
        self.cosgamma = np.zeros(len(self.datos))

        
        for i in range(len(self.datos)):
            self.gamma = mt.cos(90-self.dec1_matched[i])*mt.cos(90-self.dec2_matched[i])+mt.sin(90-self.dec1_matched[i])*mt.sin(90-self.dec2_matched[i])*mt.cos(self.ra1_matched[i]-self.ra2_matched[i])
            self.cosgamma[i] = (mt.acos(self.gamma)*3600)
        return self.cosgamma
    
    def Histograma(self, tipo):
    
        plt.figure(1)
        if tipo == 'Distancia angular':
            
            sns.displot(self.cosgamma, kde=True)
            plt.title('Histograma de distancia angular')
            title='Histograma_Distancia_Angular' + self.nombre+ '_y_'  + self.nombreMatch
            plt.xlabel('Segundos de Arco')
            plt.savefig(title+'.png')
            
        
        elif tipo == 'Jointplot':
            return sns.jointplot(x=self.KMAG, y=self.MAG_AUTO, kind="reg", truncate=False, xlim=(5, 25), ylim=(5, 25), color="m", height=7)
        
    def plot(self):
        
        dra=self.ra1_matched-self.ra2_matched
        ddec=self.dec1_matched-self.dec2_matched
        plt.plot(dra*3600, ddec*3600, '.')
        plt.plot(0,0, 'r+')
        plt.xlabel('dRA [arcsec]')
        plt.ylabel('dDEC [arcsec]')
        title= self.nombre+ '-'  + self.nombreMatch+' pos offsets'
        plt.title(title)
        plt.savefig(title+'.png')     
        
        
        
        
        
        
        
        
        
        

         
    

















