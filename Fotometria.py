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



class Catalogo():
    
    def __init__(self, nombre, archivo, RA, DEC):
        self.nombre= nombre
        self.archivo= archivo
        
        self.lineas=False
        self.Lleno=False
        self.RA=RA
        self.DEC=DEC
        self.nside=4096 
        self.maxmatch=1 
        self.radius= 1/3600.
        self.match=False
        self.mask=False
    
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
            #print(describa si ha sido leido o no, se se ha hecho un match o no
            #el numero de lineas), si he hecho una mascara
            #si hay match, print del nombre que hemos hecho
            print('La linea se ha cargado correctamente')
            
            self.Extraer_columna(self.RA)
            self.Extraer_columna(self.DEC)
           
            
            
        else:
            print('No ha sido posible cargar las lineas')

    def Match(self, ObjectCatalog):
        
        if (self.Lleno) and self.nombre!=ObjectCatalog.nombre:
            ra1_matched=self.datos[self.RA]
            dec1_matched=self.datos[self.DEC]
            
            self.datos2 = ObjectCatalog.datos
            self.RA2 = ObjectCatalog.RA
            self.DEC2 = ObjectCatalog.DEC
            ra2_matched=self.datos2[self.RA2]
            dec2_matched=self.datos2[self.DEC2]
            #print('La longitud de', self.RA, 'es', len(ra1_matched), 'y la de', self.DEC, 'es', len(dec1_matched) )
            #print('La longitud de', ObjectCatalog[2], 'es', len(ra2_matched), 'y la de', ObjectCatalog[3], 'es', len(dec2_matched), 'donde estos datos pertenecen al catalogo', ObjectCatalog[0])
            self.matches = smatch.match(ra1_matched, dec1_matched, self.radius, ra2_matched, dec2_matched, nside=self.nside, maxmatch=self.maxmatch)
            
            self.assoc1 = self.datos[self.matches['i1']]
            self.assoc2 = self.datos2[self.matches['i2']]


            print('La longitud del catalogo ya habiendo realizado el matched es', len(self.assoc1))
            self.match=True
            self.nombreMatch = ObjectCatalog.nombre 
            self.MainCatalog()
            
        else:
            print('No se ha podido hacer el match')
            
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
        
    def giveFluxes(self, MAG):
        
        
        flux= 10**((48.6 - (self.Extraer_columna(MAG[0])))/2.5)

        fluxerr=self.Extraer_columna(MAG[1]) * flux/1.086
        
        print('El valor del flujo es', flux,'y del error',fluxerr)
        self.datos['FLUX_'+MAG[0]] = flux
        
        self.datos['FLUXERR_'+MAG[0]] = fluxerr 

#sharks.Match(twomass)
#Catalogo ya tendra dentro
#astro = astrometria(sharks)

#astro.Distancia_angular()

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
        #Lleno=False
        #self.matched=False
        
    
    def Distancia_Angular(self):
        #sharks.Match(self.datos)
        self.cosgamma = np.zeros(len(self.datos))
        ra2_matched=self.datos[self.ra2]
        dec2_matched=self.datos[self.dec2]
        ra1_matched=self.datos[self.ra1]
        dec1_matched=self.datos[self.dec1]
        
        for i in range(len(self.datos)):
            self.gamma = mt.cos(90-dec1_matched[i])*mt.cos(90-dec2_matched[i])+mt.sin(90-dec1_matched[i])*mt.sin(90-dec2_matched[i])*mt.cos(ra1_matched[i]-ra2_matched[i])
            if self.gamma <= 1:
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
        


class fotometria():
    
    def __init__(self,ObjectCatalog):
        
        self.datos = ObjectCatalog.datos
                             
        
   
    def ajuste_lineal(self, magnitud1, magnitud2, ObjectCatalog, ObjectCatalog1, quitarOutliers=False):
        
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
        #plt.savefig(+title+'.png')
        print ('El valor de la pendiente del ajuste lineal es', self.p[0], 'y su ordenada en el origen', self.p[1]  )
        
        
            
    def errorfunc(self, magnitud1, magnitud2, errormagnitud1, ObjectCatalog, ObjectCatalog1, quitarOutliers=False, pendienteuno=False):
        
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
        
'''
sharks=Catalogo('Sharks', 'Sharks_sgp_e_2_cat_small.fits', 'ALPHA_J2000', 'DELTA_J2000')


sharks.LeerArchivo()
#ra = sharks.Extraer_columna('ALPHA_J2000')
sharks.mascara([['MAG_AUTO','greater_equal',12.3],['MAGERR_AUTO','greater',0.]])

sharks.estado()
sharks.MainCatalog()

twomass = Catalogo('2mass', '2mass.fit', 'RAJ2000', 'DEJ2000')
twomass.LeerArchivo()
twomass.estado()

sharks.Match(twomass)
sharks.giveFluxes(['MAG_AUTO','MAGERR_AUTO'])
#sharks.createSample(['MAG_AUTO','MAGERR_AUTO','FLUX_MAG_AUTO','FLUXERR_MAG_AUTO'],nameSample='sharks_twomass')



sharks.mascara([['MAG_AUTO','greater_equal',12.3],['MAGERR_AUTO','greater',0.], ['e_Kmag','greater',0.], ['Kmag', 'greater_equal', 12.7], ['Kmag', 'less_equal', 14.5]])
print('Despues de mascarar', len(sharks.datos))

sharks.giveFluxes(['MAG_AUTO','MAGERR_AUTO'])      


sharks_astrometria= astrometria(sharks)

sharks_astrometria.Distancia_Angular()       
 

sharks_fotometria= fotometria(sharks)
sharks_astrometria.Histograma('Distancia angular')

sharks_fotometria.ajuste_lineal('Kmag', 'MAG_AUTO', sharks, twomass)

sharks_fotometria.errorfunc('Kmag', 'MAG_AUTO', 'e_Kmag',sharks, twomass, True, True)
'''


'''
sharks_sgpe=Catalogo('Sharks_sgpe', 'sharks_sgpe.fits', 'RA', 'DEC')
two_mass=Catalogo('2mass_infield', '2mass_in_field.fits', 'RAJ2000', 'DEJ2000')

sharks_sgpe.LeerArchivo()
two_mass.LeerArchivo()

#two_mass.Extraer_columna('RAJ2000')
#sharks_sgpe.Extraer_columna('RA')

sharks_sgpe.estado()
two_mass.estado()

two_mass.Match(sharks_sgpe)

two_mass.MainCatalog()
#sharks_sgpe.mascara([['MAG_AUTO','greater_equal',12.3],['MAGERR_AUTO','greater',0.]])
two_mass.mascara([['APERMAG3','greater_equal',12.3],['APERMAG3ERR','greater',0.], ['e_Kmag','greater',0.], ['Kmag', 'greater_equal', 12.7], ['Kmag', 'less_equal', 14.5]])

two_mass.giveFluxes(['APERMAG3','APERMAG3ERR']) 

two_mass_astrometria= astrometria(two_mass)

two_mass_astrometria.Distancia_Angular()
two_mass_astrometria.Histograma('Distancia angular', two_mass, sharks_sgpe)

two_mass_fotometria= fotometria(two_mass)



two_mass_fotometria.ajuste_lineal('Kmag', 'APERMAG3', two_mass, sharks_sgpe)

two_mass_fotometria.errorfunc('Kmag', 'APERMAG3', 'e_Kmag', two_mass, sharks_sgpe, True, True)
'''

des=Catalogo('DES', 'des_in_field.fits', 'RA', 'DEC')
sharks_sgpe_des=Catalogo('Sharks_sgpe', 'sharks_sgpe.fits', 'RA', 'DEC')

#PARA AURE des=Catalogo('DES', '/home/acarnero/sharks/dr1/suplement/des_in_field.fits', 'RA', 'DEC')
#PARA AURE sharks_sgpe_des=Catalogo('Sharks_sgpe', '/home/acarnero/Documents/tfg/sharks_sgpe.fits', 'RA', 'DEC')

sharks_sgpe_des.LeerArchivo()
des.LeerArchivo()

sharks_sgpe_des.estado()
des.estado()

sharks_sgpe_des.Match(des)
 
sharks_sgpe_des.MainCatalog()
#sharks_sgpe_des.saveSample(['COADD_OBJECT_ID', 'FLUX_AUTO_G', 'FLUXERR_AUTO_G', 'FLUX_AUTO_R', 'FLUXERR_AUTO_R', 'FLUX_AUTO_I', 'FLUXERR_AUTO_I', 'FLUX_AUTO_Z', 'FLUXERR_AUTO_Z', 'FLUX_AUTO_Y', 'FLUXERR_AUTO_Y', 'PETROFLUX', 'PETROFLUXERR', 'EBV_SFD98', 'MULTIFRAMEID', 'SEQNUM'])
#sharks_sgpe_des.createSample(format='fits', nameSample='sharks_des1')

sharks_sgpe_des_astrometria=astrometria(sharks_sgpe_des)
sharks_sgpe_des_astrometria.Distancia_Angular()
sharks_sgpe_des_astrometria.Histograma('Distancia angular')#, sharks_sgpe_des, des)

      
