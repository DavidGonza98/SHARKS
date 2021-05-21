# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 22:09:02 2021

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


import catalogo_sharks, astrometria_poo, Fotometria, PhotoZ

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
'''
des=catalogo_sharks.Catalogo('DES', 'des_in_field.fits', 'RA', 'DEC')
sharks_sgpe_des=catalogo_sharks.Catalogo('Sharks_sgpe', 'sharks_sgpe.fits', 'RA', 'DEC')

#PARA AURE des=Catalogo('DES', '/home/acarnero/sharks/dr1/suplement/des_in_field.fits', 'RA', 'DEC')
#PARA AURE sharks_sgpe_des=Catalogo('Sharks_sgpe', '/home/acarnero/Documents/tfg/sharks_sgpe.fits', 'RA', 'DEC')

sharks_sgpe_des.LeerArchivo()
des.LeerArchivo()


sharks_sgpe_des.Match(des)
 
sharks_sgpe_des.MainCatalog()
#sharks_sgpe_des.saveSample(['COADD_OBJECT_ID', 'FLUX_AUTO_G', 'FLUXERR_AUTO_G', 'FLUX_AUTO_R', 'FLUXERR_AUTO_R', 'FLUX_AUTO_I', 'FLUXERR_AUTO_I', 'FLUX_AUTO_Z', 'FLUXERR_AUTO_Z', 'FLUX_AUTO_Y', 'FLUXERR_AUTO_Y', 'PETROFLUX', 'PETROFLUXERR', 'EBV_SFD98', 'MULTIFRAMEID', 'SEQNUM'])
#sharks_sgpe_des.createSample(format='fits', nameSample='sharks_des1')
sharks_sgpe_des.estado()
sharks_sgpe_des.AreaEspacial()

sharks_sgpe_des.correctExtinction(['FLUX_AUTO_G','FLUXERR_AUTO_G'],'flux','G')
sharks_sgpe_des.correctExtinction(['FLUX_AUTO_R','FLUXERR_AUTO_R'],'flux','R')
sharks_sgpe_des.correctExtinction(['FLUX_AUTO_I','FLUXERR_AUTO_I'],'flux','I')
sharks_sgpe_des.correctExtinction(['FLUX_AUTO_Z','FLUXERR_AUTO_Z'],'flux','Z')
sharks_sgpe_des.correctExtinction(['FLUX_AUTO_Y','FLUXERR_AUTO_Y'],'flux','Y')
sharks_sgpe_des.correctExtinction(['PETROFLUX','PETROFLUXERR'],'flux','Ks')

sharks_sgpe_des.saveSample(['COADD_OBJECT_ID','RA','DEC', 'FLUX_AUTO_CORRECTED_G', 'FLUXERR_AUTO_CORRECTED_G', 'FLUX_AUTO_CORRECTED_R', 'FLUXERR_AUTO_CORRECTED_R', 'FLUX_AUTO_CORRECTED_I', 'FLUXERR_AUTO_CORRECTED_I', 'FLUX_AUTO_CORRECTED_Z', 'FLUXERR_AUTO_CORRECTED_Z', 'FLUX_AUTO_CORRECTED_Y', 'FLUXERR_AUTO_CORRECTED_Y', 'FLUX_APER_CORRECTED', 'FLUXERR_APER_CORRECTED'])
sharks_sgpe_des.createSample(format='basic', nameSample='Galaxias_dessharks')
'''
'''
sharks_sgpe_des_astrometria=astrometria_poo.astrometria(sharks_sgpe_des)
sharks_sgpe_des_astrometria.Distancia_Angular()
#sharks_sgpe_des_astrometria.Histograma('Distancia angular')#, sharks_sgpe_des, des)
sharks_sgpe_des_astrometria.DefineVariables()

sharks_sgpe_des_astrometria.plot('X')
#sharks_sgpe_des_astrometria.plot('Y')
#sharks_sgpe_des_astrometria.plot('L')
#sharks_sgpe_des_astrometria.plot('B')
#sharks_sgpe_des_astrometria.plot('RA')
#sharks_sgpe_des_astrometria.plot('DEC')

sharks_sgpe_des_astrometria.plot('L')
'''

results_eazy = PhotoZ.PhotoZ('eazy', ['gals_dessharks_1.txt', 'gals_dessharks_2.txt', 'gals_dessharks_3.txt', 'gals_dessharks_4.txt', 'gals_dessharks_5.txt'], 'sharks')

results_eazy.getStats()  
results_eazy.getStats_Bin([0.,0.1,0.2,0.3,0.4,0.5])

