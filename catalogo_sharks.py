"""
Created on Tue Mar 30 12:19:08 2021

@author: lpama
"""

#from astropy.table import Table
#import numpy as np
#from astropy.io import ascii
#import smatch
#import math as mt

class Catalogo():
    """
     Class to define a Catalog.
    """
    
    def __init__(self, nombre, archivo, RA, DEC):
        """
        Instantiate a Catalog.
        
        Parameters
        ----------
        
        nombre: `int`
            The name of our catalog.
        archivo: `int`
            The completely name to bring the catalog to our code.
        RA: `int`
            Name of the right ascension of our catalog.
        DEC: `int` 
            Name of the declination of our catalog.


        Returns
        -------
        None.

        """
                
        self.nombre= nombre
        self.archivo= archivo
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
        """
        With this method we start reading the catalog.
        Also if we have read properly the catalog, we will print its length.

        Returns
        -------
        None.

        """
        from astropy.table import Table
        self.datos = Table.read(self.archivo, format= 'fits')

        self.Lleno=True
        
        if self.Lleno:
            
            self.longitud= len(self.datos)
            print('El catálogo ',self.nombre,' se ha cargado correctamente y el número de objetos será ', self.longitud)
                        
        else:
            print('No ha sido posible cargar los datos del catálogo ', self.nombre)
        
    def Extraer_columna(self, DameColumna):
        """
        Here what the code does is to extract a column, which is specified in the parameter DameColumna.
        
        First, we detect if the parameter is in our catalog, and if it does, we print its length and its maximum and minimum.
        
        If not, the code prints that it can't find that column on our catalog.
        
        Parameter
        ---------
        DameColumna: `list`
            Name of column from the catalog that we want to extract.

        
        Returns
        -------
        `list`
            Here we return the column that we have named before.

        """
        import numpy as np       
        if DameColumna in self.datos.keys():
            print('La columna', DameColumna, 'tiene una longitud de', len(self.datos[DameColumna]), 'un maximo de', np.max(self.datos[DameColumna]),'y un minimo de', np.min(self.datos[DameColumna]))
                         
            return self.datos[DameColumna]

        else:
            print ('El nombre de la columna ingresada no se encuentra en el catalogo', self.nombre)


    def Match(self, ObjectCatalog):
        """
        If we have more than one catalog and we want to know how many points their have in common we will use this method.
        
        Firts of all, what we do here is a condition, if the catalogs have been read properly and if the catalogs are not the same we will start making the matching.
        
        To make the match we import from smatch library the function match, where it needs the RA and DEC of each catalog.
        
        Finally, when the match have been completed, we add the matching values to an array for each catalog.
        
        Parameter
        ---------
        
        ObjectCatalog: `int`
            Here, we write the other catalog with which we are going to make the matching.


        Returns
        -------
        None.

        """
        import smatch
                
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
        """
        To calculate the spatial area we only need the RA and DEC that is why we don't need a parameter here.

        Returns
        -------
        None.

        """
        import numpy as np
        import math as mt        
        lado_menor1=np.max(self.datos[self.RA]) - np.min(self.datos[self.RA])
        lado_mayor1=mt.sin(np.max(self.datos[self.DEC])*mt.pi/180) - mt.sin(np.min(self.datos[self.DEC])*mt.pi/180)
        self.area= (180/mt.pi)*lado_menor1*lado_mayor1 
        

        print('El valor del area de la esfera del catalogo '+self.nombre+' es',self.area,'grados cuadrados')# y el area del catalogo '+self.nombreMatch+' es',self.area2, 'grados cuadrados')

    def MainCatalog (self):
        """
        Is here where we unify both arrays that we have obtained on Match method.
        
        To make that, we create a loop that takes the name of each column and its values.
        to a single file where it is going to have the matching data of the catalogs.
        
        If the name of two columns from differents catalogs are the same we change one of them to not lead to any failure.

        Returns
        -------
        None.

        """

        from astropy.table import Table
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
            
                          
        else:
            print('No se ha podido unificar')
        
    
    def __makeCondition(self, condicion):
        """
        This method makes directly the mask to the catalog.

        Parameters
        ----------
        condicion : `list`
            This list tell the method which column and what type of condition we want to impose to make the mask.

        Returns
        -------
        `np.ndarray`
            Depending on the condition we have entered, we will obtain an array with the values that can pass through that condition.

        """
        import numpy as np
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
        """
        With the previous method we have made the mask but if the code have read the catalog properly here we determine the mask.
        Also we differentiate here how many conditions we have.

        Parameters
        ----------
        listaCondiciones : `list`
            List that specifies the column, the condition and the values to make the mask.

        Returns
        -------
        None.

        """
        #condiciones sera una lista de listas, algo como: [['MAG_AUTO','greater',12.3], ['MAGERR_AUTO','greater',0.]]
        #Los valores aceptados de condicional son: greater, greater_equal, less, less_equal, equal, not_equal
        
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
        """
        In this method we can save certain columns from our unify catalog.

        Parameters
        ----------
        listasColumnas : `list`
            List containing the name of the columns we want to save.
            Then we extract the columns data from our unify catalog to save them with its corresponding name

        Returns
        -------
        None.

        """
        from astropy.table import Table
        nombres=[]
        arrays=[]
        for name in listasColumnas:
            nombres.append(name)
            arrays.append(self.datos[name])
       
        self.t= Table(data=arrays, names=nombres)
   
    def createSample(self,format='fits', nameSample='mysample', comment=False):
        """
        Here we create the file that we have saved it before.
        
        Also we can add some comments, change the format and introduce the name of the sample.

        Parameters
        ----------
        format : `str`, optional
            Format type. The default is 'fits'.
        nameSample : `str`, optional
            Name of the file. The default is 'mysample'.
        comment : `bool`, optional
            Adding some comments. The default is False.

        Returns
        -------
        None.

        """
        from astropy.table import Table
        from astropy.io import ascii
        #format puede ser csv, ascii o fits
        comments= '# ID F294 E294 F295 E295 F296 E296 F297 E297 F298 E298'
            
        if format=='fits':
            self.t.write(nameSample+'.fits')
       
        else:            
            ascii.write(self.t, format=format, output=nameSample+'.'+format, overwrite=True, comment=comments)
        
    def giveFluxes(self, MAG, tipo):
        """
        Method that calculates flux and its error depending on the type of flux.

        Parameters
        ----------
        MAG : `list`
            List containing the name of the columns we want to extract.
        tipo : `str`
            Type of flux.

        Returns
        -------
        None.

        """
        
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
            
    
    def correctExtinction (self, nombre_mag, flujo_o_mag, banda):
        """
        Method correcting extinction.

        Parameters
        ----------
        nombre_mag : `list`
            List containing the name of the magnitude or the flux and its error.
        flujo_o_mag : `str`
            Describing if nombre_mag is flux or mag.
        banda : `str`
            String who describes the band of our magnitude or flux.

        Returns
        -------
        None.

        """
        
        rv_pars = {'G': 3.186, 'R': 2.14, 'I': 1.569, 'Z': 1.196, 'Y': 1.048, 'Ks': 0.308}
        
        EBV = self.datos['EBV_SFD98']
        
        if flujo_o_mag == 'mag':
            mag_corrected = self.datos[nombre_mag[0]] - rv_pars[banda] * EBV
            magerr_corrected = self.datos[nombre_mag[1]]
            if banda == 'Ks':
                
                self.datos['APERMAG3_CORRECTED']= mag_corrected
                self.datos['APERMAG3ERR_CORRECTED']= magerr_corrected
                
            else:
                self.datos['MAG_AUTO_'+banda+'_CORRECTED']= mag_corrected
                self.datos['MAGERR_AUTO_'+banda+'_CORRECTED']= magerr_corrected
                
                               

        elif flujo_o_mag == 'flux':
            flux_corrected = self.datos[nombre_mag[0]] * 10**(rv_pars[banda] * EBV/2.5)
            fluxerr_corrected = self.datos[nombre_mag[1]] * 10**(rv_pars[banda] * EBV/2.5)
            if banda== 'Ks':
                self.datos['FLUX_APER_CORRECTED']= flux_corrected
                self.datos['FLUXERR_APER_CORRECTED']= fluxerr_corrected
                
            else:
                self.datos['FLUX_AUTO_'+banda+'_CORRECTED']= flux_corrected
                self.datos['FLUXERR_AUTO_'+banda+'_CORRECTED']= fluxerr_corrected
        	
    
    def estado(self):
        """
        Method who reports the status of the code

        Returns
        -------
        None.

        """
        
        
        if self.Lleno:
            
            print('El catalogo '+self.nombre+' se ha cargado correctamente')
            
            if  self.match:
                print('Se ha realizado correctamente el match entre los catalogos '+self.nombre+' y '+self.nombreMatch)
                
            else:
                print ('No se ha podido realizar el match entre los catalogos '+self.nombre+' y '+self.nombreMatch)
            
            if self.mask:
                print('Tras haber realizado correctamente la mascara tenemos ',len(self.datos),' objetos, siendo el area espacial del catalogo ',self.area)
                
            else:
                print('No se ha podido realizar la mascara correctamente')
            
        else:
            print('No ha sido posible leer el catalogo '+self.nombre+' correctamente')

        
        
           
