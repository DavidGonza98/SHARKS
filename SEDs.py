# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 15:48:01 2021

@author: lpama
"""

import numpy as np
import matplotlib.pyplot as plt


sedS0= np.loadtxt("Im_B2004a.sed")
sedSa= np.loadtxt("SB2_B2004a.sed")
sedSb= np.loadtxt("SB3_B2004a.sed")
sedSB1= np.loadtxt("Sbc_B2004a.sed")
sedSB2= np.loadtxt("Scd_B2004a.sed")
sedSc= np.loadtxt("ssp_25Myr_z008.sed")
sedSd= np.loadtxt("ssp_5Myr_z008.sed")
sedEll= np.loadtxt("El_B2004a.sed")

column1S0=sedS0[:, 0]
column2S0=sedS0[:, 1]

column1Sa=sedSa[:, 0]
column2Sa=sedSa[:, 1]

column1Sb=sedSb[:, 0]
column2Sb=sedSb[:, 1]

column1SB1=sedSB1[:, 0]
column2SB1=sedSB1[:, 1]

column1SB2=sedSB2[:, 0]
column2SB2=sedSB2[:, 1]

column1Sc=sedSc[:, 0]
column2Sc=sedSc[:, 1]

column1Sd=sedSd[:, 0]
column2Sd=sedSd[:, 1]

column1Ell=sedEll[:, 0]
column2Ell=sedEll[:, 1]



from matplotlib.pyplot import figure

figure(figsize=(12, 6), dpi=80)
plt.plot(column1Sd, np.log10(column2Sd), 'b--', label='Ssp_5', alpha=0.6)
plt.plot(column1Sc, np.log10(column2Sc), '--', label='Ssp_25', alpha=0.6)
plt.plot(column1Sa, np.log10(column2Sa), 'k', label='SB2')
plt.plot(column1SB1, np.log10(column2SB1), label='SB1')
plt.plot(column1Sb, np.log10(column2Sb), label='SB3')
plt.plot(column1S0, np.log10(column2S0), 'm', label='Im')
plt.plot(column1SB2, np.log10(column2SB2), 'r', label='Scd')
plt.plot(column1Ell, np.log10(column2Ell), label='El')
#plt.set_yscale()




plt.legend()
plt.ylim(0, 2)
plt.xlim(0, 10000)
plt.title('Plantillas de SED')
plt.xlabel('Frecuencia [${\AA}$]')
plt.ylabel( r"Flujo [$erg/s/cm^{2}/{\AA}$]")
plt.savefig('SEDs.png')




