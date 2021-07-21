import astropy.io.fits as pf
from scipy.stats import chi2
import numpy as np
import matplotlib.pyplot as plt

data_ = pf.open('master_cat_5s.fits')[1].data

chi = data_['CHI_BEST_DOY']
type_ = data_['TYPE_DOY']
nbands = data_['NBAND_USED_DOY']

mask = (chi<100) * (nbands>3) * (type_!='S')
chi = chi[mask]
type_ = type_[mask]
nbands = nbands[mask]

nbands3_A = len(nbands[(nbands==3)]) / float(len(nbands)) 
nbands4_A = len(nbands[(nbands==4)]) / float(len(nbands))
nbands5_A = len(nbands[(nbands==5)]) / float(len(nbands))
nbands6_A = len(nbands[(nbands==6)]) / float(len(nbands))

bins = np.linspace(0.001, 50, 100)

chi21, chi22, chi23, chi24 = [], [], [], []
CHI2ALL = []
for i in bins:
    chi21.append(chi2.pdf(i, 1))
    chi22.append(chi2.pdf(i, 2))
    chi23.append(chi2.pdf(i, 3))
    chi24.append(chi2.pdf(i, 4))

for i,kk in enumerate(chi21):
    CHI2ALL.append(nbands3_A*chi21[i]+ nbands4_A*chi22[i]+nbands5_A*chi23[i]+nbands6_A*chi24[i])


fig, ax = plt.subplots()

ax.hist(chi, bins=1000, density=True)
plt.xlim([0,30])
plt.title('DES Only Poggianti')

print (len(chi[(chi>20.)]))

plt.legend()
plt.xlabel('$\chi^{2}$ Lephare')
plt.ylabel('Density')
ax.plot(bins,CHI2ALL,'r',label='Expected $\chi^{2}$')

plt.legend(fontsize=10,frameon=False)
fig.tight_layout()

plt.savefig('chi2_DOY.png',bbox_inches='tight')
