import numpy as np
import matplotlib.pyplot as pl
#import pidly
#idl=pidly.IDL('/opt/ioa/idl/idl82/bin/idl')
#idl('.run ../localpopcalc.pro')
import pickle
#import pyfits as pf
from scipy.constants import h,c,k
import time
import uuid
import os
import sys
sys.path.insert(0, '/sma/SMAusers/lmatra/CO_excitation')
import elrovibpopcalc_CO
reload(elrovibpopcalc_CO)
import elrovibpopcalc_CO_noISRFUV
reload(elrovibpopcalc_CO_noISRFUV)
import scipy.ndimage
from scipy.ndimage.filters import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.ndimage.interpolation import shift

rcParams['font.family'] = 'Times New Roman'
font={'size':25}	
rc('font', **font)
rc('axes', linewidth=2)
rc('xtick.major', width=2, size=8)
rc('xtick.minor', width=2, size=4)
rc('ytick.major', width=2, size=8)
rc('ytick.minor', width=2, size=4)
rc('xtick', labelsize=22)
rc('ytick', labelsize=22)
matplotlib.rcParams.update({'figure.autolayout': True})

meanintfilename=None

#Read star flux from Grant's sdb file.
starfluxfilename='./HD172555_starfluxJyfromEarth.npy'
import json
import matplotlib.pyplot as pl
import numpy as np

filename='./SED_HD172555.json'
with open(filename, 'r') as f:
	datastore=json.load(f)
starspecwav=np.asarray(datastore['model_spectra'][0]['wavelength']) #micron
starspecfnujy=np.asarray(datastore['model_spectra'][0]['fnujy']) #Jy
#starspecfnujy[(starspecwav<0.3)]=0.0
np.save(starfluxfilename, [starspecwav, starspecfnujy])
#Check
pl.figure()
pl.plot(starspecwav, starspecfnujy)
pl.xscale('log')
pl.yscale('log')


#Measured upper limit in Jy km/s
measureduplim=0.120 #25.0*1e-3

#Define electron density and temperature arrays to loop over
ncollps=10**(np.arange(16)-5.0)
temps=np.asarray([20.0,50.0,100.0,169.0,250.0])	#10**(np.arange(12)/11.0*(np.log10(1000.0)-np.log10(5.0))+np.log10(5))

masslims=np.zeros((ncollps.size, temps.size))
fpoplev1=np.zeros((ncollps.size, temps.size))
fpoplev2=np.zeros((ncollps.size, temps.size))
pops=np.zeros((ncollps.size, temps.size,540))
for i in np.arange(ncollps.size):
	for j in np.arange(temps.size):
		masslims[i,j]=elrovibpopcalc_CO.elrovibpopcalc_COfunc(tgas=temps[j], ncollp=ncollps[i], datafile='/sma/SMAusers/lmatra/CO_excitation/COvibrot_Lambda_nelec2_nvib9_nrot30.dat', jnufilestar=starfluxfilename, linewithntrcoll=114067, distpc=28.3, distau=7.5, levupint=2, si=0, fluxJykmsinp=measureduplim)
		#fpoplev2[i,j]=elrovibpopcalc_CO.elrovibpopcalc_COfunc(tgas=temps[j], ncollp=ncollps[i], datafile='/sma/SMAusers/lmatra/CO_excitation/COvibrot_Lambda_nelec2_nvib9_nrot30.dat', jnufilestar=starfluxfilename, linewithntrcoll=114067, distpc=132.1, distau=58.4, levupint=2)
		#fpoplev1[i,j]=elrovibpopcalc_CO.elrovibpopcalc_COfunc(tgas=temps[j], ncollp=ncollps[i], datafile='/sma/SMAusers/lmatra/CO_excitation/COvibrot_Lambda_nelec2_nvib9_nrot30.dat', jnufilestar=starfluxfilename, linewithntrcoll=114067, distpc=132.1, distau=58.4, levupint=1)
print 'NB No disk continuum radiation field considered at the moment!'
print masslims.min(), masslims.max()

plotmasslims=1
if plotmasslims:
	pl.figure(figsize=(12,7))
	#pl.plot(ncollps, masslims, ls='dotted', lw=2.0)
	for i in np.arange(temps.size):
		pl.plot(ncollps, masslims[:,i], linewidth=3.0, ls='solid', label=r'T$_{\rm k}$ = '+str(int(temps[i]))+' K', alpha=0.4)
	pl.xscale('log')
	pl.yscale('log')
	pl.xlabel(r'e$^{-}$ density (cm$^{-3}$)')
	pl.ylabel(r'CO mass (M$_{\oplus}$)')
	pl.legend(frameon=0)
	ax=pl.gca()
	pl.text(0.04,0.12, 'HD172555', transform=ax.transAxes, fontsize=28)
	pl.savefig('HD172555_21_COmass_vs_ncoll_Tk.pdf')
	pl.close()
	os.system('evince HD172555_21_COmass_vs_ncoll_Tk.pdf &')
