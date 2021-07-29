import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
import os
import sys
sys.path.append('/Users/lmatra/radmc-3d/version_0.41/python')
import problem_setup_12CO_powlaw
import spectres
import run_radmc_tautest_fct
c = 299792458.0	   #m/s

#Beauty plot parameters
#mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif']=['Arial']
font={'size':18}
mpl.rc('font', **font)
mpl.rc('axes', linewidth=2)
mpl.rc('xtick.major', width=2, size=10)
mpl.rc('xtick.minor', width=2, size=5)
mpl.rc('ytick.major', width=2, size=10)
mpl.rc('ytick.minor', width=2, size=5)
mpl.rc('xtick', labelsize=16)
mpl.rc('ytick', labelsize=16)
mpl.rcParams.update({'figure.autolayout': True})
mpl.rcParams['pdf.fonttype']=42

#Read in data
freq_obs, sp_obs_mJy = np.load('./unbinned_spectrum_observed.npy')
linefreq=2.30538 #100s of GHz

v_obs=-(freq_obs-linefreq)/linefreq*c/1e3 
vstar=2.3 #km/s
rms_spectrum=np.sqrt(np.mean(sp_obs_mJy[np.abs(v_obs)>20.0]**2.0))

run_models=False
if run_models:
    #Within emcee function:
    logtemp=np.arange(200)/199.0*0.9+1.4
    logmass=np.arange(200)/199.0*5.0-5.5
    chisquare=np.zeros((logtemp.size, logmass.size))
    for i in np.arange(logtemp.size):
         for j in np.arange(logmass.size):
            print(i,j)
            modelvels, modelspec = run_radmc_tautest_fct.run_radmc_tautest_fct((10**logtemp[i]), logmass[j])

            modelspec_rebinned = spectres.spectres(v_obs[::-1], modelvels[::-1]-vstar, modelspec[::-1]*1000.0,spec_errs=None, fill=0.0, verbose=False)
            chisquare[i,j] = np.sum(((sp_obs_mJy-modelspec_rebinned)/rms_spectrum)**2.0)/2.667

    np.save('logtemplogmass_200temp_200mass_highres.npy', [logtemp,logmass])
    np.save('chisquare_200temp_200mass_highres.npy', chisquare)
else:
    logtemp, logmass = np.load('logtemplogmass_200temp_200mass_highres.npy')
    chisquare = np.load('chisquare_200temp_200mass_highres.npy')   

fig = pl.figure(figsize=(17,7.3))
ax = fig.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [0.9, 1.1]})
im = ax[1].imshow(((chisquare-np.min(chisquare))), origin='lower', extent=(logmass[0],logmass[-1],logtemp[0], logtemp[-1]), aspect=5.0, cmap='inferno', vmax=30.0)
im.cmap.set_over('white')
ax[1].plot(np.log10(1.4e-5),np.log10(169.0),'o',color='red')
ax[1].plot(np.log10(1.5e-2),np.log10(35.0),'o',color='cyan')




#Normalised probability map (where P \propto exp(-0.5*chi), so Prel=P/Pmin=exp(-chi+chi_min)=exp(-deltachi2)), so peak probability is 1
levels=1.0-np.exp(-0.5*np.asarray([1.0,2.0,3.0,4.0])**2)
normprob=np.exp(-0.5*(chisquare-np.min(chisquare)))
H=normprob
Hflat=H.flatten()
inds=np.argsort(Hflat)[::-1]
Hflat=Hflat[inds]
sm=np.cumsum(Hflat)
sm/=sm[-1]
V=np.empty(len(levels))
for i, v0 in enumerate(levels):
    try:
        V[i] = Hflat[sm<=v0][-1]
    except IndexError:
        V[i] = Hflat[0]
V.sort()
m = np.diff(V) == 0
while np.any(m):
    V[np.where(m)[0][0]]*= 1.0-1e-4
    m = np.diff(V) == 0
V.sort()    

cts = ax[1].contour(H, V, origin='lower', linewidths=3.0, alpha=0.7, linestyles='dotted', extent=(logmass[0],logmass[-1],logtemp[0], logtemp[-1]), colors='white')
cbar = fig.colorbar(im, ax=ax[1], label=r'$\Delta \chi^2$', fraction=0.05)

ax[1].set_yticks([np.log10(25),np.log10(50),np.log10(75),np.log10(100),np.log10(125),np.log10(150),np.log10(175),np.log10(200)])
ax[1].set_yticklabels(['25','50','75','100','125','150','175','200'])
ax[1].set_xticklabels(['dum',r'10$^{-5}$',r'10$^{-4}$',r'10$^{-3}$',r'10$^{-2}$',r'10$^{-1}$'])
ax[1].set_xlabel(r'CO Mass [M$_{\oplus}$]')
ax[1].set_ylabel('Temperature [K]')
ax[1].text(np.log10(1.5e-1),np.log10(175.0), 'b', fontweight='bold')


newchansize=0.5 #km/s
newvelarr=v_obs[0]-0.5*(v_obs[1]-v_obs[0])-np.arange(np.floor(-(v_obs[-1]-v_obs[0])/newchansize))*newchansize

modelvels, modelspec_warm = run_radmc_tautest_fct.run_radmc_tautest_fct(169.0, np.log10(1.4e-5))
modelvels, modelspec_cold = run_radmc_tautest_fct.run_radmc_tautest_fct(35.0, np.log10(1.5e-2))
modelspec_warm_rebinned = spectres.spectres(newvelarr[::-1], modelvels[::-1]-vstar, modelspec_warm[::-1]*1000.0,spec_errs=None, fill=0.0, verbose=False)
modelspec_cold_rebinned = spectres.spectres(newvelarr[::-1], modelvels[::-1]-vstar, modelspec_cold[::-1]*1000.0,spec_errs=None, fill=0.0, verbose=False)
dataspec=spectres.spectres(newvelarr[::-1], v_obs[::-1], sp_obs_mJy[::-1],spec_errs=None, fill=0.0, verbose=True)
dataspec=dataspec[::-1]


ax[0].plot(newvelarr,dataspec, drawstyle='steps-mid', alpha=0.4, color='black', linewidth=2.5, label='Data')
ax[0].fill_between(newvelarr[::-1]+vstar, modelspec_warm_rebinned, alpha=0.15, color='red', label=r'Warm, low mass, optically thin model', linewidth=0.0)
ax[0].fill_between(newvelarr[::-1]+vstar, modelspec_cold_rebinned, alpha=0.15, color='blue', label=r'Cold, high mass, optically thick model', linewidth=0.0)
ax[0].plot(newvelarr[::-1]+vstar, modelspec_warm_rebinned, alpha=0.5, drawstyle='steps-mid', color='red', linestyle='dashed', linewidth=2.5)#, label=r'RADMC-3D, T = 169 K, M$_{\rm CO}$=10$^{-5}$ M$_{\oplus}$')
ax[0].plot(newvelarr[::-1]+vstar, modelspec_cold_rebinned, alpha=0.5, drawstyle='steps-mid', color='blue', linestyle='dotted', linewidth=2.5)#, label=r'RADMC-3D, T = 35 K, M$_{\rm CO}$=2$\times$10$^{-3}$ M$_{\oplus}$')
ax[0].set_ylabel('Flux (mJy)')
ax[0].set_xlabel('Velocity (km/s)')
ax[0].set_xlim(-28,30.3)
ax[0].set_ylim(-15,20)
ax[0].legend(frameon=False, loc=4)
ax[0].text(-26.0,18.0, 'a', fontweight='bold')
fig.savefig('Modelspectracomparison_pluschisquaremap_sans.pdf')
os.system('open Modelspectracomparison_pluschisquaremap_sans.pdf')
