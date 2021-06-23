import emcee
import math
import statistics
import time
from gas_ring_fitting import *
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
import pickle

ndim,nwalkers = 5,50
nsteps = 2000

pi = np.pi

ring_mid = 8# * u.AU
ring_wid = 2# * u.AU
inclination = math.radians(103.5)
stellar_velocity = 2000# * u.m/u.s #m/s

hd172555 = gasringfit(1.71,2.12)

cube_file = "../../data/ALMA_raw/hd172555/hd172555_CO_BARY_aftstatwt_nofarants_natural_cleaned.fits"
mom0_file = "../../data/ALMA_raw/hd172555/hd172555_CO_BARY_aftstatwt_nofarants_natural_cleaned_mom0.fits"

hd172555.return_alma_data(cube_file,mom0_file)

loboundtarg = hd172555.veltofreq(15000,hd172555.line_freq.value,stellar_velocity)
hiboundtarg = hd172555.veltofreq(-15000,hd172555.line_freq.value,stellar_velocity)
lobound = hd172555.find_nearest(hd172555.data_freqarray,loboundtarg)
hibound = hd172555.find_nearest(hd172555.data_freqarray,hiboundtarg)

yerr_1 = (statistics.stdev(hd172555.data_intensities[0:lobound])+statistics.stdev(hd172555.data_intensities[hibound:-1]))/2
yerr = np.array([yerr_1 for i in hd172555.data_intensities])
print("The yerr is " + str(yerr_1))
print(" ")

def log_likelihood(theta,x,data,yerr):
    i,r_mid,r_wid,v_s,scal = theta
    model = hd172555.keplerianring(i,r_mid,r_wid,v_s,scal)
    sigma2 = yerr**2
    return(-0.5*np.sum((data - model)**2/sigma2 + np.log(2*pi*sigma2)))

def log_prior_gaussi(theta):
    i,r_mid,r_wid,v_s,scal = theta

    if not (0.25<r_mid<12 and \
            0.05<r_wid<2*r_mid and \
            -5e5<v_s<5e5 and \
            0<scal<1e7 and \
            i>math.radians(90)):
        return(-np.inf)
    mu = math.radians(103.5)
    sigma = math.radians(7)
    return(np.log(1.0/(np.sqrt(2*pi)*sigma)) - 0.5*(i-mu)**2/sigma**2)

def log_prior_flati(theta):
    i,r_mid,r_wid,v_s,scal = theta
    if not (math.radians(90)<i<math.radians(180) and \
            0.25<r_mid<12 and \
            0.05<r_wid<2*r_mid and \
            -5e5<v_s<5e5 and \
            0<scal<1e7):
        return(-np.inf)
    return(0.0)

def log_probability_gaussi(theta,x,data,yerr):
    lp = log_prior_gaussi(theta)
    if not np.isfinite(lp):
        return(-np.inf)
    return(lp+log_likelihood(theta,x,data,yerr))

def log_probability_flati(theta,x,data,yerr):
    lp = log_prior_flati(theta)
    if not np.isfinite(lp):
        return(-np.inf)
    return(lp+log_likelihood(theta,x,data,yerr))

"""
Integrated flux scaling MCMC runs
"""
print("Beginning the integrated flux scaling MCMC runs")

initial = np.array([inclination,ring_mid,ring_wid,stellar_velocity,hd172555.data_integflux])
pos = [initial+1e-6*np.random.randn(ndim) for i in range(nwalkers)]

print("integrated flux - gaussian i:")

sampler_gi = emcee.EnsembleSampler(nwalkers,
                                   ndim,
                                   log_probability_gaussi,
                                   args = (hd172555.data_freqarray,hd172555.data_intensities,yerr),
                                   threads = 7)
start_time = time.time()
for i, result in enumerate(sampler_gi.sample(pos, iterations=nsteps)):
    if (i+1) % 100 == 0:
        print("{0:5.1%}".format(float(i+1) / nsteps) + \
              "     Minutes Elapsed: "+str((time.time() - start_time)/60))
pickle.dump(sampler_gi.chain,
            open("mcmc_outputs/5paramfit_fluxscaling_hd172555_gaussiani_500x500_unbinned.p", "wb"))

print("integrated flux - flat i:")

sampler_fi = emcee.EnsembleSampler(nwalkers,
                                   ndim,
                                   log_probability_flati,
                                   args = (hd172555.data_freqarray,hd172555.data_intensities,yerr),
                                   threads = 7)
start_time = time.time()
for i, result in enumerate(sampler_fi.sample(pos, iterations=nsteps)):
    if (i+1) % 100 == 0:
        print("{0:5.1%}".format(float(i+1) / nsteps) + \
              "     Minutes Elapsed: "+str((time.time() - start_time)/60))
pickle.dump(sampler_fi.chain,
            open("mcmc_outputs/5paramfit_fluxscaling_hd172555_flati_500x500_unbinned.p", "wb"))
