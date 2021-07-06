def run_radmc_tautest_fct(temp, logmass):
	import numpy as np
	import matplotlib.pyplot as pl
	import os
	import sys
	sys.path.append('/Users/lmatra/radmc-3d/version_0.41/python')
	import radmc3dPy
	import astropy.io.fits as pf
	import problem_setup_12CO_powlaw
	import uuid
	from scipy.ndimage.filters import gaussian_filter1d
	c = 299792458.0	   #m/s
	
	### THIS GOES OUTSIDE MODEL RUN
	# Write dummy wavelength file
	lambda1 = 0.1e0
	lambda2 = 7.0e0
	lambda3 = 25.e0
	lambda4 = 1.0e4
	n12     = 20
	n23     = 100
	n34     = 30
	lam12   = lambda1 * (lambda2/lambda1)**(np.arange(n12, dtype=float)/(1.e0*n12))
	lam23   = lambda2 * (lambda3/lambda2)**(np.arange(n23, dtype=float)/(1.e0*n23))
	lam34   = lambda3 * (lambda4/lambda3)**(np.arange(n34, dtype=float)/(1.e0*(n34-1.e0)))
	lambd = np.concatenate((lam12,lam23,lam34))
	nlam    = len(lambd)
	wfile = open('wavelength_micron.inp', 'w')
	wfile.write('%d\n'%nlam)
	for ilam in range(nlam): wfile.write('%.9e\n'%lambd[ilam])
	wfile.close()
	
	# Set up frequency channels for RADMC run
	#Here if you want to set it up yourself
	deltafreqwanted=244.153e3
	nchanwanted=150
	linefreq=230.538e9#329.3305525e9
	freqc=linefreq-nchanwanted/2e0*deltafreqwanted+deltafreqwanted*np.arange(0,nchanwanted)+deltafreqwanted/2.0
	#Here if you just take it from the data
	#vstarlsrk=2.0
	#freqc=f-(-vstarlsrk/(c/1e3)*linefreq)
	#nchanwanted=f.shape[0]
	#deltafreqwanted=np.abs(f[1]-f[0])
	
	#Oversample spectrally if needed
	#overspecres=4e0 #Keep even!
	#freqi=np.append(freqc-deltafreqwanted/2.0, freqc[freqc.size-1]+deltafreqwanted/2.0)
	#deltafreqover=deltafreqwanted/overspecres
	#for iu in arange(0,freqi.size-1):
	#		if iu==0: freqo=freqi[iu]+np.arange(0,overspecres)*deltafreqover
	#		if iu!=0: freqo=np.concatenate((freqo,freqi[iu]+np.arange(0,overspecres)*deltafreqover))
	#freqo+=deltafreqover/2.0 											#Center velocity value of new overresolved spectral channels
	
	#Otherwise:
	freqo=freqc
	
	#Then get wavelegths
	if (freqo[1]-freqo[0])>0:
		lamo=c/freqo[::-1]*1e6	#in microns
	else:
		lamo=c/freqo*1e6	#in microns
	# image size parameters
	pxsz=0.0062
	dxy=pxsz*np.pi/180.0/3600.0
	npix=128
	imszarcsec=npix*pxsz
	dist=28.55	 #pc
	sizeau=imszarcsec*dist
	
	
	
	##############
	## READ OR DEFINE MODEL PARAMETERS, SET UP MODEL IN NEW INDEP. FOLDER
	##############
	
	p=[5.7, 9.1, 0.0, temp, 0.0, 14.0, logmass, -0.03, -0.02, 1.71, 2.0, 76.5, 100.0, 0.529] #-5.02
	#-4.72
	mstar=p[9]	#Solar masses
	levintdwn=2 #Should be level number AS IN LAMDA FILE
	levintup=3	#Should be level number AS IN LAMDA FILE
	#If running Gaussian model:
	#rmid=p[0]	#AU
	#sigma=p[1]	#AU
	#or if running pow law model:
	rminpowlaw=p[0]
	rmaxpowlaw=p[1]
	gamma=p[2]
	t100au=p[3]	#K
	gammat=p[4]	#Unitless
	mu=p[5]	#Unitless		
	Mdisk=10**(p[6])	#MEarth
	#other free parameters, to run radmc but not to go in problem_setup:
	incl=p[11]#p[6]	#in degrees
	posang=p[12]#p[7]	#in degrees
	#other free parameters, to shift and reweight image:
	dRA=p[7]
	dDec=p[8]
	velstar=p[10]
	wtfact=p[13]
	
	
	#if star velocity is free parameter
	if (freqo[1]-freqo[0])>0:
		lamo=c/(freqo-(-velstar/(c/1e3)*linefreq))[::-1]*1e6	#in microns
	else:
		lamo=c/(freqo-(-velstar/(c/1e3)*linefreq))*1e6	#in microns
	
	ciao = str(uuid.uuid4())
	os.system('mkdir '+'/Users/lmatra/HD172555/radmctautest/'+ciao)
	os.chdir('/Users/lmatra/HD172555/radmctautest/'+ciao)
	os.system('cp '+'/Users/lmatra/HD172555/radmctautest/molecule_12CO.inp .')
	os.system('cp '+'/Users/lmatra/HD172555/radmctautest/wavelength_micron.inp .')
		
		
	problem_setup_12CO_powlaw.problem_setup([mstar,levintdwn,levintup,rminpowlaw, rmaxpowlaw, gamma,t100au,gammat,mu,Mdisk])
	
	##############
	## RUN RADMC AND READ OUTPUT
	##############

	os.system('/Users/lmatra/radmc-3d/version_0.41/srcnoprint/radmc3d image lambdarange '+str(lamo[0])+' '+str(lamo[lamo.size-1])+' nlam '+str(lamo.size)+' incl '+str(incl)+' posang '+str(-posang)+' sizeau '+str(sizeau)+' npix '+str(npix)+' imageunform')
	# Read
	imag     =     radmc3dPy.image.readImage(binary=True)
	cube	=	np.zeros([npix,npix,freqc.size])
	spectrum=np.zeros(freqc.size)
	
	#If needed, Re-average cube to produce desired spectral axis
	#for iw in arange(0,freqc.size): 
	#	cube[:,:,iw]=np.mean(imag.imageJyppix[:,:,(int(overspecres)*iw):(int(overspecres)*(iw+1))],2)/(dist**2e0)
	#	#cube[:,:,iw]=np.mean(imag.image[:,:,(int(overspecres)*iw):(int(overspecres)*(iw+1))],2)#/(dist**2e0)
	#Otherwise:
	cube=imag.imageJyppix/(dist**2e0)
	#If need tau:
	#cube=imag.image
	
	cube=cube[:,:,::-1]
	
	#If need to compare with observations:
	cube=gaussian_filter1d(cube, sigma=2e0/(2e0*np.sqrt(2e0*np.log(2))))		#Convolve spectra with Hanning smoothing (approx. by Gaussian with FWHM=2 chan)
	#print('NB Convolving spectrum by Hanning smoothing!')
	#Pad model so that image is large enough to sample shortest visibility spacings
	#nxy=256
	#modpad=np.zeros((nxy,nxy,freqc.size))
	#modpad[(nxy/2-32):(nxy/2+32),(nxy/2-32):(nxy/2+32),:]+=cube	
	
	os.chdir('..')	
	os.system('rm -r '+ciao)
	
	#CREATE 1D spectrum:
	for iw in np.arange(0,freqc.size):
		spectrum[iw]=np.sum(cube[:,:,iw])

	
	plot=0
	if plot:
		#PLOT MOMENT-0 IF NEEDED
		pl.figure()
		pl.imshow(np.sum(cube, axis=2), origin='lower')
		#pl.imshow(np.amax(cube, axis=2), origin='lower')
		pl.colorbar()
		#pl.savefig('HD172555_mom0_newmass_169K_h0p022.pdf')
		
		#PLOT 1D SPECTRUM IF NEEDED
		#pl.figure()
		#pl.plot(-(freqo-linefreq)/linefreq*c/1e3, spectrum, drawstyle='steps-mid')
		#pl.savefig('HD172555_1dmodeltausum_newmass_169K_h0p022.pdf')
		#np.save('HD172555_1dspec_169K_h0p022.npy', [-(freqo-linefreq)/linefreq*c/1e3, spectrum])
		#pl.plot(freqo, spectrum, drawstyle='steps-mid')
		print('Integrated model line flux is '+str(np.sum(spectrum*deltafreqwanted/linefreq*c/1e3))) #Jy km/s
	
		print(cube.max())
		pl.show()
		
		
	return -(freqo-linefreq)/linefreq*c/1e3, spectrum