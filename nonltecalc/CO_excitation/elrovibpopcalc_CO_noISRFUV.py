def elrovibpopcalc_COfunc(tgas=None, ncollp=None, jnufilestar='', nlinesjnufilestar=None, jnufilerot=None, datafile='',  linewithntrcoll=None, distpc=None, distau=None, massmearth=None, levupint=None, si=0, fluxJykmsinp=None):

	### Packages
	import matplotlib.pyplot as pl
	import numpy as np
	import sys
	import scipy
	from scipy.interpolate import interp1d
	import calc_electcollcoeff
	reload(calc_electcollcoeff)
	
	### Constants
	from scipy.constants import h,c,k  
	
	### Inputs
	#distau: distance of gas from the star in AU (for fluorescence due to stellar emission calculation)
	#distpc: distance of star from Earth in pc
	#ncollp: density of collisional partners (electrons only for now)
	#jnufilestar: numpy save file containing stellar flux as observed in Jy from Earth which needs continuous wavelength (micron) coverage between 0.088 and 5.5 micron, and that can be called with: wavuvum, fluxuvJy=np.load(jnufilestar)
	#nlinesjnufilestar: number of lines in flux file; not needed if file given is numpy array.
	#jnufilerot: file containing disk mean intensity in Jy/sr at all 261 frequencies listed in file 'freq_rottransCO_GHz.txt'. Needs to have size 261 or routine breaks
	#tgas: kinetic temperature of the gas
	#levupint: upper level of the transition of interest, only needed if massmearth or fluxJykmsinp are specified
	#massmearth: mass of gas in Earth masses - only if you want a flux calculation
	#fluxJykmsinp: integrated line flux of gas in Jy km/s - only if you want the mass calculated. NB: the si flag in the input commands only works when flux is the output and not the input. So always give Jy km/s as input for now.
	#datafile: file containing radiative and collisional transitions with collisional partner of choice in LAMDA format
	#linewithntrcoll: line in LAMDA molecular data file containing the number of collisional transitions considered

	##tgas and ncollp are always needed so put a control straight away
	if (not tgas) or (not ncollp) or (not levupint):
		print 'Error!' 
		print 'Gas kinetic temperature (tgas) or collisional partner density (ncollp) or upper level of interest (levup) were not set!'
		sys.exit()
			

	### Read in line transition data 
	fname=datafile
	#fname='COvibrot_Lambda_nelec2_nvib9_nrot30.dat'
	#fname='COrotonly.dat'
	try:
		rfile = open(fname, 'r')
	except:
		print 'Error!' 
		print 'Excitation file was not found!'
		sys.exit() 
	allfilelines=rfile.readlines() 	
	rfile.close()
	rfile = open(fname, 'r')
	
	
	### Read energy levels with quantum numbers
	dum=rfile.readline()
	species=rfile.readline().split()[0]
	dum=rfile.readline()
	mwt=float(rfile.readline().split()[0])
	dum=rfile.readline()
	n_levels=int(rfile.readline().split()[0])
	dum=rfile.readline()
	
	levid=np.zeros(n_levels, dtype=int)
	energycmmin1=np.zeros(n_levels, dtype=float)
	wt=np.zeros(n_levels, dtype=float)
	qnel=np.zeros(n_levels, dtype=int)
	qnvib=np.zeros(n_levels, dtype=int)
	qnrot=np.zeros(n_levels, dtype=int)
	
	for i in np.arange(n_levels):
		levid[i], energycmmin1[i], wt[i], qnel[i], qnvib[i], qnrot[i] = rfile.readline().split()
		#levid[i], energycmmin1[i], wt[i], qnrot[i] = rfile.readline().split()    #If using rotational transitions only
	
	
	### Read radiative transition data
	dum=rfile.readline()
	ntr=int(rfile.readline().split()[0])	
	dum=rfile.readline()
	
	trid=np.zeros(ntr, dtype=int)
	levup=np.zeros(ntr, dtype=int)
	levdwn=np.zeros(ntr, dtype=int)
	einstA=np.zeros(ntr, dtype=float)
	freqghz=np.zeros(ntr, dtype=float)
	eupink=np.zeros(ntr, dtype=float)
	
	for i in np.arange(ntr):
		trid[i], levup[i], levdwn[i], einstA[i], freqghz[i], eupink[i] = rfile.readline().split()
	rfile.close()
	
	### Read collisional transition data in LAMDA format for collider of interest
	if not linewithntrcoll:
		print 'Error!' 
		print 'Line containing number of collisional transitions for collider of interest was not found!'
		sys.exit() 
	ntrcollelect=int(allfilelines[linewithntrcoll-1])
	neleccolltemps=int(allfilelines[linewithntrcoll+1])
	eleccolltemps=np.zeros(neleccolltemps, dtype=float)
	for i in np.arange(neleccolltemps):
		eleccolltemps[i]=allfilelines[linewithntrcoll+3].split()[i]
	
	trcollelectid=np.zeros(ntrcollelect, dtype=int)
	colllevup=np.zeros(ntrcollelect, dtype=int)
	colllevdwn=np.zeros(ntrcollelect, dtype=int)
	collrates=np.zeros((ntrcollelect,neleccolltemps), dtype=float)
	
	# Commented out below here because tabulated electron-CO collision rates are wrong
	#for t in np.arange(neleccolltemps):
	#	for i in np.arange(ntrcollelect):
	#		if t==0:
	#			trcollelectid[i]=allfilelines[linewithntrcoll+5+i].split()[0]
	#			colllevup[i]=allfilelines[linewithntrcoll+5+i].split()[1]
	#			colllevdwn[i]=allfilelines[linewithntrcoll+5+i].split()[2]
	#		collrates[i,t]=allfilelines[linewithntrcoll+5+i].split()[3+t]
	#tgas=float(tgas)
	#collratesattgas=np.zeros(ntrcollelect, dtype=float)
	#for i in np.arange(ntrcollelect):
	#	interpcollrates=interp1d(eleccolltemps,collrates[i,:],kind='linear')
	#	collratesattgas[i]=interpcollrates(tgas)
	#ncoll=ncollp #in cm^-3
	#np.save(collratesav, 'collratesav.npy')
	#collratesav=collrates[np.intersect1d(np.where(colllevup==11)[0],np.where(colllevdwn==10)[0]),:]
	#np.save('collratesav.npy', collratesav)
	#np.save('tempscollsav.npy', eleccolltemps)
	
	
	# Workaround here to take external electron-CO rate coefficients rather than tabulated (wrong) ones
	#for i in np.arange(ntrcollelect):
	#	for t in np.arange(neleccolltemps):			
	#		if t==0:
	#			trcollelectid[i]=allfilelines[linewithntrcoll+5+i].split()[0]
	#			colllevup[i]=allfilelines[linewithntrcoll+5+i].split()[1]
	#			colllevdwn[i]=allfilelines[linewithntrcoll+5+i].split()[2]
	#		if colllevup[i]-colllevdwn[i]==1:
	#			collrates[i,t]=calc_electcollcoeff.calc_electcollcoeff(eleccolltemps[t],colllevup[i]-1,colllevdwn[i]-1)			
	#tgas=float(tgas)
	#collratesattgas=np.zeros(ntrcollelect, dtype=float)
	#for i in np.arange(ntrcollelect):
	#	interpcollrates=interp1d(eleccolltemps,collrates[i,:],kind='linear')
	#	collratesattgas[i]=interpcollrates(tgas)
	#ncoll=ncollp #in cm^-3
	#np.save(collratesav, 'collratesav.npy')
	#collratesav=collrates[np.intersect1d(np.where(colllevup==11)[0],np.where(colllevdwn==10)[0]),:]
	#np.save('collratesav.npy', collratesav)
	#np.save('tempscollsav.npy', eleccolltemps)

	#Better workaround which calculates rates for temperature provided without need for interpolation
	n=0
	tgas=float(tgas)
	collratesattgas=np.zeros(ntrcollelect, dtype=float)
	for i in np.arange(ntrcollelect):
		trcollelectid[i]=allfilelines[linewithntrcoll+5+i].split()[0]
		colllevup[i]=allfilelines[linewithntrcoll+5+i].split()[1]
		colllevdwn[i]=allfilelines[linewithntrcoll+5+i].split()[2]
		#print colllevup[i]-colllevdwn[i]
		if (colllevup[i]-colllevdwn[i])==1 and (colllevup[i]<31):
			#print colllevup[i]-collevdwn[i]
			n+=1
			collratesattgas[i]=calc_electcollcoeff.calc_electcollcoeff(tgas,colllevup[i]-1,colllevdwn[i]-1)	
	ncoll=ncollp #in cm^-3
	#print n


	### Read in UV flux from ASCII table
	##flname='bPicNormFlux1AU.dat'
	#if jnufilestar:
	#	flname=jnufilestar
	#	#nlines=261549
	#	if not nlinesjnufilestar:
	#		print 'Error!' 
	#		print 'Stellar flux file name was set but the number of lines in this file was not set!'
	#		sys.exit()
	#	nlines=int(nlinesjnufilestar)
	#	try:
	#		rfile = open(flname, 'r')
	#	except:
	#		print 'Error!' 
	#		print 'Stellar flux file name was set but the file was not found!'
	#		sys.exit()
	#	if (not distpc) or (not distau):
	#		print 'Error!' 
	#		print 'Stellar flux file name (containing observed fluxes in Jy) was set but either the distpc or distau were not set!'
	#		sys.exit()
	#	wavuvum=np.zeros(nlines)
	#	fluxuvJy=np.zeros(nlines)
	#	for i in np.arange(nlines):
	#		wavuvum[i], fluxuvJy[i]=rfile.readline().split()
	#	
	#	#wavuvum=wavuv*1e-4
	#	freqjnughz=c/(wavuvum*1e-6)/1e9
	#	fluxuvSI1AU=fluxuvJy*1e-26*((distpc*2.0626*1e5)**2.0)
	#	#(10**logfluxuv)/2.99792458e21*(wavuv**2.0)
	#	fluxuvSIat85au=fluxuvSI1AU/(distau**2.0)
	#	
	#	freqjnughz=freqjnughz[::-1]
	#	fluxuvSIat85au=fluxuvSIat85au[::-1]
	#	
	#	### Check that flux makes sense as seen from Earth in Jy
	#	#pl.figure()
	#	#pl.plot(wavuvum, fluxuvSI1AU*1e26/((distpc*2.0626*1e5)**2.0))
	#	#pl.yscale('log')
	#	#pl.xscale('log')
	#
	### Read in stellar flux from numpy array
	if jnufilestar:
		wavuvum, fluxuvJy=np.load(jnufilestar)
		freqjnughz=c/(wavuvum*1e-6)/1e9
		fluxuvSI1AU=fluxuvJy*1e-26*((distpc*2.0626*1e5)**2.0)#/1000.0
		fluxuvSIat85au=fluxuvSI1AU/(distau**2.0)
		freqjnughz=freqjnughz[::-1]
		fluxuvSIat85au=fluxuvSIat85au[::-1]
	### Read in disk mean intensity for rotational levels from numpy array
	if jnufilerot:
		freqjnughzrotonly=freqghz[:261]
		fluxrotonlyJysr=np.load(jnufilerot)
		if (fluxrotonlyJysr.size != freqjnughzrotonly.size):
			print 'Error!' 
			print 'Number of input mean intensities for rotational transitions (far-IR to mm) does not match the number of frequencies for these transitions provided!'
			sys.exit()
					

	
	#Set up all matrices in preparation of solving statistical equilibrium
	A=np.zeros((n_levels, n_levels))
	BJstim=np.zeros((n_levels, n_levels))
	BJabs=np.zeros((n_levels, n_levels))
	K=np.zeros((n_levels, n_levels))
	Jstim=np.zeros((n_levels, n_levels))
	Jabs=np.zeros((n_levels, n_levels))
	T_cmb		=	float(2.72548)			#temperature of CMB
	jnu=2.*h*((freqghz*1e9)**3.0)/((c**2.0)*(np.exp(h*freqghz*1e9/(k*T_cmb))-1.))	#Bnu(CMB)
	if jnufilestar:
		fint=interp1d(freqjnughz, fluxuvSIat85au, kind='linear')
		print 'Make sure radiation field covers well the UV out to '+str(c/freqghz.max()/1e9*1e6)+' micron - right now it goes out to '+str(c/freqjnughz.max()/1e9*1e6)+' micron'
		jnu[np.intersect1d(np.where(freqghz>freqjnughz.min())[0],np.where(freqghz<freqjnughz.max())[0])]+=1.0/4.0/np.pi*fint(freqghz[np.intersect1d(np.where(freqghz>freqjnughz.min())[0],np.where(freqghz<freqjnughz.max())[0])]) #Note division by 4pi here as we want mean intensity rather than flux
	if jnufilerot:
		jnu[:261]+=fluxrotonlyJysr*1e-26	

	#pl.figure()
	#pl.plot(c/freqghz/1e9*1e6, jnu, '+')
	#pl.xscale('log')
	#pl.yscale('log')
	#print np.max(colllevup)
	### Prepare matrices for radiative transitions
	for i in np.arange(ntr):
		if (einstA[i] !=0.0):
			A[np.where(levid==levdwn[i])[0], np.where(levid==levup[i])[0]]=einstA[i]	#Spontaneous emission into the level given by the row
			BJstim[np.where(levid==levdwn[i])[0], np.where(levid==levup[i])[0]]=einstA[i]*(c**2.0)/(2.0*h*((freqghz[i]*1e9)**3.0))*jnu[i]	#Stimulated emission into the level given by the row
			BJabs[np.where(levid==levup[i])[0], np.where(levid==levdwn[i])[0]]=wt[np.where(levid==levup[i])[0]]/wt[np.where(levid==levdwn[i])[0]]*einstA[i]*(c**2.0)/(2.0*h*((freqghz[i]*1e9)**3.0))*jnu[i]	#Absorption into the level given by the row (different from one above)
	
	#print 'ciao', ntr
	### Prepare matrices for collisional transitions
	for i in np.arange(ntrcollelect):
		if (collratesattgas[i] != 0.0):
			K[np.where(levid==colllevdwn[i])[0], np.where(levid==colllevup[i])[0]]=collratesattgas[i]*ncoll		#Collisional stimulated emission into the level given by the row
			K[np.where(levid==colllevup[i])[0], np.where(levid==colllevdwn[i])[0]]=collratesattgas[i]*ncoll*wt[np.where(levid==colllevup[i])[0]]/wt[np.where(levid==colllevdwn[i])[0]]*np.exp(-float(1.43877736)*np.abs(energycmmin1[np.where(levid==colllevup[i])[0]]-energycmmin1[np.where(levid==colllevdwn[i])[0]])/tgas)		#Collisional absorption into the level given by the row (different from one above)
	
	#pl.figure()
	#pl.imshow(np.log10(np.abs(K)), origin='lower')	

	#Now take into account transitions out of level (along diagonal)
	for i in np.arange(n_levels):
		if i!=0: A[i,i]=-np.sum(A[:,i])
		if i!=0: BJstim[i,i]=-np.sum(BJstim[:,i])
		if (i!=(n_levels-1)): BJabs[i,i]=-np.sum(BJabs[:,i])
		K[i,i]=-np.sum(K[:,i])
	
	S=A+BJstim+BJabs+K
	S[0,:]=1.0 #NB: for some reason if I set -1 instead of 0 as the normalising equation, it doesnt work.
	rhs=np.zeros((n_levels))
	rhs[0]=1.0
	#return S, rhs
	#print 'ciao'
	fracpoplev=scipy.linalg.lu_solve(scipy.linalg.lu_factor(S), np.transpose(rhs))
	#fracpoplev=scipy.linalg.solve(S, np.transpose(rhs))
	#return
	#print 'ciao'
	
	### Calculate flux from input mass (not needed, can just return populations)
	#massmearth=1e-5
	stimem=0
	if levupint=='all':
		fluxwm2=np.zeros(ntr)
		fluxjykms=np.zeros(ntr)
		fracpoplevup=np.zeros(ntr)
		mass=np.zeros(ntr)
		if massmearth:
			for i in np.arange(ntr):
				fracpoplevup[i]=fracpoplev[np.where(levid==levup[i])[0]]
				if stimem: fluxwm2[i]=massmearth*5.9736e24*h*freqghz[i]*1e9*einstA[i]*(1.0+(c**2.0)/(2.0*h*((freqghz[i]*1e9)**3.0))*jnu[i])*fracpoplevup[i]/(4.0*np.pi*1.67e-27*mwt*((distpc*3.0857e16)**2.0))
				else: fluxwm2[i]=massmearth*5.9736e24*h*freqghz[i]*1e9*einstA[i]*(1.0)*fracpoplevup[i]/(4.0*np.pi*1.67e-27*mwt*((distpc*3.0857e16)**2.0))
			if si:
				#pl.figure()
				#pl.plot(c/freqghz[np.intersect1d(np.where(levup>0)[0],np.where(levup<5)[0])]/1e9*1e6, fluxwm2[np.intersect1d(np.where(levup>0)[0],np.where(levup<5)[0])], 'o')
				return fluxwm2
			else:
				fluxjykms=fluxwm2/freqghz/1e9*c/1e3*1e26
				return fluxjykms
		else: 
			if fluxJykmsinp:
				for i in np.arange(ntr):
					fracpoplevup[i]=fracpoplev[np.where(levid==levup[i])[0]]
					if stimem: mass[i]=fluxJykmsinp*freqghz*1e9/c*1e3/1e26/(5.9736e24*h*freqghz[i]*1e9*einstA[i]*(1.0+(c**2.0)/(2.0*h*((freqghz[i]*1e9)**3.0))*jnu[i])*fracpoplevup[i]/(4.0*np.pi*1.67e-27*mwt*((distpc*3.0857e16)**2.0)))
					else: mass[i]=fluxJykmsinp*freqghz*1e9/c*1e3/1e26/(5.9736e24*h*freqghz[i]*1e9*einstA[i]*(1.0)*fracpoplevup[i]/(4.0*np.pi*1.67e-27*mwt*((distpc*3.0857e16)**2.0)))
				return mass
			else:
				return fracpoplev				
	else:
		if massmearth:
			#print (c**2.0)/(2.0*h*((freqghz[np.where(levup==(levupint+1))[0]]*1e9)**3.0))*jnu[np.where(levup==(levupint+1))[0]]
			if stimem: fluxwm2=massmearth*5.9736e24*h*freqghz[np.where(levup==(levupint+1))[0]]*1e9*einstA[np.where(levup==(levupint+1))[0]]*(1.0+(c**2.0)/(2.0*h*((freqghz[np.where(levup==(levupint+1))[0]]*1e9)**3.0))*jnu[np.where(levup==(levupint+1))[0]])*fracpoplev[levupint]/(4.0*np.pi*1.67e-27*mwt*((distpc*3.0857e16)**2.0))
			else: fluxwm2=massmearth*5.9736e24*h*freqghz[np.where(levup==(levupint+1))[0]]*1e9*einstA[np.where(levup==(levupint+1))[0]]*(1.0)*fracpoplev[levupint]/(4.0*np.pi*1.67e-27*mwt*((distpc*3.0857e16)**2.0))
			if si:
				return fluxwm2
			else:
				fluxjykms=fluxwm2/freqghz[np.where(levup==(levupint+1))[0]]/1e9*c/1e3*1e26
				return fluxjykms
		else:
			if fluxJykmsinp:
				if stimem: mass=fluxJykmsinp*freqghz[np.where(levup==(levupint+1))[0]]*1e9/c*1e3/1e26/(5.9736e24*h*freqghz[np.where(levup==(levupint+1))[0]]*1e9*einstA[np.where(levup==(levupint+1))[0]]*(1.0+(c**2.0)/(2.0*h*((freqghz[np.where(levup==(levupint+1))[0]]*1e9)**3.0))*jnu[np.where(levup==(levupint+1))[0]])*fracpoplev[levupint]/(4.0*np.pi*1.67e-27*mwt*((distpc*3.0857e16)**2.0)))
				else: mass=fluxJykmsinp*freqghz[np.where(levup==(levupint+1))[0]]*1e9/c*1e3/1e26/(5.9736e24*h*freqghz[np.where(levup==(levupint+1))[0]]*1e9*einstA[np.where(levup==(levupint+1))[0]]*(1.0)*fracpoplev[levupint]/(4.0*np.pi*1.67e-27*mwt*((distpc*3.0857e16)**2.0)))		
				#print 'ciao', mass.shape
				return mass
			else:
				return fracpoplev[levupint]
