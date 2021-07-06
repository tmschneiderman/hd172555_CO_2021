def problem_setup(theta):

	import numpy as np
	import matplotlib.pyplot as pl

	mstar=theta[0]	#Solar masses
	levintdwn=theta[1] #Should be level number AS IN LAMDA FILE
	levintup=theta[2]	#Should be level number AS IN LAMDA FILE
	rminpowlaw=theta[3]	#AU
	rmaxpowlaw=theta[4]	#AU
	gamma=theta[5]
	t100au=theta[6]	#K
	gammat=theta[7]	#Unitless
	mu=theta[8]	#Unitless
	Mdisk=theta[9]	#MEarth
	atomicmassisotope=28.0#30.0 #30 for C18O

	#
	# Some natural constants
	#
	AU  = 1.49598e13     # Astronomical Unit       [cm]
	ME  = 5.97219e27 # Earth mass			[g]	

	###################################################################################
	############################# MAKE POSITIONAL GRID #######################
	####################################################################################
	rmin       =     5.0*AU
	rmax        =     10.0*AU
	nr		= 20.0	
	nth		= 38.0 
	nphi		= 3.0
	
	vertcutoff	= 5e0
	be=1.0
	H1AU=0.1#0.0076#0.022#0.1
	
	logscalconst	=	((rmax/rmin)**(1e0/(1e0*nr)))-1
	totalthetacutoff=2e0*abs(np.arcsin((np.min([0.95,vertcutoff*H1AU*np.asarray([(rmin/AU)**(be-1.0), (rmax/AU)**(be-1.0)]).max()]))))

	thmin				=	np.pi/2e0-totalthetacutoff/2e0
	thmax				=	np.pi/2e0+totalthetacutoff/2e0

	closeto0	=	1e-2
	logscalconstth	=	((totalthetacutoff/2e0/np.radians(closeto0))**(1e0/(1e0*nth/2e0)))-1

	phimin	=	0.0
	phimax	=	2e0*np.pi
	sizephi   = 	(phimax-phimin)/(nphi)

	ri       = rmin*(logscalconst+1)**(np.arange(nr+1, dtype=float))	
	#ri       = (np.arange(nr+1, dtype=float))/nr*(rmax-rmin)+rmin
	thi      = np.radians(closeto0)*(logscalconstth+1)**np.arange(nth/2.0+1, dtype=float)
	thi 	 = np.concatenate([-thi[::-1],[0.0],thi])+np.pi/2e0
	ti=thi

	pi       = phimin + sizephi*np.arange(nphi+1, dtype=float)

	nr	 = len(ri)-1
	dr	 = (np.subtract(ri[1:],ri[:-1]))
	rc       = np.asarray(ri[:-1]+dr/2e0)
	
	nth	  = len(thi)-1
	dt	  = (np.subtract(thi[1:],thi[:-1]))
	tc       = np.asarray(thi[:-1]+dt/2e0)
	#print tc
	#print dt

	nphi = len(pi)-1
	dp = (np.subtract(pi[1:],pi[:-1]))
	pc = np.asarray(pi[:-1]+dp/2e0)
	
	wfile=open('./amr_grid.inp', 'w')
	wfile.write('%d\n'%1)                      # Format number
	wfile.write('%d\n'%0)                    # AMR self.style (0=regular self. NO AMR)
	wfile.write('%d\n'%150)                  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical
	wfile.write('%d\n'%0)                    # Gridinfo
	wfile.write('%d %d %d \n'%(1,1,1))       # Which dimension is active
	wfile.write('%d %d %d \n'%(rc.size, tc.size, pc.size))    # Grid size (x,y,z or r,phi,theta, or r,phi,z)
	for i in range(rc.size+1): wfile.write('%.9e\n'%ri[i])			# X coordinates (cell walls)
	for i in range(tc.size+1): wfile.write('%.9e\n'%ti[i])		# Y coordinates (cell walls)
	for i in range(pc.size+1): wfile.write('%.9e\n'%pi[i])			# Z coordinates (cell walls)
	wfile.close()							

	
	###################################################################################
	############################# MAKE TEMPERATURE GRID #######################
	####################################################################################

	rr=np.reshape(np.repeat(rc, tc.size*pc.size), (rc.size, tc.size, pc.size))
	drdr=np.reshape(np.repeat(dr, tc.size*pc.size), (rc.size, tc.size, pc.size))
	dpdp=np.reshape(np.repeat(dp, rc.size*tc.size), (rc.size, tc.size, pc.size))
	tt=np.reshape(np.repeat(tc, rc.size*pc.size), (rc.size, tc.size, pc.size))
	dtdt=np.reshape(np.repeat(dt, rc.size*pc.size), (rc.size, tc.size, pc.size))
	Tkgrid=t100au*(rr/AU/100.0)**gammat
	#pl.figure()
	#pl.plot(rr[:,0,0]/AU, Tkgrid[:,0,0])
	#print(Tkgrid)
	
	wfile = open('./gas_temperature.inp', 'w')
	wfile.write('%d\n'%1)                   		# Format number
	wfile.write('%d\n'%(rc.size*tc.size*pc.size)) 		# Nr of cells
	for ip in range(pc.size):
	   for ith in range(tc.size):
	      for ir in range(rc.size):
	         wfile.write('%.9e\n'%(Tkgrid[ir,ith,ip]))
	wfile.close()	

	###################################################################################
	############################# MAKE NUMBER DENSITY GRID #######################
	####################################################################################

	M	=	Mdisk*ME #IN GRAMS!!
	
	ncogrid	= np.zeros([rc.size, tc.size, pc.size])
	rhodtp		= np.zeros([rc.size, tc.size])	

	#sd = np.exp(-((rc/AU-rmid/AU)**2.0/(2.0*(std/AU)**2.0))) # for Gaussian
	sd = (rc/(rminpowlaw*AU))**gamma #for pow law
	sd[rc>rmaxpowlaw*AU]=0.0 #for pow law
	sd[rc<rminpowlaw*AU]=0.0 #for pow law

	surfmassdens	= np.reshape(np.repeat(sd*M/np.sum(sd*2.0*np.pi*rc*dr),pc.size), (rc.size,pc.size))
	#print np.sum(surfmassdens[:,0]*2.0*np.pi*rc*dr/ME)
	#pl.figure()
	#pl.plot(rc/AU, surfmassdens)
	
	#Now calculate scale height in SI then convert to cm given temperature and mu at a given radius (and mstar)
    #Assume hydrostatic equilibrium and calculate from temperature, radius, molecular mass, and stellar mass.
	H = np.sqrt(1.38e-23*Tkgrid[:,0,0]*((rc/1e2)**3.0)/mu/1.67e-27/6.67e-11/mstar/1.98892e30)*1e2

	
	for i in range(rc.size):
		for j in range(tc.size):
			
			#HERE USE VERTICAL GAUSSIAN:
			rhodtp[i,j]    = np.exp(-(rc[i]*np.sin(((np.pi/2e0)-tc[j]))/(np.sqrt(2.0)*H[i]))**2e0)/(np.sqrt(2.0*np.pi)*H[i]) 
		
		rhodtp[i,:] = 	rhodtp[i,:]/np.sum(rhodtp[i,:]*rc[i]*(dt))
		#print np.sum(rhodtp[i,:]*rc[i]*(dt))
		for k in range(pc.size): ncogrid[i,:,k] = rhodtp[i,:]*surfmassdens[i,k]
	#print np.sin(tc)
	#print 'mass in the model is '+str(np.sum(ncogrid*(rr**2)*np.sin(tt)*dtdt*drdr*dpdp)/ME)+' Earth masses'
	#Turn to number density and not volume density:
	ncogrid/=(atomicmassisotope*1.67e-24) #e-24 perche si parla di grammi!!
		
		
	wfile=open('./numberdens_12CO.inp', 'w')
	wfile.write('%d\n'%1)					# Format
	wfile.write('%d\n'%(rc.size*tc.size*pc.size)) 		# Nr of cells
	for ip in range(pc.size):
	   for ith in range(tc.size):
	      for ir in range(rc.size):
	      	wfile.write('%.9e\n'%(ncogrid[ir,ith,ip]))
	wfile.close()

	
	
	
	
	###############################################################
	####### CALCULATE GAS KEPLERIAN VELOCITIES GIVEN Mstar ########
	###############################################################
	
	#Calculate Keplerian velocity in cm/s along phi direction (counterclockwise in orbital frame)
	
	vkep_phi = +29.8e5*np.sqrt(mstar)/np.sqrt(rc/AU)	# Keplerian velocities at each cell centre rc in cm/s, IN PHI DIRECTION.
	#print vkep_phi.shape
	vkepgrid_phi = np.resize(vkep_phi, (tc.size,pc.size,vkep_phi.size))
	vkepgrid_phi = np.swapaxes(vkepgrid_phi,1,2)
	vkepgrid_phi = np.swapaxes(vkepgrid_phi,0,1)
	
	#velocity is zero in the r direction (outwards) and in the theta direction (downwards)
	vkepgrid_r   = np.zeros(vkepgrid_phi.shape)	
	vkepgrid_theta  = np.zeros(vkepgrid_phi.shape)
	
	wfile = open('./gas_velocity.inp', 'w')
	wfile.write('%d\n'%1)                   		# Format number
	wfile.write('%d\n'%(rc.size*tc.size*pc.size)) 		# Nr ofnewpi[i] cells
	for ip in range(pc.size):
	   for ith in range(tc.size):
	      for ir in range(rc.size):
	         wfile.write('%.9e %.9e %.9e\n'%(vkepgrid_r[ir,ith,ip],vkepgrid_theta[ir,ith,ip],vkepgrid_phi[ir,ith,ip]))
	wfile.close()		
	
	
	###############################################################
	#### PRINT RADMC-3D FILES THAT ARE INDEP. OF ANY PARAMETER ####
	###############################################################
	
	#
	# Write the lines.inp file
	#
	wfile=open('./lines.inp', 'w')
	wfile.write('%d\n'%2)                   		# Format number
	wfile.write('%d\n'%1)					# Nr of species to be modelled
	wfile.write('%s %s %d %d %d\n'%('12CO','leiden',1,2,0))	# Watch out: the final 0 means LTE... but here we are giving it the populations, so RADMC3D shouldn't ask - though it works. 1 means we are providing it the subset of levels we are interested in, and 2 is the number of levels we are interested in
	wfile.write('%d %d\n'%(levintdwn,levintup))		# Levels we are interested in
	wfile.close()
	#
	# Write the radmc3d.inp file
	#lines_autosubset = 0
	wfile = open('./radmc3d.inp', 'w')
	#wfile.write('%s %d\n'%('nphot'+' =',nphot))
	#wfile.write('%s %s\n'%('scattering_mode_max'+' =','0'))
	#wfile.write('%s %s\n'%('iranfreqmode '+' =','1'))
	#wfile.write('%s %s\n'%('istar_sphere '+' = ','0' ))
	wfile.write('%s %s\n'%('lines_mode '+' = ','1' ))
	wfile.write('%s %s\n'%('incl_dust '+' = ','0' ))
	wfile.write('%s %s\n'%('incl_lines '+' = ','1' ))
	wfile.write('%s %s\n'%('lines_autosubset '+' = ','0' ))	# IMPORTANT: This tells it that we are selecting levels of interest for the given species manually in the lines.inp file
	wfile.close()


	

	return 0.0
