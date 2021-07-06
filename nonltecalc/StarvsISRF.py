import matplotlib.pyplot as pl

#Plotting
font={'family':'Times New Roman', 'size':23}	
rc('font', **font)
rc('axes', linewidth=2)
rc('xtick.major', width=2, size=8)
rc('xtick.minor', width=2, size=4)
rc('ytick.major', width=2, size=8)
rc('ytick.minor', width=2, size=4)
rc('xtick', labelsize=20)
rc('ytick', labelsize=20)
matplotlib.rcParams.update({'figure.autolayout': True})


starfluxfilename='./HD172555_starfluxJyfromEarth.npy'
wavuvum, fluxuvJy=np.load(starfluxfilename)
distpc=28.3
fluxuvSI1AU=fluxuvJy*1e-26*((distpc*2.0626*1e5)**2.0)#/1000.0
distau=7.5
fluxuvSIat85au=fluxuvSI1AU/(distau**2.0)
fluxuvJyat85au=fluxuvSIat85au*1e26

#Read interstellar Draine field flux file
rfile = open('/sma/SMAusers/lmatra/CO_excitation/ISRF.dat', 'r')
dum=rfile.readline()
dum=rfile.readline()
dum=rfile.readline()
dum=rfile.readline()
nlinesisrf=1909
wavisrfnm=np.zeros(nlinesisrf)
fluxisrfweird=np.zeros(nlinesisrf)
#Read wavelengths in nm and fluxes in photons cm^-2 s^-1 nm^-1
for i in np.arange(nlinesisrf):
	wavisrfnm[i], fluxisrfweird[i] = rfile.readline().split()
wavisrfang=wavisrfnm*10.0
fluxisrf=fluxisrfweird/10.0 # Now in photons cm^-2 s^-1 A^-1
#print wavisrfang.min(), wavisrfang.max()
#print wavcrsectang.min(), wavcrsectang.max()
fluxisrfjy=fluxisrf*6.63e-4*wavisrfang

#Read in CO uv line wavelengths for overplotting
datafile='/sma/SMAusers/lmatra/CO_excitation/COvibrot_Lambda_nelec2_nvib9_nrot30.dat'
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




#Now plot it all
pl.figure()
#Choose random 100 frequency so as to not plot them all
randfreq=np.random.choice(freqghz, 1000)
for i in np.arange(1000):
	pl.axvline(2.9979e8/(randfreq[i]*1e9)*1e6, alpha=0.005, color='red')
pl.axvline(920e-4, alpha=0.3, color='black')
pl.axvline(1080e-4, alpha=0.3, color='black')
pl.plot(wavuvum, fluxuvJyat85au, label='Star')
pl.plot(wavisrfang*1e-4, fluxisrfjy, label='ISRF')
pl.xscale('log')
pl.yscale('log')
pl.ylim(1e2,5e13)
pl.xlim(4e-2,1e1)
pl.xlabel('Wavelength (micron)')
pl.ylabel('Flux @ 7.5 au (Jy)')
ax=pl.gca()
pl.text(0.68,0.88, 'HD172555', transform=ax.transAxes, fontsize=22)
pl.legend(frameon=False)
pl.savefig('StarvsISRFforphotodissandfluo.pdf')
pl.close()
import os
os.system('evince StarvsISRFforphotodissandfluo.pdf &')
