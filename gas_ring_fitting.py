import math
import astropy.units as u
import astropy.constants as const
import numpy as np

from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter1d

class gasringfit():
    def __init__(self,
                 stellar_mass):
        """
        Inputs:
            -   stellar_mass:   The value of the stellar mass in terms of
                                solar masses
            -   mask_size_corr_factor:  mask size is the correction factor times
                                        the HWHM
        """
        self.pi = math.pi
        self.lightspeed = (const.c.value*const.c.unit).to(u.m/u.s)
        self.gravconst = (const.G.value *const.G.unit)
        self.solar_mass = (const.M_sun.value*const.M_sun.unit).to(u.kg)
        self.stellar_mass = stellar_mass#
        self.mask_size_corr_factor = 2.12
        self.omega = 0

    def return_alma_beamsize(self,filename):
        cube,header_cube = fits.getdata(filename,0,header=True)
        self.bmaj = header_cube['BMAJ']
        self.bmin = header_cube['BMIN']


    def return_alma_data(self,filename,beaminfo_filename):
        self.return_alma_beamsize(beaminfo_filename)
        cube,header_cube = fits.getdata(filename,0,header=True)
        cube = cube[0,:,:,:]
        cube[np.isnan(cube)] = 0  #changes all nan values to 0

        #gets the pixel to arcsecond conversion factor
        self.pixtoarcsecfactor = np.abs(header_cube['CDELT1']*3600.)
        #gets the array of frequencies observed
        self.data_freqarray = header_cube['CRVAL3']+\
                    header_cube['CDELT3']*np.arange(0,header_cube['NAXIS3'])
        #gets the channel width
        self.chanwidth = abs(header_cube['CDELT3'])
        #gets the frequency of the CO line/the line in question
        self.line_freq = header_cube['RESTFRQ']*u.Hz
        #sets a mask radius for the image. Takes the average radius of the beam,
        #converts that to a pixel value, and scales it upwards. BMAJ and BMIN are
        #the FWHM of the major and minor axes of the beam.
        mask_size = int((((self.bmaj+self.bmin)/4)*3600/self.pixtoarcsecfactor)*self.mask_size_corr_factor)
        beamarea = (self.pi/(4*math.log(2)))*(self.bmaj*self.bmin)

        #creates a mask for the ALMA data
        mask = np.zeros([len(cube[0,:,:]),len(cube[0,:,:])])
        midpoint = int(len(mask[0,:])/2) - 1
        for x_index in range(0,len(mask[0,:])):
            for y_index in range(0,len(mask[0,:])):
                dist = math.sqrt((midpoint - x_index)**2 + (midpoint - y_index)**2 )
                #print(dist)
                if dist<=mask_size:
                    mask[x_index,y_index] = 1

        #This section sums the intensities at a given frequency. It then rescales the value
        #to Jy/frequency/pixel, rather than Jy/frequency/beam size
        pixelarea = np.abs(header_cube['CDELT1'])**2
        self.pixelsperbeam = beamarea/pixelarea

        freq_is = np.zeros(len(cube[:,0,0]))
        for index in range(0,len(cube[:,0,0])):
            temp_img = cube[index,:,:]
            freq_is[index] = temp_img[mask>0].sum()
            freq_is[index] = freq_is[index]/self.pixelsperbeam
        self.data_intensities = freq_is
        self.data_integflux = sum(self.data_intensities*self.chanwidth)


    def veltofreq(self,vel,freq_o,stellar_v):
        """
        Inputs:
            -   vel is the velocity to convert in m/s
            -   freq_o is the central frequency of the line
            -   stellar_v is the stellar velocity with respect to the viewer
        Outputs:
            -   freq, the Doppler shifted frequency corresponding to the velocity
                in question
        """
        freq = freq_o * ((self.lightspeed.value - stellar_v)/(self.lightspeed.value + vel))
        return(freq)


    def freqtovel(self,freq,freq_o,stellar_v):
        """
        Inputs:
            -   freq is the frequency to convert
            -   freq_o is the central frequency of the line
            -   stellar_v is the stellar velocity with respect to the viewer
        Outputs:
            -   vel, the Doppler shifted velocity corresponding to the frequency in
                question in m/s
        """
        vel = ((freq_o/freq)*(self.lightspeed.value - stellar_v)) - self.lightspeed.value
        return(vel)


    def normalize(self,array):
        """
        Normalizes an array between the values of 0 and 1
        """
        minimum = min(array)
        maximum = max(array)
        normarray = np.zeros(len(array))
        for index in range(0,len(array)):
            normarray[index] = (array[index] - minimum)/(maximum - minimum)
        return(normarray)


    def find_nearest(self,array,value):
        """
        Inputs:
            -   array is the array you're searching for a value in
            -   value is the value you're searching for
        Returns:
            -   idx, the index of the array entry closest to the value
        """
        temp = 1e20
        for i in range(0,len(array)):
            if  np.abs(array[i]-value) <= temp:
                temp = np.abs(array[i]-value)
                idx = i
        return(idx)


    def keplerianring(self, inc, rmid, rwid, stellar_v, intflux):
        """
        Returns binned intensity values for a ring of material orbiting  in a
        keplerian orbit around a star. Normalizes the intensity by the integrated
        flux across the bands in consideration. Takes as inputs:
            -   inc: inclination in radians
            -   rmid: distance of the center point of the ring in AU
            -   rwid: width of the ring
            -   stellar_v: the stellar velocity of the star wrt to observer in m/s
            -   intflux: the integrated flux across the binned intensity values
        """

        #rmid = rmid.to(u.AU).value
        #rwid = rwid.to(u.AU).value
        #stellar_v = stellar_v.to(u.m/u.s).value
        numradelements = 500
        numphielements = 500

        rinner = rmid - rwid/2
        router = rmid + rwid/2

        radarray = np.linspace(rinner,router,numradelements)
        phiarray = np.linspace(0,2*self.pi,numphielements)

        modelgasloc = np.zeros((numradelements,numphielements))
        for radindex in range(0,numradelements):
            rad = radarray[radindex]
            if rad>=rinner and rad<=router:
                modelgasloc[radindex,:] = 1

        normv = math.sqrt(self.gravconst.value * self.solar_mass.value / ((const.au).value))
        kepv = lambda r: normv * math.sqrt(self.stellar_mass/r)

        sky_zvel = np.zeros((numradelements,numphielements))
        conv_freqs = np.zeros((numradelements,numphielements))

        #the following calculation makes the assumption that the z-velocity
        #in the orbital plane is zero. Furthermore, it eliminates matrix
        #multiplication and calculates only sky_zvelocity
        for radindex in range(0,numradelements):
            for phiindex in range(0,numphielements):
                kepvel = kepv(radarray[radindex])
                s_zv = (kepvel * math.sin(phiarray[phiindex]))*\
                        math.sin(self.omega)*math.sin(inc) + \
                       (kepvel * math.cos(phiarray[phiindex]))*\
                        math.cos(self.omega) * math.sin(inc)
                sky_zvel[radindex,phiindex] = s_zv
                conv_freqs[radindex,phiindex] = self.veltofreq(s_zv,self.line_freq.value,stellar_v)

        binedges = np.zeros(len(self.data_freqarray)+1)
        sep = self.chanwidth
        for index in range(0,len(binedges) - 1):
            binedges[index] = self.data_freqarray[index] - 0.5*sep
        binedges[-1] = self.data_freqarray[-1]+0.5*sep
        binedges = binedges#[::-1]
        binnedmodelcounts,binnedmodeledges = np.histogram(conv_freqs[modelgasloc>0],
            bins = binedges,
            weights = 1/np.reshape(np.repeat(radarray,len(phiarray)),
                                   conv_freqs.shape)[modelgasloc>0])
        normbinnedmodel = self.normalize(binnedmodelcounts)#[::-1]
        scalval = sum(normbinnedmodel*sep)#chanwidth)
        normbinnedmodel = normbinnedmodel*intflux/scalval
        convolvedmodel = gaussian_filter1d(normbinnedmodel,
                                           2/(2*math.sqrt(2*math.log(2))))
        return(convolvedmodel)
