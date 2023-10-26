import numpy as np
import mpmath as mp
from LinkedAttribute import *

class ComputeTurbulenceStrength(object):
    """Definition of the attributes that need to be taken from another instance so that they are automatically updated if the instance has them changed.
    For more info read comments in LinkedAtrribute class"""
    c2n = LinkedAttribute('refractiveIndexObject')

    def __init__(self, refractiveIndexObject,  wavelength = 1550e-9, elevation=90., geometry = 'uplink', ALT_altitude=4e3, SLT_altitude=600e3, W_0 = None, F_0 = None ):
        # ACCESSED PARAMETERS
        self.refractiveIndexObject = refractiveIndexObject  # ComputeRefractiveIndexStructure instance
        self.c2n = refractiveIndexObject.c2n                # function,  function of c2n (in m**(-2/3)) vs altitude of height ,accessed from ComputeRefractiveIndexStructure instance

        # INPUT PARAMETERS
        self.wavelength = wavelength                        # float,     beam wavelength (in m), by default 1550e-9 m
        self.elevation = elevation                          # float,     elevation angle (in degrees), by default 90 degrees
        self.geometry = geometry                            # string,    link geometry, by default 'downlink'
        self.ALT_altitude = ALT_altitude                    # float,     ALT altitude (in m), by default 4e3 m
        self.SLT_altitude = SLT_altitude                    # float,     SLT altitude (in m), by default 600e3 m
        self.W_0 = W_0                                      # float,     beam waist at transmitter (in m)
        self.F_0 = F_0                                      # float,     phase front radius of curvature (in m)

        # COMPUTED PARAMETERS
        self.zenith_angle = 90 - elevation                  # float,     zenith angle (in degrees)
        self.zenith_angle_rad = self.compute_zenith_angle_rad()#float,   zenith angle (in rad)
        self.k = (2 * np.pi) / self.wavelength              # float,     beam wavenumber (in m**-1)
        self.R = None                                       # float,     distance of propagation ALT -SLT (in m)
        self.W = None                                       # float,     diffractive beam radius at the receiver (in m)
        self.cap_lambda_0 = None                            # float,     adimensional beam parameter at the transmitter, by default None
        self.cap_lambda = None                              # float,     adimensional beam parameter at the receiver, by default None
        self.cap_theta_0 = None                             # float,     adimensional beam parameter at the transmitter, by default None
        self.cap_theta_0_bar = None                         # float,     adimensional beam parameter at the transmitter, by default None
        self.cap_theta = None                               # float,     adimensional beam parameter at the receiver, by default None
        self.cap_theta_bar = None                           # float,     adimensional beam parameter at the receiver, by default None
        self.turbulence_strength = None                     # string,    turbulence strength, by default None
        self.rytov_variance = None                          # float,     Rytov variance, by default None
        self.r_0 = None                                     # float,     Fried parameter (in m), by default None

    def compute_zenith_angle_rad(self):
        self.zenith_angle = 90. - self.elevation
        self.zenith_angle_rad = np.deg2rad(self.zenith_angle)
        return self.zenith_angle_rad

    def compute_propagation_distance(self):
        return (self.SLT_altitude - self.ALT_altitude)/np.cos(self.zenith_angle_rad)

    def compute_r_0(self):
        """" Computes Fried parameter for vertical link
        """" Andrews p.522 eq. 89"""
        rho_0 = (1.45*float(mp.quad(lambda h: self.c2n(h), [self.ALT_altitude,self.SLT_altitude]))*self.k**2)**(-3/5)*np.cos(self.zenith_angle_rad)**(3/5)
        self.r_0 = 2.1*rho_0


    def compute_gaussian_parameters(self):
        """"Computes all necessary parameters for gaussian treatment"""
        """" Andrews p.488 eq. 8"""
        self.cap_lambda_0 = float((2 * self.R) / (self.k * self.W_0 ** 2))

        """" Andrews p.488 eq. 8"""
        self.cap_theta_0 = 1 - self.R/self.F_0
        self.cap_theta_0_bar = 1 - self.cap_theta_0

        """" Andrews p.489 eq. 9"""
        self.cap_lambda = float(self.cap_lambda_0 / (self.cap_theta_0**2 + self.cap_lambda_0**2))

        """" Andrews p.489 eq. 9"""
        self.cap_theta = float(self.cap_theta_0/(self.cap_theta_0**2 + self.cap_lambda_0**2))
        self.cap_theta_bar = 1-self.cap_theta

        """" Andrews p.489 eq. 9"""
        self.W = float(self.W_0*mp.sqrt(self.cap_theta_0**2 + self.cap_lambda_0**2))


    def compute_UL_rytov(self):
        """Computes the rytov variance in the case of an uplink geometry"""
        """" Andrews p.490 eq.14"""
        xi = lambda h: 1 - (h - self.ALT_altitude) / (self.SLT_altitude - self.ALT_altitude)

        """" Andrews p.503 eq. 55 and p.504 eq.58"""
        mu_3u = float(np.real(mp.quad(lambda h:
                              self.refractiveIndexObject.c2n(h)*(xi(h)**(5./6.)*(self.cap_lambda*xi(h) + 1j*(1-self.cap_theta_bar*xi(h)))**(5./6.)-self.cap_lambda**(5./6.)*xi(h)**(5./3.)),
                              [self.ALT_altitude, self.SLT_altitude])))

        rytov_UL = 8.70*mu_3u*self.k**(7 / 6.) * (self.SLT_altitude-self.ALT_altitude)**(5 / 6.)*mp.sec(self.zenith_angle_rad)**(11 / 6.)

        return rytov_UL

    def compute_DL_rytov(self):
        """Computes the rytov variance in the case of a downlink geometry"""
        """" Andrews p.522 eq. 92"""
        mp.mp.dps = 50
        mu_3d = 0.26* float(mp.quad(lambda h:  self.c2n(h) * ((h - self.ALT_altitude)/(self.SLT_altitude - self.ALT_altitude))**(5 / 6.), [self.ALT_altitude, self.SLT_altitude]))
        rytov_DL = 8.70 * mu_3d* self.k**(7 / 6.) * (self.SLT_altitude - self.ALT_altitude)**(5/6)* mp.sec(self.zenith_angle_rad)**(11 / 6.)
        return rytov_DL

    def compute_turbulence_strength(self):
        """Initializes the parameters, computes the rytov variance and determines if the turbulence is weak or strong"""
        self.zenith_angle_rad = self.compute_zenith_angle_rad()
        self.R = self.compute_propagation_distance()
        self.compute_gaussian_parameters()
        self.compute_r_0()
        if self.geometry == 'uplink':
            self.rytov_variance= float(self.compute_UL_rytov())
            if self.rytov_variance < 1 and self.rytov_variance*self.cap_lambda**(5/6) <1:
                self.turbulence_strength = 'weak'
            else:
                self.turbulence_strength = 'strong'

        elif self.geometry == 'downlink':
            self.rytov_variance = float(self.compute_DL_rytov())
            if self.rytov_variance < 1:
                self.turbulence_strength = 'weak'
            else:
                self.turbulence_strength = 'strong'
        else:
            raise Exception("Please enter: 'uplink' or 'downlink' as geometry parameter")



