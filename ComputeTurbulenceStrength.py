import numpy as np
import mpmath as mp
from LinkedAttribute import *

class ComputeTurbulenceStrength(object):
    """
    This class calculates the turbulence strength for optical links between ground stations and satellites.
    It utilizes the refractive index structure parameter to assess the impact of atmospheric turbulence on the link.
    """

    # Linking the refractive index attribute from another instance
    c2n = LinkedAttribute('refractiveIndexObject')

    def __init__(self, refractiveIndexObject, wavelength=1550e-9, elevation=90., geometry='uplink', ALT_altitude=4e3, SLT_altitude=600e3, W_0=None, F_0=None):
        # Initializing input and accessed parameters
        self.refractiveIndexObject = refractiveIndexObject  # Instance of ComputeRefractiveIndexStructure
        self.c2n = refractiveIndexObject.c2n                # Function of c2n vs altitude, accessed from ComputeRefractiveIndexStructure

        # Setting default and input parameters
        self.wavelength = wavelength                        # Beam wavelength in meters
        self.elevation = elevation                          # Elevation angle in degrees
        self.geometry = geometry                            # Link geometry ('uplink' or 'downlink')
        self.ALT_altitude = ALT_altitude                    # ALT altitude in meters
        self.SLT_altitude = SLT_altitude                    # SLT altitude in meters
        self.W_0 = W_0                                      # Beam waist at transmitter in meters
        self.F_0 = F_0                                      # Phase front radius of curvature in meters

        # Initializing computed parameters
        self.zenith_angle = 90 - elevation                  # Zenith angle in degrees
        self.zenith_angle_rad = self.compute_zenith_angle_rad() # Zenith angle in radians
        self.k = (2 * np.pi) / self.wavelength              # Beam wavenumber
        # The following parameters are to be computed
        self.R = None                                       # Propagation distance ALT - SLT
        self.W = None                                       # Diffractive beam radius at receiver
        # Adimensional beam parameters at the transmitter and receiver
        self.cap_lambda_0 = None
        self.cap_lambda = None
        self.cap_theta_0 = None
        self.cap_theta_0_bar = None
        self.cap_theta = None
        self.cap_theta_bar = None
        # Turbulence parameters
        self.turbulence_strength = None
        self.rytov_variance = None
        self.r_0 = None

    def compute_zenith_angle_rad(self):
        """ Computes the zenith angle in radians. """
        self.zenith_angle = 90. - self.elevation
        self.zenith_angle_rad = np.deg2rad(self.zenith_angle)
        return self.zenith_angle_rad

    def compute_propagation_distance(self):
        """ Computes the distance of propagation between ALT and SLT. """
        return (self.SLT_altitude - self.ALT_altitude) / np.cos(self.zenith_angle_rad)

    def compute_r_0(self):
        """
        Computes the Fried parameter for a vertical link.
        This parameter characterizes the atmospheric coherence length.
        Reference: Andrews, Laser Beam Propagation through Random Media (Eq. 89 on p.522)
        """
        rho_0 = (1.45 * float(mp.quad(lambda h: self.c2n(h), [self.ALT_altitude, self.SLT_altitude])) * self.k**2)**(-3/5) * np.cos(self.zenith_angle_rad)**(3/5)
        self.r_0 = 2.1 * rho_0

    def compute_gaussian_parameters(self):
        """
        Computes all necessary Gaussian beam parameters.
        These include adimensional parameters at the transmitter and receiver, and the beam radius at the receiver.
        References to equations from Andrews are included in each calculation step.
        """
        # Computing adimensional beam parameter at the transmitter
        self.cap_lambda_0 = float((2 * self.R) / (self.k * self.W_0 ** 2)) # Eq. 8 on p.488

        # Computing adimensional parameters related to phase front radius
        self.cap_theta_0 = 1 - self.R / self.F_0
        self.cap_theta_0_bar = 1 - self.cap_theta_0

        # Computing adimensional parameters at the receiver
        self.cap_lambda = float(self.cap_lambda_0 / (self.cap_theta_0**2 + self.cap_lambda_0**2)) # Eq. 9 on p.489
        self.cap_theta = float(self.cap_theta_0 / (self.cap_theta_0**2 + self.cap_lambda_0**2))
        self.cap_theta_bar = 1 - self.cap_theta

        # Computing the diffractive beam radius at the receiver
        self.W = float(self.W_0 * mp.sqrt(self.cap_theta_0**2 + self.cap_lambda_0**2))

    def compute_UL_rytov(self):
        """
        Computes the Rytov variance for an uplink geometry.
        The Rytov variance is a measure of the strength of turbulence-induced scintillations in the link.
        Reference: Andrews, Laser Beam Propagation through Random Media (Eq. 14 on p.490)
        """
        xi = lambda h: 1 - (h - self.ALT_altitude) / (self.SLT_altitude - self.ALT_altitude)
        mu_3u = float(np.real(mp.quad(lambda h: self.refractiveIndexObject.c2n(h) * (xi(h)**(5./6.) * (self.cap_lambda * xi(h) + 1j * (1 - self.cap_theta_bar * xi(h)))**(5./6.) - self.cap_lambda**(5./6.) * xi(h)**(5./3.)), [self.ALT_altitude, self.SLT_altitude])))

        rytov_UL = 8.70 * mu_3u * self.k**(7 / 6.) * (self.SLT_altitude - self.ALT_altitude)**(5 / 6.) * mp.sec(self.zenith_angle_rad)**(11 / 6.)
        return rytov_UL

    def compute_DL_rytov(self):
        """
        Computes the Rytov variance for a downlink geometry.
        Similar to the uplink, it assesses turbulence-induced scintillations but for a downlink scenario.
        Reference: Andrews, Laser Beam Propagation through Random Media (Eq. 92 on p.522)
        """
        mp.mp.dps = 50  # Setting the decimal precision
        mu_3d = 0.26 * float(mp.quad(lambda h: self.c2n(h) * ((h - self.ALT_altitude) / (self.SLT_altitude - self.ALT_altitude))**(5 / 6.), [self.ALT_altitude, self.SLT_altitude]))
        rytov_DL = 8.70 * mu_3d * self.k**(7 / 6.) * (self.SLT_altitude - self.ALT_altitude)**(5/6) * mp.sec(self.zenith_angle_rad)**(11 / 6.)
        return rytov_DL

    def compute_turbulence_strength(self):
        """
        Initializes parameters, computes the Rytov variance, and determines if the turbulence is weak or strong.
        This classification is crucial for understanding the impact of atmospheric turbulence on the optical link.
        """
        # Initializing and computing necessary parameters
        self.zenith_angle_rad = self.compute_zenith_angle_rad()
        self.R = self.compute_propagation_distance()
        self.compute_gaussian_parameters()
        self.compute_r_0()

        # Computing Rytov variance and assessing turbulence strength
        if self.geometry == 'uplink':
            self.rytov_variance = float(self.compute_UL_rytov())
            if self.rytov_variance < 1 and self.rytov_variance * self.cap_lambda**(5/6) < 1:
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
            raise Exception("Invalid geometry parameter. Please use 'uplink' or 'downlink'.")
