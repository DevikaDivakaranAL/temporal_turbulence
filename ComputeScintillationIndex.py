import mpmath as mp
import numpy as np
from LinkedAttribute import *

class ComputeScintillationIndex(object):
    """
    This class computes the scintillation index, which quantifies the intensity fluctuations of a laser beam
    due to atmospheric turbulence. It utilizes linked attributes from instances of refractive index,
    turbulence strength, and beam effect objects.
    """

    # Linking attributes from other classes for direct access
    c2n = LinkedAttribute('refractiveIndexObject')
    turbulence_strength = LinkedAttribute('turbulenceStrengthObject')
    rytov_variance = LinkedAttribute('turbulenceStrengthObject')
    wavelength = LinkedAttribute('turbulenceStrengthObject')
    k = LinkedAttribute('turbulenceStrengthObject')
    zenith_angle_rad = LinkedAttribute('turbulenceStrengthObject')
    geometry = LinkedAttribute('turbulenceStrengthObject')
    ALT_altitude = LinkedAttribute('turbulenceStrengthObject')
    SLT_altitude = LinkedAttribute('turbulenceStrengthObject')
    R = LinkedAttribute('turbulenceStrengthObject')
    W_0 = LinkedAttribute('turbulenceStrengthObject')
    W = LinkedAttribute('turbulenceStrengthObject')
    cap_theta = LinkedAttribute('turbulenceStrengthObject')
    r_0 = LinkedAttribute('turbulenceStrengthObject')
    W_eff = LinkedAttribute('beamEffectObject')
    r2_c = LinkedAttribute('beamEffectObject')

    def __init__(self, refractiveIndexObject, turbulenceStrengthObject, beamEffectObject, C_r, D):
        """
        Initializes the ComputeScintillationIndex class.

        Parameters:
        - refractiveIndexObject: Instance of ComputeRefractiveIndexStructure.
        - turbulenceStrengthObject: Instance of ComputeTurbulenceStrength.
        - beamEffectObject: Instance of ComputeBeamEffects.
        - C_r: Scaling constant, typically ranges from 1 (worst case) to 2π (best case).
        - D: Receiver aperture diameter in meters.
        """
        # Accessing parameters from other instances
        self.refractiveIndexObject = refractiveIndexObject  # Instance for refractive index calculations
        self.c2n = refractiveIndexObject.c2n                # Function of c2n vs altitude

        # Turbulence strength parameters
        self.turbulenceStrengthObject = turbulenceStrengthObject
        self.turbulence_strength = turbulenceStrengthObject.turbulence_strength
        self.rytov_variance = turbulenceStrengthObject.rytov_variance
        self.wavelength = turbulenceStrengthObject.wavelength
        self.k = turbulenceStrengthObject.k
        self.zenith_angle_rad = turbulenceStrengthObject.zenith_angle_rad
        self.geometry = turbulenceStrengthObject.geometry
        self.ALT_altitude = turbulenceStrengthObject.ALT_altitude
        self.SLT_altitude = turbulenceStrengthObject.SLT_altitude
        self.R = turbulenceStrengthObject.R
        self.W_0 = turbulenceStrengthObject.W_0
        self.W = turbulenceStrengthObject.W
        self.cap_theta = turbulenceStrengthObject.cap_theta
        self.r_0 = turbulenceStrengthObject.r_0

        # Beam effect parameters
        self.beamEffectObject = beamEffectObject
        self.W_eff = beamEffectObject.W_eff
        self.r2_c = beamEffectObject.r2_c

        # Input parameters
        self.C_r = C_r                                      # Scaling constant
        self.D = D                                          # Receiver aperture diameter

        # Parameters to be computed
        self.sigma2_pe = None                               # Beam wander induced pointing error variance
        self.alpha_pe = None                                # Normalized pointing error
        self.scintillation_index = None                     # Scintillation index (untracked in uplink case)
        self.scintillation_index_averaged = None            # Averaged scintillation index over the aperture
        self.scintillation_index_tracked = None             # Scintillation index without beam wander (computed in uplink case)

    def compute_sigma2_pe(self):
        """
        Computes the beam wander induced pointing error variance and its normalized value by the distance.
        This calculation is important for understanding the impact of beam wander on pointing accuracy.
        Reference: Andrews, Laser Beam Propagation through Random Media, P.503 (Eq. 53)
        """
        num = self.C_r**2 * self.W_0**2 / self.r_0**2
        den = 1 + num
        self.sigma2_pe = float(7.25 * (self.SLT_altitude - self.ALT_altitude)**2 * mp.sec(self.zenith_angle_rad)**3 * self.W_0**(-1/3) * float(mp.quad(lambda h: self.c2n(h) * (1 - (h - self.ALT_altitude) / (self.SLT_altitude - self.ALT_altitude))**2, [self.ALT_altitude, self.SLT_altitude])) * (1 - (num / den)**(1./6.)))
        self.alpha_pe = float(mp.sqrt(self.sigma2_pe) / self.R)

    def compute_avg_scint_index(self):
        """
        Computes the scintillation index normalized by the receiver aperture for downlink scenarios.
        This calculation helps assess the impact of receiver size on the intensity fluctuations of the received signal.
        Reference: Andrews, Laser Beam Propagation through Random Media, p.496 (Eq. 39)
        """
        scintillation_index_averaged = (8.70 * self.k**(7 / 6.) * (self.SLT_altitude - self.ALT_altitude)**(5/6) * mp.sec(self.zenith_angle_rad)**(11 / 6.) * float(mp.re(mp.quad(lambda h: self.c2n(h) * (((self.k * self.D**2) / (16 * self.R) + 1j * (h - self.ALT_altitude) / (self.SLT_altitude - self.ALT_altitude))**(5./6.) - ((self.k * self.D**2) / (16 * self.R))**(5./6.)), [self.ALT_altitude, self.SLT_altitude]))))
        return float(scintillation_index_averaged)

    def compute_scintillation_index(self):
        """
        Computes the scintillation index depending on the link geometry and turbulence strength.
        This index quantifies the intensity fluctuations of the beam due to atmospheric turbulence.
        Different calculations are applied for uplink and downlink scenarios, and for weak and strong turbulence regimes.
        """
        if self.geometry == 'uplink':
            if self.turbulence_strength == 'weak':
                # TRACKED
                """Andrews, Laser Beam Propagation through Random Media p.524 eq.99"""
                self.scintillation_index_tracked =  float(mp.expm1((0.49 * self.rytov_variance) / (
                            1 + (1 + self.cap_theta) * 0.56 * self.rytov_variance ** (12. / 10.)) ** (7. / 6.) + (
                                   0.51 * self.rytov_variance) / (1 + 0.69 * self.rytov_variance ** (12. / 10.)) ** (
                                   5. / 6.)))

                # UNTRACKED + beam wander
                self.compute_sigma2_pe()

                """" Andrews, Laser Beam Propagation through Random Media p. 503 eq. 54 and following"""
                self.scintillation_index = float(5.95*(self.SLT_altitude-self.ALT_altitude)**2*mp.sec(self.zenith_angle_rad)**2*((2*self.W_0)/self.r_0)**(5./3.)*(self.alpha_pe/self.W)**2 + self.scintillation_index_tracked)

            else:
                """" Andrews, Laser Beam Propagation through Random Media p. 506 eq. 60 and eq. 61"""
                # TRACKED
                self.scintillation_index_tracked =  float(mp.expm1((0.49*self.rytov_variance)/(1+(1 + self.cap_theta)*0.56*self.rytov_variance**(12./10.))**(7./6.) + (0.51*self.rytov_variance)/(1+0.69*self.rytov_variance**(12./10.))**(5./6.)))

                # UNTRACKED + beam wander
                self.compute_sigma2_pe()
                term1 = 5.95 * (self.SLT_altitude - self.ALT_altitude)**2 * mp.sec(self.zenith_angle_rad) ** 2 * (2 * self.W_0 / self.r_0)**(5. / 3.) * (self.alpha_pe / self.W) ** 2

                self.scintillation_index = float(term1+ self.scintillation_index_tracked)

        elif self.geometry == 'downlink':
            if self.turbulence_strength == 'weak':
                """" Andrews, Laser Beam Propagation through Random Media p. 495 eq. 38"""
                self.scintillation_index = float(self.rytov_variance)
                """" In the weak case we can also compute the scintillation index averaged over the aperture"""
                self.scintillation_index_averaged = float(self.compute_avg_scint_index())
            else:
                """" Andrews, Laser Beam Propagation through Random Media p. 497 eq. 40"""
                self.scintillation_index =  float(mp.expm1((0.49*self.rytov_variance)/(1+1.11*mp.sqrt(self.rytov_variance)**(12./5.))**(7./6.) + (0.51*self.rytov_variance)/(1+0.69*mp.sqrt(self.rytov_variance)**(12./5.))**(5./6.)))