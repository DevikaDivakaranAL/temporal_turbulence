import mpmath as mp
import numpy as np
from LinkedAttribute import *

class ComputeScintillationIndex(object):
    """Definition of the attributes that need to be taken from another instance so that they are automatically updated if the instance has them changed.
    For more info read comments in LinkedAttribute class"""
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


    def __init__(self,  refractiveIndexObject,turbulenceStrengthObject,beamEffectObject, C_r = 1, D = 0):
        # ACCESSED PARAMETERS
        self.refractiveIndexObject = refractiveIndexObject                                       # ComputeRefractiveIndexStructure instance
        self.c2n = refractiveIndexObject.c2n                                                     # function,  function of c2n (in m**(-2/3)) vs altitude, taken from ComputeRefractiveIndexStructure instance

        self.turbulenceStrengthObject = turbulenceStrengthObject                                 # ComputeTurbulenceStrength instance
        self.turbulence_strength = turbulenceStrengthObject.turbulence_strength                 # string,    turbulence strength, taken from ComputeTurbulenceStrength instance
        self.rytov_variance = turbulenceStrengthObject.rytov_variance                           # float,     Rytov variance, taken from ComputeTurbulenceStrength instance
        self.wavelength = turbulenceStrengthObject.wavelength                                    # float,     laser wavelength (in m), taken from ComputeTurbulenceStrength instance
        self.k = turbulenceStrengthObject.k                                                     # float,     laser wavenumber (in m**-1), taken from ComputeTurbulenceStrength instance
        self.zenith_angle_rad = turbulenceStrengthObject.zenith_angle_rad                        # float,      zenith angle (in rad), taken from ComputeTurbulenceStrength instance
        self.geometry = turbulenceStrengthObject.geometry                                        # string,    link geometry, taken from ComputeTurbulenceStrength instance
        self.ALT_altitude = turbulenceStrengthObject.ALT_altitude                                # float,     ALT altitude (in m), taken from ComputeTurbulenceStrength instance
        self.SLT_altitude = turbulenceStrengthObject.SLT_altitude                                # float,     SLT altitude (in m), taken from ComputeTurbulenceStrength instance
        self.R = turbulenceStrengthObject.R                                                      #  float,     distance of propagation ALT -SLT (in m), taken from ComputeTurbulenceStrength instance
        self.W_0 = turbulenceStrengthObject.W_0                                                 # float,     beam waist at transmitter (in m), taken from ComputeTurbulenceStrength instance
        self.W = turbulenceStrengthObject.W                                                     # float,     diffractive beam radius at the receiver (in m), taken from ComputeTurbulenceStrength instance
        self.cap_theta = turbulenceStrengthObject.cap_theta                                     # float,     adimensional beam parameter at the receiver, , taken from ComputeTurbulenceStrength instance
        self.r_0 = turbulenceStrengthObject.r_0                                                 # float,     Fried parameter (in m), taken from ComputeTurbulenceStrength instance

        self.beamEffectObject = beamEffectObject                                                 # ComputeBeamEffects instance
        self.W_eff = beamEffectObject.W_eff                                                      # float,     long term beam radius at the receiver (in m) including optical turbulence effects, taken from ComputeBeamEffects instance
        self.r2_c = beamEffectObject.r2_c                                                        # float,     beam wander variance at receiver (in m), taken from ComputeBeamEffects instance

        # INPUT PARAMETERS
        self.C_r = C_r                                                                           # Scaling constant, typically in the range of 1 (worst case) to 2pi (best case), by default 1
        self.D = D                                                                               # receiver aperture diameter in m, by default 0 (point receiver)

        # COMPUTED PARAMETERS
        self.sigma2_pe = None                                                                    # float,     beam wander induced pointing error (variance), by default None
        self.alpha_pe = None                                                                     # float,     beam wander induced pointing error normalized by the distance, by default None
        self.scintillation_index = None                                                          # float,     scintillation index (untracked in uplink case), by default None
        self.scintillation_index_averaged = None                                                 # float,     scintillation index averaged by the aperture, by default None
        self.scintillation_index_tracked = None                                                  # float,     scintillation index without beam wander (computed in uplink case), by default None

    def compute_sigma2_pe(self):
        """Computes the beam wander induced pointing variance and its value normalized by the distance"""
        """" Andrews P.503 eq. 53"""
        num = self.C_r**2 *self.W_0**2/self.r_0**2
        den = 1 + num
        self.sigma2_pe = float( 7.25*(self.SLT_altitude - self.ALT_altitude)**2*mp.sec(self.zenith_angle_rad)**3*self.W_0**(-1/3)*float(mp.quad(lambda h:
                                                      self.c2n(h) * (1 - (h - self.ALT_altitude)/(self.SLT_altitude - self.ALT_altitude))**2,
                                                      [self.ALT_altitude, self.SLT_altitude]))*(1 - (num / den)**(1./6.)))

        self.alpha_pe = float(mp.sqrt(self.sigma2_pe)/self.R)

    def compute_avg_scint_index(self):
        """In case of downlink, it computes the scintillation index normalized by the receiver aperture"""
        mp.mp.dps = 50
        # """" Andrews, p.496 eq. 39 """
        scintillation_index_averaged = (8.70 * self.k**(7 / 6.)*(self.SLT_altitude - self.ALT_altitude)**(5/6) * mp.sec(self.zenith_angle_rad)**(11 / 6.)*float(mp.re(mp.quad(lambda h:
                              self.c2n(h)*(((self.k *self.D**2)/(16*self.R) + 1j*(h - self.ALT_altitude)/(self.SLT_altitude -self.ALT_altitude))**(5./6.) - ((self.k *self.D**2)/(16*self.R))**(5./6.)),
                              [self.ALT_altitude, self.SLT_altitude]))))

        return float(scintillation_index_averaged)

    def compute_scintillation_index(self):
        """calls the different functions to compute the scintillation index depending on the case"""
        if self.geometry == 'uplink':
            if self.turbulence_strength == 'weak':
                # TRACKED
                """Andrews p.524 eq.99"""
                self.scintillation_index_tracked =  float(mp.expm1((0.49 * self.rytov_variance) / (
                            1 + (1 + self.cap_theta) * 0.56 * self.rytov_variance ** (12. / 10.)) ** (7. / 6.) + (
                                   0.51 * self.rytov_variance) / (1 + 0.69 * self.rytov_variance ** (12. / 10.)) ** (
                                   5. / 6.)))

                # UNTRACKED + beam wander
                self.compute_sigma2_pe()

                """" Andrews p. 503 eq. 54 and following"""
                self.scintillation_index = float(5.95*(self.SLT_altitude-self.ALT_altitude)**2*mp.sec(self.zenith_angle_rad)**2*((2*self.W_0)/self.r_0)**(5./3.)*(self.alpha_pe/self.W)**2 + self.scintillation_index_tracked)

            else:
                """" Andrews p. 506 eq. 60 and eq. 61"""
                # TRACKED
                self.scintillation_index_tracked =  float(mp.expm1((0.49*self.rytov_variance)/(1+(1 + self.cap_theta)*0.56*self.rytov_variance**(12./10.))**(7./6.) + (0.51*self.rytov_variance)/(1+0.69*self.rytov_variance**(12./10.))**(5./6.)))

                # UNTRACKED + beam wander
                self.compute_sigma2_pe()
                term1 = 5.95 * (self.SLT_altitude - self.ALT_altitude)**2 * mp.sec(self.zenith_angle_rad) ** 2 * (2 * self.W_0 / self.r_0)**(5. / 3.) * (self.alpha_pe / self.W) ** 2

                self.scintillation_index = float(term1+ self.scintillation_index_tracked)

        elif self.geometry == 'downlink':
            if self.turbulence_strength == 'weak':
                """" Andrews p. 495 eq. 38"""
                self.scintillation_index = float(self.rytov_variance)
                """" In the weak case we can also compute the scintillation index averaged over the aperture"""
                self.scintillation_index_averaged = float(self.compute_avg_scint_index())
            else:
                """" Andrews p. 497 eq. 40"""
                self.scintillation_index =  float(mp.expm1((0.49*self.rytov_variance)/(1+1.11*mp.sqrt(self.rytov_variance)**(12./5.))**(7./6.) + (0.51*self.rytov_variance)/(1+0.69*mp.sqrt(self.rytov_variance)**(12./5.))**(5./6.)))