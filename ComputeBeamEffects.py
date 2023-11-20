import mpmath as mp
from LinkedAttribute import *


class ComputeBeamEffects(object):
    """
    Computes the beam effects including effective beam radius and beam wander due to atmospheric turbulence.
    Utilizes linked attributes from refractive index and turbulence strength objects.
    """
    c2n = LinkedAttribute('refractiveIndexObject')
    geometry = LinkedAttribute('turbulenceStrengthObject')
    ALT_altitude = LinkedAttribute('turbulenceStrengthObject')
    SLT_altitude = LinkedAttribute('turbulenceStrengthObject')
    zenith_angle_rad = LinkedAttribute('turbulenceStrengthObject')
    wavelength = LinkedAttribute('turbulenceStrengthObject')
    W_0 = LinkedAttribute('turbulenceStrengthObject')
    W = LinkedAttribute('turbulenceStrengthObject')
    r_0 = LinkedAttribute('turbulenceStrengthObject')
    cap_lambda = LinkedAttribute('turbulenceStrengthObject')
    k = LinkedAttribute('turbulenceStrengthObject')
    turbulence_strength = LinkedAttribute('turbulenceStrengthObject')

    def __init__(self, refractiveIndexObject,turbulenceStrengthObject):
        """
               Initializes the ComputeBeamEffects class.

               Parameters:
               - refractiveIndexObject: Instance of ComputeRefractiveIndexStructure.
               - turbulenceStrengthObject: Instance of ComputeTurbulenceStrength.
               """
        # ACCESSED PARAMETERS
        self.refractiveIndexObject = refractiveIndexObject                  # ComputeRefractiveIndexStructure instance
        self.c2n = refractiveIndexObject.c2n                                # function,  function of c2n (in m**(-2/3)) vs altitude of height, taken from ComputeRefractiveIndexStructure instance

        self.turbulenceStrengthObject = turbulenceStrengthObject           # ComputeTurbulenceStrength instance
        self.geometry = turbulenceStrengthObject.geometry                  # string,    link geometry, taken from ComputeTurbulenceStrength instance
        self.ALT_altitude = turbulenceStrengthObject.ALT_altitude          # float,     ALT altitude (in m), taken from ComputeTurbulenceStrength instance
        self.SLT_altitude = turbulenceStrengthObject.SLT_altitude          # float,     SLT altitude (in m), taken from ComputeTurbulenceStrength instance
        self.zenith_angle_rad = turbulenceStrengthObject.zenith_angle_rad  #float,     zenith angle (in rad), taken from ComputeTurbulenceStrength instance
        self.wavelength = turbulenceStrengthObject.wavelength              # float,     laser wavelength (in m), taken from ComputeTurbulenceStrength instance
        self.W_0 = turbulenceStrengthObject.W_0                            # float,     beam waist at transmitter (in m), taken from ComputeTurbulenceStrength instance
        self.W = turbulenceStrengthObject.W                                # float,     diffractive beam radius at the receiver (in m), taken from ComputeTurbulenceStrength instance
        self.r_0 = turbulenceStrengthObject.r_0                            # float,     Fried parameter (in m), taken from ComputeTurbulenceStrength instance
        self.cap_lambda = turbulenceStrengthObject.cap_lambda              # float,     adimensional beam parameter at the receiver, taken from ComputeTurbulenceStrength instance
        self.k = turbulenceStrengthObject.k                                # float,     laser wavenumber (in m**-1), taken from ComputeTurbulenceStrength instance
        self.turbulence_strength = turbulenceStrengthObject.turbulence_strength  # string,    turbulence strength, taken from ComputeTurbulenceStrength instance

        # COMPUTED PARAMETERS
        self.W_eff = None                                                   # float,    long term beam radius at the receiver (in m) including optical turbulence effects, by default None
        self.r2_c = None                                                    # float,    beam wander variance at receiver (in m), by default None

    def compute_r2_c(self):
        """Computes the beam wander displacement variance """

        """" Andrews et al.: Strehl ratio and scintillation theory, eq. 20 """
        self.r2_c = float(0.54 * (self.SLT_altitude - self.ALT_altitude)**2 * mp.sec(self.zenith_angle_rad)**2 * (
                    self.wavelength / (2 * self.W_0))**2 * (2 * self.W_0 / self.r_0)**(5./3.))


    def compute_W_eff(self):
        """Computes the long-term spot size for the uplink channel """

        # Strong and weak case
        """" Andrews p.500 eq. 48 """
        D_0 = mp.sqrt(8) * self.W_0
        self.W_eff = float(self.W*(1+ (D_0/self.r_0)**(5/3))**(3/5))

    def compute_beam_effects(self):
        """
        Determines the beam effects based on the communication geometry and turbulence strength.
        Assigns computed values to the appropriate attributes.
        """
        if self.geometry == 'uplink':
            self.compute_W_eff()
            if self.turbulence_strength == 'weak':
                self.compute_r2_c()
        else:
            self.W_eff = self.W # For downlink, the effective beam radius is the same as the diffractive beam radius
