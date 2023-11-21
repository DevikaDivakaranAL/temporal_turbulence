import mpmath as mp
from LinkedAttribute import *

class ComputeBeamEffects(object):
    """
    This class computes the beam effects, including effective beam radius and beam wander, due to atmospheric turbulence.
    It utilizes linked attributes from instances of refractive index and turbulence strength objects to calculate these effects.
    """

    # Linking attributes from other classes for direct access
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

    def __init__(self, refractiveIndexObject, turbulenceStrengthObject):
        """
        Initializes the ComputeBeamEffects class.

        Parameters:
        - refractiveIndexObject: Instance of ComputeRefractiveIndexStructure.
        - turbulenceStrengthObject: Instance of ComputeTurbulenceStrength.
        """
        # Accessing parameters from other instances
        self.refractiveIndexObject = refractiveIndexObject                  # Instance for refractive index calculations
        self.c2n = refractiveIndexObject.c2n                                # Function of c2n vs altitude

        self.turbulenceStrengthObject = turbulenceStrengthObject            # Instance for turbulence strength calculations
        self.geometry = turbulenceStrengthObject.geometry                   # Link geometry (uplink or downlink)
        self.ALT_altitude = turbulenceStrengthObject.ALT_altitude           # ALT altitude in meters
        self.SLT_altitude = turbulenceStrengthObject.SLT_altitude           # SLT altitude in meters
        self.zenith_angle_rad = turbulenceStrengthObject.zenith_angle_rad   # Zenith angle in radians
        self.wavelength = turbulenceStrengthObject.wavelength               # Laser wavelength in meters
        self.W_0 = turbulenceStrengthObject.W_0                             # Beam waist at transmitter in meters
        self.W = turbulenceStrengthObject.W                                 # Diffractive beam radius at receiver in meters
        self.r_0 = turbulenceStrengthObject.r_0                             # Fried parameter in meters
        self.cap_lambda = turbulenceStrengthObject.cap_lambda               # Adimensional beam parameter at receiver
        self.k = turbulenceStrengthObject.k                                 # Laser wavenumber
        self.turbulence_strength = turbulenceStrengthObject.turbulence_strength  # Turbulence strength category

        # Parameters to be computed
        self.W_eff = None                                                   # Long-term beam radius at receiver, including optical turbulence effects
        self.r2_c = None                                                    # Beam wander variance at receiver

    def compute_r2_c(self):
        """
        Computes the beam wander displacement variance.
        This variance characterizes the lateral displacement of the beam centroid due to atmospheric turbulence.
        Reference: Andrews et al., Strehl ratio and scintillation theory (Eq. 20)
        """
        self.r2_c = float(0.54 * (self.SLT_altitude - self.ALT_altitude)**2 * mp.sec(self.zenith_angle_rad)**2 *
                          (self.wavelength / (2 * self.W_0))**2 * (2 * self.W_0 / self.r_0)**(5./3.))

    def compute_W_eff(self):
        """
        Computes the long-term spot size for the uplink channel.
        The effective spot size takes into account the spreading of the beam due to atmospheric turbulence.
        Reference: Andrews, Laser Beam Propagation through Random Media (Eq. 48 on p.500)
        """
        D_0 = mp.sqrt(8) * self.W_0
        self.W_eff = float(self.W * (1 + (D_0 / self.r_0)**(5/3))**(3/5))

    def compute_beam_effects(self):
        """
        Determines the beam effects based on the communication geometry and turbulence strength.
        This includes calculating the effective beam radius and beam wander variance for uplink scenarios.
        For downlink, the effective beam radius is assumed to be the same as the diffractive beam radius.
        """
        if self.geometry == 'uplink':
            self.compute_W_eff()
            if self.turbulence_strength == 'weak':
                self.compute_r2_c()
        else:
            # In downlink geometry, the effective beam radius is assumed to be the same as the diffractive beam radius
            self.W_eff = self.W
