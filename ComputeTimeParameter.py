import numpy as np
import mpmath as mp
from LinkedAttribute import *


class ComputeTimeParameters(object):
    """Definition of the attributes that need to be taken from another instance so that they are automatically updated if the instance has them changed.
    For more info read comments in LinkedAttribute class"""
    c2n=LinkedAttribute('refractiveIndexObject')
    rms_wind_speed = LinkedAttribute('refractiveIndexObject')
    turbulence_strength = LinkedAttribute('turbulenceStrengthObject')
    k=LinkedAttribute('turbulenceStrengthObject')
    zenith_angle_rad=LinkedAttribute('turbulenceStrengthObject')
    geometry=LinkedAttribute('turbulenceStrengthObject')
    ALT_altitude=LinkedAttribute('turbulenceStrengthObject')
    SLT_altitude=LinkedAttribute('turbulenceStrengthObject')
    R=LinkedAttribute('turbulenceStrengthObject')
    cap_lambda=LinkedAttribute('turbulenceStrengthObject')
    cap_theta_bar = LinkedAttribute('turbulenceStrengthObject')
    scintillation_index=LinkedAttribute('scintillationIndexObject')
    scintillation_index_averaged=LinkedAttribute('scintillationIndexObject')
    scintillation_index_tracked=LinkedAttribute('scintillationIndexObject')
    alpha = LinkedAttribute('pdfObject')
    beta = LinkedAttribute('pdfObject')
    alpha_avg = LinkedAttribute('pdfObject')
    beta_avg = LinkedAttribute('pdfObject')
    alpha_tracked = LinkedAttribute('pdfObject')
    beta_tracked = LinkedAttribute('pdfObject')
    mean_intensity = LinkedAttribute('pdfObject')
    F_T = LinkedAttribute('intensityParametersObject')
    S_T = LinkedAttribute('intensityParametersObject')
    fade_probability = LinkedAttribute('intensityParametersObject')
    fade_probability_averaged = LinkedAttribute('intensityParametersObject')
    fade_probability_tracked = LinkedAttribute('intensityParametersObject')
    surge_probability = LinkedAttribute('intensityParametersObject')
    surge_probability_averaged = LinkedAttribute('intensityParametersObject')
    surge_probability_tracked = LinkedAttribute('intensityParametersObject')
    compute_only_fades = LinkedAttribute('intensityParametersObject')


    def __init__(self,  refractiveIndexObject,turbulenceStrengthObject, scintillationIndexObject, pdfObject, intensityParametersObject ):
        # ACCESSED PARAMETERS
        self.refractiveIndexObject =refractiveIndexObject                            # ComputeRefractiveIndexStructure instance
        self.c2n = refractiveIndexObject.c2n                                         # function,  function of c2n (in m**(-2/3)) vs altitude, taken from ComputeRefractiveIndexStructure instance
        self.rms_wind_speed = refractiveIndexObject.wind_speed                       #float,    rms of wind data (in m/s), taken from ComputeRefractiveIndexStructure instance

        self.turbulenceStrengthObject =turbulenceStrengthObject                      # ComputeTurbulenceStrength instance
        self.turbulence_strength = turbulenceStrengthObject.turbulence_strength      # string,    turbulence strength, taken from ComputeTurbulenceStrength instance
        self.k = turbulenceStrengthObject.k                                          # float,     laser wavenumber (in m**-1), taken from ComputeTurbulenceStrength instance
        self.zenith_angle_rad = turbulenceStrengthObject.zenith_angle_rad            # float,      zenith angle (in rad), taken from ComputeTurbulenceStrength instance
        self.geometry = turbulenceStrengthObject.geometry                            # string,    link geometry, taken from ComputeTurbulenceStrength instance
        self.ALT_altitude = turbulenceStrengthObject.ALT_altitude                    # float,     ALT altitude (in m), taken from ComputeTurbulenceStrength instance
        self.SLT_altitude = turbulenceStrengthObject.SLT_altitude                    # float,     SLT altitude (in m), taken from ComputeTurbulenceStrength instance
        self.R = turbulenceStrengthObject.R                                          # float,     distance of propagation ALT -SLT (in m), taken from ComputeTurbulenceStrength instance
        self.cap_lambda = turbulenceStrengthObject.cap_lambda                        # float,     adimensional beam parameter at the receiver, taken from ComputeTurbulenceStrength instance
        self.cap_theta_bar = turbulenceStrengthObject.cap_theta_bar                  # float,     adimensional beam parameter at the receiver, taken from ComputeTurbulenceStrength instance

        self.scintillationIndexObject = scintillationIndexObject                     # ComputeScintillationIndex instance
        self.scintillation_index = scintillationIndexObject.scintillation_index      # float,     scintillation index (untracked in uplink case), taken from ComputeScintillationIndex instance
        self.scintillation_index_averaged = scintillationIndexObject.scintillation_index_averaged  # float,     scintillation index averaged by the aperture, taken from ComputeScintillationIndex instance
        self.scintillation_index_tracked = scintillationIndexObject.scintillation_index_tracked  # float,     scintillation index without beam wander (computed in uplink case), taken from ComputeScintillationIndex instance

        self.pdfObject = pdfObject                                              # ComputeProbabilityDensityFunction instance
        self.alpha = pdfObject.alpha                                            # float,     parameter of the gamma-gamma distribution, taken from ComputeProbabilityDensityFunction instance
        self.beta = pdfObject.beta                                              # float,     parameter of the gamma-gamma distribution, taken from ComputeProbabilityDensityFunction instance
        self.alpha_avg = pdfObject.alpha_avg                                    # float,     parameter of the gamma-gamma distribution computed using averaged scintillation index,taken from ComputeProbabilityDensityFunction instance
        self.beta_avg = pdfObject.beta_avg                                      # float,     parameter of the gamma-gamma distribution computed using averaged scintillation index,taken from ComputeProbabilityDensityFunction instance
        self.alpha_tracked = pdfObject.alpha_tracked                            # float,     parameter of the gamma-gamma distribution computed using tracked scintillation index,taken from ComputeProbabilityDensityFunction instance
        self.beta_tracked = pdfObject.beta_tracked                              # float,     parameter of the gamma-gamma distribution computed using tracked scintillation index,taken from ComputeProbabilityDensityFunction instance
        self.mean_intensity = pdfObject.mean_intensity                          # mean intensity at the center of the receiver aperture when turbulence is considered,taken from ComputeProbabilityDensityFunction instance

        self.intensityParametersObject = intensityParametersObject              # ComputeIntensityParameters instance
        self.F_T = intensityParametersObject.F_T                                # float, fade threshold in dB, taken from ComputeIntensityParameters instance
        self.S_T = intensityParametersObject.S_T                                # float, surge threshold in dB, taken from ComputeIntensityParameters instance
        self.fade_probability = intensityParametersObject.fade_probability      # computed probability associated to F_T, taken from ComputeIntensityParameters instance.
        self.fade_probability_averaged = intensityParametersObject.fade_probability_averaged    # computed probability associated to F_T_averaged, taken from ComputeIntensityParameters instance .
        self.fade_probability_tracked = intensityParametersObject.fade_probability_tracked      # computed probability associated to F_T_tracked,  taken from ComputeIntensityParameters instance.
        self.fade_probability_multiple = intensityParametersObject.fade_probability_multiple    # computed probability associated to F_T_multiple,  taken from ComputeIntensityParameters instance.
        self.surge_probability = intensityParametersObject.surge_probability                    # computed probability associated to S_T, taken from ComputeIntensityParameters instance.
        self.surge_probability_averaged = intensityParametersObject.surge_probability_averaged  # computed probability associated to S_T_averaged, taken from ComputeIntensityParameters instance .
        self.surge_probability_tracked = intensityParametersObject.surge_probability_tracked    # computed probability associated to S_T_tracked,  taken from ComputeIntensityParameters instance.
        self.compute_only_fades = intensityParametersObject.compute_only_fades                  # parameter to determine if we want to compute only fade statistics or surges too, taken from ComputeIntensityParameters instance.


        #COMPUTED PARAMETERS
        self.transverse_wind = self.rms_wind_speed*mp.sin(np.pi/2 - self.zenith_angle_rad) # float, average transverse wind speed in m/s
        self.nu_0 = None                                                        # float, quasi-frequency
        self.avg_number_fades = None                                            # float, average number of fades per second
        self.avg_number_fades_aperture = None                                   # float, average number of fades per second for aperture averaged scintillation index
        self.avg_number_fades_tracked = None                                    # float, average number of fades per second for tracked scintillation index

        self.avg_fade_duration = None                                           # float, average duration of a fade in seconds
        self.avg_fade_duration_aperture = None                                  # float, average duration of a fade in seconds for aperture averaged scintillation index
        self.avg_fade_duration_tracked =None                                    # float, average duration of a fade in seconds for tracked scintillation index

        self.avg_number_surges = None                                           # float, average number of surges per second
        self.avg_number_surges_aperture = None                                  # float, average number of surges per second for aperture averaged scintillation index
        self.avg_number_surges_tracked = None                                   # float, average number of surges per second for aperture averaged scintillation index

        self.avg_surges_duration = None                                         # float, average duration of a surge in seconds
        self.avg_surges_duration_aperture = None                                # float, average duration of a surge in seconds for aperture averaged scintillation index
        self.avg_surges_duration_tracked = None                                 # float, average duration of a surge in seconds for tracked scintillation index


    def compute_nu_0_DL(self):
        # The quasi frequency is best characterized as the width
        # (standard deviation) of the normalized irradiance power spectrum.
        # "for 'tau = 0', the temporal covariance function reduces to the scintillation index",
        # we are in a situation where tau == 0. Andrews (1995), eq36/37/38
        b_i = self.scintillation_index

        """Andrews, p.514 eq. 76 and 77"""
        mu_5d = float(mp.quad(lambda h: self.c2n(h) * ((self.SLT_altitude - self.ALT_altitude) / (h - self.ALT_altitude)) ** (1 / 6.),
            [self.ALT_altitude, self.SLT_altitude]))
        b_i_deriv = -3.63*self.k**(7/6)*mu_5d*(self.SLT_altitude - self.ALT_altitude)**(5/6)*mp.sec(self.zenith_angle_rad)**(11/6)*(self.k*self.transverse_wind**2/self.R)

        """Andrews, p.514 eq. 73"""
        self.nu_0 =  1 / (2 * np.pi) * (- b_i_deriv/b_i )**(1/2)

    def compute_nu_0_UL(self):
        # The quasi frequency is best characterized as the width
        # (standard deviation) of the normalized irradiance power spectrum.
        # "for 'tau = 0', the temporal covariance function reduces to the scintillation index",
        # we are in a situation where tau == 0. Andrews (1995), eq36/37/38
        b_i = self.scintillation_index        # the temporal covariance for tau = 0 reduces to the scintillation index.

        """Andrews, p.519 eq. 82 and 84"""
        xi = lambda h: 1 - (h - self.ALT_altitude) / (self.SLT_altitude - self.ALT_altitude)
        mu_5u = float(np.real( mp.quad(lambda h: self.c2n(h)*(1/(self.cap_lambda**(1/6)*xi(h)**(1/3))-1/(xi(h)**(1/6)*(self.cap_lambda*xi(h)+1j*(1-self.cap_theta_bar*xi(h)))**(1/6))) ,
            [self.ALT_altitude, self.SLT_altitude])))
        b_i_deriv = -3.63*self.k**(7/6)*mu_5u*(self.SLT_altitude - self.ALT_altitude)**(5/6)*mp.sec(self.zenith_angle_rad)**(11/6)*(self.k*self.transverse_wind**2/self.R)

        """Andrews, p.514 eq. 73"""
        self.nu_0 = 1 / (2 * np.pi) * (-  b_i_deriv/b_i)**(1/2)

    def compute_avg_number_fades(self, alpha,beta,scintillation_index):
        if self.geometry == 'downlink':
            self.compute_nu_0_DL()
        else:
            self.compute_nu_0_UL()
        I_T = self.mean_intensity * 10 ** (-self.F_T / 10)

        """Andrews, p.514 eq. 74.
        Due to numerical problems (overflow) the logarithmic method is used to calculate the preceding eq."""
        # avg_number_fades = float(2 * np.sqrt( 2 * np.pi * self.alpha * self.beta * self.scintillation_index) * self.nu_0 / (   mp.gamma(self.alpha) * mp.gamma(self.beta)) * (    self.alpha * self.beta * I_T / self.mean_intensity) ** (  (self.alpha + self.beta - 1) / 2) * mp.besselk(self.alpha - self.beta,  2 * mp.sqrt( self.alpha * self.beta * I_T / self.mean_intensity)))
        log_avg_number_fades = mp.log(2) +1/2*mp.log(2)+1/2*mp.log(np.pi)+1/2*mp.log(alpha)+1/2*mp.log(beta)+mp.log(self.nu_0)+1/2*mp.log(scintillation_index)-mp.log(mp.gamma(alpha))-mp.log(mp.gamma(beta))+(alpha+beta-1)/2*mp.log(alpha)+(alpha+beta-1)/2*mp.log(beta)+(alpha+beta-1)/2*mp.log(I_T)-(alpha+beta-1)/2*mp.log(self.mean_intensity)+mp.log(mp.besselk(self.alpha - self.beta,  2 * mp.sqrt( self.alpha * self.beta * I_T / self.mean_intensity)))
        return float(mp.exp(log_avg_number_fades))

    def compute_avg_number_surges(self,alpha,beta,scintillation_index):
        I_T = self.mean_intensity * 10 ** (self.S_T / 10)

        """Andrews, p.514 eq. 74.
        Due to numerical problems (overflow) the logarithmic method is used to calculate the preceding eq."""
        # avg_number_surges =float(2 * np.sqrt(2 * np.pi * self.alpha * self.beta * self.scintillation_index) * self.nu_0 / (mp.gamma(self.alpha) * mp.gamma(self.beta)) * ( self.alpha * self.beta * self.mean_intensity / I_T) ** ((self.alpha + self.beta - 1) / 2) * mp.besselk(self.alpha - self.beta, 2 * mp.sqrt(self.alpha * self.beta * self.mean_intensity / I_T)))
        log_avg_number_surges = mp.log(2) + 1 / 2 * mp.log(2) + 1 / 2 * mp.log(np.pi) + 1 / 2 * mp.log( alpha) + 1 / 2 * mp.log(beta) + mp.log(self.nu_0) + 1 / 2 * mp.log(scintillation_index) - mp.log( mp.gamma(alpha)) - mp.log(mp.gamma(beta)) + (alpha + beta - 1) / 2 * mp.log(alpha) + (alpha + beta - 1) / 2 * mp.log(beta) + (alpha + beta - 1) / 2 * mp.log( I_T) - (alpha + beta - 1) / 2 * mp.log(self.mean_intensity) + mp.log(  mp.besselk(self.alpha - self.beta, 2 * mp.sqrt(self.alpha * self.beta * I_T / self.mean_intensity)))
        return float(mp.exp(log_avg_number_surges))

    def compute_time_parameters(self):
        self.avg_number_fades = self.compute_avg_number_fades(self.alpha,self.beta,self.scintillation_index)
        """Andrews, p.514 eq. 78."""
        self.avg_fade_duration = float(self.fade_probability / self.avg_number_fades)

        if self.scintillation_index_averaged != None:
           self.avg_number_fades_aperture = self.compute_avg_number_fades(self.alpha_avg,self.beta_avg,self.scintillation_index_averaged)
           """Andrews, p.514 eq. 78."""
           self.avg_fade_duration_aperture = float(self.fade_probability_averaged / self.avg_number_fades_aperture)

        if self.scintillation_index_tracked != None:
            self.avg_number_fades_tracked = self.compute_avg_number_fades(self.alpha_tracked, self.beta_tracked,self.scintillation_index_tracked)
            """Andrews, p.514 eq. 78."""
            self.avg_fade_duration_tracked = float(self.fade_probability_tracked / self.avg_number_fades_tracked)

        if self.compute_only_fades != True:

            self.avg_number_surges = self.compute_avg_number_surges(self.alpha, self.beta, self.scintillation_index)
            """Andrews, p.514 eq. 78."""
            self.avg_surges_duration = float(self.surge_probability / self.avg_number_surges)

            if self.scintillation_index_averaged != None:
                self.avg_number_surges_aperture = self.compute_avg_number_surges(self.alpha_avg, self.beta_avg,
                                                                                 self.scintillation_index_averaged)
                """Andrews, p.514 eq. 78."""
                self.avg_surges_duration_aperture = float(  self.surge_probability_averaged / self.avg_number_surges_aperture)

            if self.scintillation_index_tracked != None:
                self.avg_number_surges_tracked = self.compute_avg_number_surges(self.alpha_tracked, self.beta_tracked,
                                                                                self.scintillation_index_tracked)
                """Andrews, p.514 eq. 78."""
                self.avg_surges_duration_tracked = float(
                    self.surge_probability_tracked / self.avg_number_surges_tracked)

