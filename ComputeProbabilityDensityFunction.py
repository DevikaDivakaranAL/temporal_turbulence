import numpy as np
import mpmath as mp
from LinkedAttribute import *


class ComputeProbabilityDensityFunction(object):
    """Definition of the attributes that need to be taken from another instance so that they are automatically updated if the instance has them changed.
    For more info read comments in LinkedAttribute class"""
    geometry = LinkedAttribute('turbulenceStrengthObject')
    rytov_variance = LinkedAttribute('turbulenceStrengthObject')
    W_0 = LinkedAttribute('turbulenceStrengthObject')
    W = LinkedAttribute('turbulenceStrengthObject')
    turbulence_strength = LinkedAttribute('turbulenceStrengthObject')
    zenith_angle_rad = LinkedAttribute('turbulenceStrengthObject')
    ALT_altitude = LinkedAttribute('turbulenceStrengthObject')
    SLT_altitude = LinkedAttribute('turbulenceStrengthObject')
    cap_theta = LinkedAttribute('turbulenceStrengthObject')
    r_0 = LinkedAttribute('turbulenceStrengthObject')
    W_eff = LinkedAttribute('beamEffectObject')
    scintillation_index = LinkedAttribute('scintillationIndexObject')
    scintillation_index_averaged = LinkedAttribute('scintillationIndexObject')
    alpha_pe = LinkedAttribute('scintillationIndexObject')
    scintillation_index_tracked =  LinkedAttribute('scintillationIndexObject')
    RescaledIntensity = LinkedAttribute('timeSeriesObject')
    def __init__(self, turbulenceStrengthObject,beamEffectObject,scintillationIndexObject, timeSeriesObject):
        # ACCESSED PARAMETERS
        self.turbulenceStrengthObject = turbulenceStrengthObject                                # ComputeTurbulenceStrength instance
        self.geometry = turbulenceStrengthObject.geometry                                       # string,    link geometry, taken from ComputeTurbulenceStrength instance
        self.rytov_variance = turbulenceStrengthObject.rytov_variance                           # float,     Rytov variance, taken from ComputeTurbulenceStrength instance
        self.W_0 = turbulenceStrengthObject.W_0                                                 # float,     beam waist at transmitter (in m), taken from ComputeTurbulenceStrength instance
        self.W = turbulenceStrengthObject.W                                                     # float,     diffractive beam radius at the receiver (in m),  taken from ComputeTurbulenceStrength instance
        self.turbulence_strength = turbulenceStrengthObject.turbulence_strength                 # string,    turbulence strength, taken from TurbulenceStrength instance
        self.zenith_angle_rad = turbulenceStrengthObject.zenith_angle_rad                       # float,     zenith angle (in rad), taken from ComputeTurbulenceStrength instance
        self.ALT_altitude = turbulenceStrengthObject.ALT_altitude                               # float,     ALT altitude (in m), taken from ComputeTurbulenceStrength instance
        self.SLT_altitude = turbulenceStrengthObject.SLT_altitude                               # float,     SLT altitude (in m), taken from ComputeTurbulenceStrength instance
        self.cap_theta = turbulenceStrengthObject.cap_theta                                     # float,     adimensional beam parameter at the receiver, , taken from ComputeTurbulenceStrength instance
        self.r_0 = turbulenceStrengthObject.r_0                                                 # float,     Fried parameter (in m), taken from ComputeTurbulenceStrength instance

        self.beamEffectObject = beamEffectObject                                                # computeBeamEffects instance
        self.W_eff = beamEffectObject.W_eff                                                     # float,     diffractive beam radius at the receiver (in m), taken from ScintillationIndex instance

        self.scintillationIndexObject = scintillationIndexObject                                # computeScintillationIndex instance
        self.scintillation_index = scintillationIndexObject.scintillation_index                 # float,     scintillation index (untracked in uplink case), taken from computeScintillationIndex instance
        self.scintillation_index_averaged = scintillationIndexObject.scintillation_index_averaged# float,    scintillation index averaged by the aperture, taken from computeScintillationIndex instance
        self.alpha_pe = scintillationIndexObject.alpha_pe                                       # float,     beam wander induced pointing error normalized by the distance, taken from computeScintillationIndex instance
        self.scintillation_index_tracked = scintillationIndexObject.scintillation_index_tracked # float,     scintillation index without beam wander (computed in uplink case), taken from computeScintillationIndex instance

        self.timeSeriesObject = timeSeriesObject
        self.RescaledIntensity = timeSeriesObject.RescaledIntensity
        # COMPUTED PARAMETERS
        self.mean_intensity = None                                                              # mean intensity at the center of the receiver aperture when turbulence is considered
        self.pdf = None                                                                         # function,  probability density as a function of intensity
        self.pdf_averaged = None                                                                # function,  probability density as a function of intensity using averaged scintillation index
        self.pdf_tracked = None                                                                 # function,  probability density as a function of intensity using tracked scintillation index
        self.alpha = None                                                                       # float,     parameter of the gamma-gamma distribution
        self.beta = None                                                                        # float,     parameter of the gamma-gamma distribution
        self.alpha_avg = None                                                                   # float,     parameter of the gamma-gamma distribution computed using averaged scintillation index
        self.beta_avg = None                                                                    # float,     parameter of the gamma-gamma distribution computed using averaged scintillation index
        self.alpha_tracked = None                                                               # float,     parameter of the gamma-gamma distribution computed using tracked scintillation index
        self.beta_tracked = None                                                                # float,     parameter of the gamma-gamma distribution computed using tracked scintillation index
        self.mean_diffraction_only = None                                                       # float,     mean intensity at the center of the receiver aperture considering only diffraction effects (no turbulence)
        self.mean_change_loss = None                                                            # float,     losses due to the movement of the mean intensity when turbulence is considered in the model (losses due to mean decrease)

    def compute_mean_intensity(self):
        """Computes the mean intensity at the center of the receiver aperture (r=0) when turbulence is considered"""
        """" Andrews, p.494 eq. 32 """
        self.mean_intensity = self.W_0**2/self.W_eff**2

        """Computes the mean intensity at the center of the receiver aperture (r=0) when only diffraction is considered"""
        self.mean_diffraction_only = self.W_0**2/self.W**2

        """Computes the losses due to the movement of the mean when turbulence is considered"""
        self.mean_change_loss = 10*np.log10(self.mean_diffraction_only/self.mean_intensity)


    def compute_gammagamma_pdf(self,alpha,beta):
        """Assigns the function to calculate a gamma-gamma to the pdf attribute
        :param alpha: The alpha parameter of a gamma-gamma distribution
        :param beta: The beta parameter of a gamma-gamma distribution
        """

        """" Andrews, p.410 eq. 67 """
        # return lambda i :  2*(alpha*beta)**((alpha+beta)/2)/(mp.gamma(alpha)*mp.gamma(beta)*i)*(i/self.mean_intensity)**((alpha+beta)/2 )*mp.besselk(alpha-beta,2*mp.sqrt(alpha*beta*i/self.mean_intensity))

        """" Due to numerical problems (overflow) the logarithmic method is used to calculate the preceding eq."""
        def log_gamma_gamma(i):
            mp.mp.dps = 50
            return mp.log(2)+(alpha+beta)/2*mp.log(alpha)+(alpha+beta)/2*mp.log(beta)-mp.log(mp.gamma(alpha))-mp.log(mp.gamma(beta))-mp.log(i)+(alpha+beta)/2*mp.log(i)-(alpha+beta)/2*mp.log(self.mean_intensity) + mp.log(mp.besselk(alpha-beta,2*mp.sqrt(alpha*beta*i/self.mean_intensity)))

        return lambda i: mp.exp(log_gamma_gamma(i))

    def compute_pdf(self):
        """Computes the mean intensity and alpha and beta parameters to assign the correct pdf function depending on the configuration."""
        self.compute_mean_intensity()

        if self.geometry == 'downlink':
            """"Andrews, p.511 eq. 68, uses mp.expm1 to reduce truncation error, mp.expm1() = mp.exp() -1 """
            self.alpha = (mp.expm1((0.49 * self.rytov_variance) / (1 + 1.11 * self.rytov_variance ** (6. / 5.)) ** (7. / 6.)))**(-1)
            self.beta =  (mp.expm1((0.51 * self.rytov_variance) / (1 + 0.69 * self.rytov_variance ** (6. / 5.)) ** (5. / 6.)))**(-1)

            if self.scintillation_index_averaged != None:
                """" In the downlink and weak case we can also compute the scintillation index averaged over the aperture. 
                Andrews, p.513 middle paragraph."""
                self.alpha_avg = (0.49 * self.scintillation_index_averaged)**(-1)
                self.beta_avg = (0.51 * self.scintillation_index_averaged) ** (-1)
                self.pdf_averaged = self.compute_gammagamma_pdf(self.alpha_avg,self.beta_avg)

        else:
            """" For the uplink case: Andrews, p.517 eq. 80 and eq. 81 """
            sigma2_x =  5.95*(self.SLT_altitude-self.ALT_altitude)**2*mp.sec(self.zenith_angle_rad)**2*((2*self.W_0)/self.r_0)**(5/3)*(self.alpha_pe/self.W)**2 + mp.expm1((0.49*self.rytov_variance)/(1+(1 + self.cap_theta)*0.56*self.rytov_variance**(12./10.))**(7./6.) )
            sigma2_y = mp.expm1((0.51*self.rytov_variance)/(1+0.69*self.rytov_variance**(12./10.))**(5./6.))
            self.alpha = 1/sigma2_x
            self.beta = 1/sigma2_y

            sigma2_x_tracked =  mp.expm1((0.49*self.rytov_variance)/(1+(1 + self.cap_theta)*0.56*self.rytov_variance**(12./10.))**(7./6.) )
            sigma2_y_tracked = mp.expm1((0.51*self.rytov_variance)/(1+0.69*self.rytov_variance**(12./10.))**(5./6.))
            self.alpha_tracked = 1 / sigma2_x_tracked
            self.beta_tracked = 1 / sigma2_y_tracked
            self.pdf_tracked = self.compute_gammagamma_pdf(self.alpha_tracked,self.beta_tracked)

        self.pdf = self.compute_gammagamma_pdf(self.alpha, self.beta)
