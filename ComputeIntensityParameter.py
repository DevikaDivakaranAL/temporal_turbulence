import numpy as np
import mpmath as mp
from LinkedAttribute import *
from scipy.integrate import quad
from scipy.optimize import fsolve


class ComputeIntensityParameters(object):
    """Definition of the attributes that need to be taken from another instance so that they are automatically updated if the instance has them changed.
    For more info read comments in LinkedAttribute class"""
    elevation = LinkedAttribute('turbulenceStrengthObject')
    geometry = LinkedAttribute('turbulenceStrengthObject')
    scintillation_index_averaged = LinkedAttribute('scintillationIndexObject')
    scintillation_index_tracked = LinkedAttribute('scintillationIndexObject')
    mean_intensity = LinkedAttribute('pdfObject')
    pdf = LinkedAttribute('pdfObject')
    pdf_averaged = LinkedAttribute('pdfObject')
    pdf_tracked= LinkedAttribute('pdfObject')

    def __init__(self, turbulenceStrengthObject,scintillationIndexObject,pdfObject, probability =1e-2,integration_step_multiplier =1,nr_transmitters =1,compute_only_fades =True):
        #ACCESSED PARAMETERS
        self.turbulenceStrengthObject= turbulenceStrengthObject                         # ComputeTurbulenceStrength instance
        self.elevation = turbulenceStrengthObject.elevation                             # float,     elevation angle (in degrees), taken from ComputeTurbulenceStrength instance
        self.geometry = turbulenceStrengthObject.geometry                               # string,    link geometry, taken from ComputeTurbulenceStrength instance

        self.scintillationIndexObject =scintillationIndexObject                         # ComputeScintillationIndex instance
        self.scintillation_index_averaged = scintillationIndexObject.scintillation_index_averaged # float, scintillation index averaged by the aperture, taken from ComputeScintillationIndex instance
        self.scintillation_index_tracked = scintillationIndexObject.scintillation_index_tracked # float, scintillation index without beam wander (computed in uplink case), taken from ComputeScintillationIndex instance

        self.pdfObject = pdfObject                          # ComputeProbabilityDensityFunction instance
        self.mean_intensity = pdfObject.mean_intensity      # float,     mean intensity at the center of the receiver aperture, taken from ComputeProbabilityDensityFunction instance
        self.pdf = pdfObject.pdf                            # function,  probability density as a function of intensity, taken from ComputeProbabilityDensityFunction instance
        self.pdf_averaged = pdfObject.pdf_averaged          # function,  probability density as a function of intensity using averaged scintillation index, taken from ComputeProbabilityDensityFunction instance
        self.pdf_tracked = pdfObject.pdf_tracked            # function,  probability density as a function of intensity using tracked scintillation index, taken from ComputeProbabilityDensityFunction instance


        # INPUT PARAMETERS
        self.probability = probability                      # probability of interest. We want to know the fade level associated to this probability
        self.integration_step_multiplier = integration_step_multiplier # factor that multiplies the integration step. Can be increased to make computation faster (or decreased for more precision), by default = 1.
        self.nr_transmitters = nr_transmitters              # Number of independent trasmitters N, by default = 1.
        self.compute_only_fades = compute_only_fades        # Parameter to determine if we want to compute only fade statistics or surges too.

        # COMPUTED PARAMETERS
        self.F_T = None                                     # Fade threshold associated to the probability of interest.
        self.F_T_tracked = None                             # Fade threshold associated to the probability of interest for the tracked case
        self.F_T_averaged = None                            # Fade threshold associated to the probability of interest for the aperture averaged scint index case
        self.F_T_multiple = None                            # Fade threshold associated to the probability of interest for N transmitters

        self.S_T = None                                    # Surge threshold associated to the probability of interest.
        self.S_T_tracked = None                            # Surge threshold associated to the probability of interest for the tracked case
        self.S_T_averaged = None                           # Surge threshold associated to the probability of interest for the aperture averaged scint index case
        self.S_T_multiple = None                           # Surge threshold associated to the probability of interest for N transmitters

        self.fade_probability = None                       # computed probability associated to F_T, should be very close to the probability of interest (check parameter).
        self.fade_probability_tracked = None               # computed probability associated to F_T_tracked, should be very close to the probability of interest (check parameter).
        self.fade_probability_averaged = None              # computed probability associated to F_T_averaged, should be very close to the probability of interest (check parameter).
        self.fade_probability_multiple = None              # computed probability associated to F_T_multiple, should be very close to the probability of interest (check parameter).

        self.surge_probability = None                      # computed probability associated to S_T, should be very close to the probability of interest (check parameter).
        self.surge_probability_averaged = None             # computed probability associated to S_T_tracked, should be very close to the probability of interest (check parameter).
        self.surge_probability_tracked = None              # computed probability associated to S_T_averaged, should be very close to the probability of interest (check parameter).
        self.surge_probability_multiple = None             # computed probability associated to S_T_multiple, should be very close to the probability of interest (check parameter).


    def converge_to_fade_intensity_of_interest(self, pdf):
        """calculates the intensity level (I_T_fade) where the cumulative distribution function (CDF) of a given probability density function (pdf) is closest to a specified probability of interest."""

        di = self.mean_intensity * self.probability / 10 *self.elevation *self.integration_step_multiplier #choose an integration step proportional to the elevation and the probability of interest
        I_T_fade = 0
        p = 0


        while p <= self.probability:  # iterate until the integral of pdf is close to the probability of interest
            I_T_fade += di
            p += (pdf(I_T_fade + di) + pdf(I_T_fade)) /2* di  # uses a numerical integration approach (Trapezoidal rule)
            # to integrate the pdf and find the corresponding intensity level.

        return I_T_fade, p

    def converge_to_surge_intensity_of_interest(self, pdf):
        """Given a pdf it computes the intensity S_T for which the cumulative distribution function is the closest to (1 - probability) using optimization methods"""

        def objective_func_surge(I_T_surge):
            """Define the objective function to minimize, in this case: int_0^inf(pdf(i))di - (1-probability)"""
            integral = quad(self.pdf, 0, I_T_surge)[0]
            return (1 - self.probability) - integral

        vfunc_surge = np.vectorize(objective_func_surge)                      # transforms functions which are not numpy-aware into functions that can operate on (and return) numpy arrays.
        I_T_surge_guess = self.mean_intensity                                 # Intensity guess used to start the iterations to minimize the objective function
        I_T_surge = float(fsolve(vfunc_surge, I_T_surge_guess)[0])            # Solves for the minimum (or zero) of the objective function

        return I_T_surge, (1 - quad(pdf, 0, I_T_surge)[0])

    def find_nearest(self, array, value):
        """
        From https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
        :param array: a list of 1d array of value
        :param value: the value to approach
        :return: the closest value in array
        """
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx], idx

    def compute_multiple_transmitters(self,nr_transmitters, pdf, fadeProb):
        """Given a pdf it computes N times the convolution with itself (pdf of the sum of N i.i.d. variables) and returns fade and surge statistics on the resulting convolved pdf."""

        pmf_array = []  # array of probability mass functions
        i_initial = self.mean_intensity * fadeProb  # smallest intensity for which the pdf is evaluated
        i_final = 4 * self.mean_intensity           # biggest intensity for which the pdf is evaluated
        n_points = 1000
        di = (i_final - i_initial) / n_points
        i_range = np.arange(i_initial, i_final, di) # array of intensity values for which the pdf is evaluated

        for i in i_range:
            pmf_array.append(pdf(i) * di)           # discretize the pdf in pmf

        pmf = pmf_array                             # define the pmf to convolve with itself N times.
        convolved_i_initial = i_initial             # smallest intensity for which the convolved pdf will be evaluated
        convolved_i_final = i_final                 # biggest intensity for which the convolved pdf will be evaluated

        for i in range(1, nr_transmitters):
            pmf = np.convolve(pmf_array, pmf)       #convolve pmf with itself N times.
            convolved_i_initial += i_initial        # smallest intensity for which the convolved pdf will be evaluated is N*i_initial
            convolved_i_final += i_final - di       # biggest intensity for each convolution of functions of length M is 2*M - 1

        pdf_cont = pmf / di                         # make the convolved pmf a pdf again

        convolved_i_range = np.arange(convolved_i_initial, convolved_i_final, di) # array of intensity values for which the convolved pdf is evaluated

        """Nestor D. Chatzidiamantis, On the Distribution of the Sum of Gamma-Gamma Variates and Applications in RF and Optical Wireless Communications"""
        convolved_mean_intensity = nr_transmitters * self.mean_intensity



        i_fade_range = [i for i in convolved_i_range if i <= convolved_mean_intensity] # define range of intensities below the mean (fade range)
        proba_fade_vector = []                                                          # define the vector of probabilities that will be associated with each intensity drop below the mean

        i_surge_range = [i for i in convolved_i_range if i > convolved_mean_intensity] # define range of intensities below the mean (fade range)
        proba_surge_vector = []                                                        # define the vector of probabilities that will be associated with each intensity drop below the mean


        for i in range(len(i_fade_range)):                                              # compute the vector of probabilities that will be associated with each intensity drop below the mean
            proba_fade_i = np.trapz(pdf_cont[:i], i_fade_range[:i])
            proba_fade_vector.append(proba_fade_i)

        p_mean = proba_fade_vector[-1]                                                  # probability associated to mean intensity corresponds to the last entry of the vector of fade probabilities

        for i in range(len(i_surge_range)):
            proba_surge_i = np.trapz(pdf_cont[:i], i_surge_range[:i])                   # compute the vector of probabilities that will be associated with each intensity above the mean (surge)
            proba_surge_vector.append(proba_surge_i)

        p_fade, index_fade = self.find_nearest(proba_fade_vector, fadeProb)             # find the probability of a fade that has the closest value to the probability of interest and the correspective index
        p_surge, index_surge = self.find_nearest(proba_surge_vector, (1-fadeProb))      # find the probability of a surge that has the closest value to the probability of interest and the correspective index

        I_T = float(convolved_i_range[index_fade])                                      #via the index determine the fade intensity associated to the closest probability to the probability of interest
        F_T = float(10 * np.log10(convolved_mean_intensity / I_T))                      #compute fade threshold in dB

        I_T_surge = float(convolved_i_range[index_surge])                               #via the index determine the surge intensity associated to the closest probability to the probability of interest
        S_T = float(10 * np.log10(I_T_surge/convolved_mean_intensity))                #compute surge threshold in dB

        return F_T, float(p_fade), S_T, float(p_surge)

    def compute_intensity_parameters(self):
        """calls the different functions to compute the fade/surge threshold associated to the probability of interest."""
        I_T_fade, fade_proba = self.converge_to_fade_intensity_of_interest(self.pdf)
        self.F_T = 10 * np.log10(self.mean_intensity / I_T_fade)
        self.fade_probability = float(fade_proba)

        if self.geometry == 'uplink' and self.scintillation_index_tracked != None:
                """In case of uplink, it computes the fade statistics with no beam wander"""
                I_T_fade_tracked, fade_proba_tracked = self.converge_to_fade_intensity_of_interest(self.pdf_tracked)

                self.F_T_tracked = 10 * np.log10(self.mean_intensity / I_T_fade_tracked)
                self.fade_probability_tracked= float(fade_proba_tracked)
                # print('proof: p_fade_tracked', mp.quad(self.pdf_tracked, [0, I_T_fade_tracked]))

        if self.geometry == 'downlink' and self.scintillation_index_averaged != None:
            """In case of downlink and weak turbulence, it computes the fade statistics normalized by the receiver aperture"""
            I_T_fade_avg, fade_proba_avg = self.converge_to_fade_intensity_of_interest(self.pdf_averaged)
            self.F_T_averaged = 10 * np.log10(self.mean_intensity / I_T_fade_avg)
            self.fade_probability_averaged = float(fade_proba_avg)
            # print('proof: p_avg', mp.quad(self.pdf_averaged, [0, I_T_fade_avg]))

        if self.compute_only_fades != True:
            I_T_surge, surge_proba = self.converge_to_surge_intensity_of_interest(self.pdf)
            self.S_T = - 10 * np.log10(self.mean_intensity / I_T_surge)
            self.surge_probability = float(surge_proba)

            if self.geometry == 'uplink' and self.scintillation_index_tracked != None:
                """In case of uplink, it computes the surge statistics with no beam wander"""

                I_T_surge_tracked, surge_proba_tracked = self.converge_to_surge_intensity_of_interest(self.pdf_tracked)

                self.S_T_tracked = -10 * np.log10(self.mean_intensity / I_T_surge_tracked)
                self.surge_probability_tracked = float(surge_proba_tracked)
                # print('proof: p_surge_tracked', mp.quad(self.pdf_tracked, [0, I_T_surge_tracked]))

            if self.geometry == 'downlink' and self.scintillation_index_averaged != None:
                """In case of downlink and weak turbulence, it computes the surge statistics normalized by the receiver aperture"""

                I_T_surge_avg, surge_proba_avg = self.converge_to_surge_intensity_of_interest(self.pdf_averaged)
                self.S_T_averaged = -10 * np.log10(self.mean_intensity / I_T_surge_avg)
                self.surge_probability_averaged = float(surge_proba_avg)
                # print('proof: p_surge_averaged', mp.quad(self.pdf_averaged, [0, I_T_surge_avg]))

        if self.nr_transmitters > 1:
            fade_multiple, multiple_fade_proba, surge_multiple, multiple_surge_proba = self.compute_multiple_transmitters(self.nr_transmitters, self.pdf, self.fade_probability)
            self.F_T_multiple = fade_multiple
            self.fade_probability_multiple = multiple_fade_proba
            self.S_T_multiple = surge_multiple
            self.surge_probability_multiple = multiple_surge_proba
            mp.mp.dps = 50

        proba_calculated_mpmath = mp.quad(self.pdf, [0, I_T_fade])
        #rel_probability_error = abs((proba_calculated_mpmath - fade_proba)/proba_calculated_mpmath)
        # print('proof: p', mp.quad(self.pdf, [0, I_T_fade]))
        # print('relative error', rel_probability_error)
        # print('proof: p_surge', mp.quad(self.pdf, [0, I_T_surge]))
