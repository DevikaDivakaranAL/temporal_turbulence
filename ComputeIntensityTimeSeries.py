import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j1
from scipy.integrate import quad
from scipy.fft import ifft
from LinkedAttribute import *
class ComputeIntensityTimeSeries:
    c2n_value = LinkedAttribute('refractiveIndexObject')
    scintillation_index_tracked = LinkedAttribute('scintillationIndexObject')
    rms_wind_speed = LinkedAttribute('windSpeedObject')
    wavelength = LinkedAttribute('turbulenceStrengthObject')
    def __init__(self,refractiveIndexObject, scintillationIndexObject, windSpeedObject, turbulenceStrengthObject,
                 lower_limit,upper_limit, Transmission_losses, aperture_diameter):
        #accessed parameters
        self.refractiveIndexObject = refractiveIndexObject  # ComputeRefractiveIndexStructure instance
        self.c2n_value = refractiveIndexObject.c2n_value
        self.scintillationIndexObject = scintillationIndexObject  # ComputeScintillationIndex instance
        self.scintillation_index = scintillationIndexObject.scintillation_index
        self.windSpeedObject = windSpeedObject  # ComputeRefractiveIndexStructure instance
        self.V = windSpeedObject.rms_wind_speed
        self.turbulenceStrengthObject = turbulenceStrengthObject  # ComputeRefractiveIndexStructure instance
        self.wavelength = turbulenceStrengthObject.wavelength
        #input parameters
        self.l0 = lower_limit
        self.L0 = upper_limit
        self.Tr = Transmission_losses
        self.D = aperture_diameter
        #computed parameters
        self.We_f_values = None
        self.RescaledIntensity = None
        self.time_values = None

        self.k_m = 5.92 / self.l0
        self.k_0 = 2 * np.pi / self.L0
       # self.frequency = 2.99e8 / self.wavelength
        self.frequencies = np.linspace(1, 10e4, 100000)
    def Ws(self, k):
        """power spectral component described by the von Kármán spectrum"""
        """" M.Toyoshima P.4 eq. 11"""
        Ws_k = np.sqrt((0.033 * self.c2n_value * np.exp(-k ** 2 / self.k_m ** 2)) / (k ** 2 + self.k_0 ** 2) ** (11 / 6))
        return Ws_k
    def Wa(self, k):
        """spatial spectrum"""
        """" M.Toyoshima P.4 eq. 13"""
        Wa_k= (np.pi * self.D ** 2) / 4 * 2 * j1(np.pi * self.D * k) / (np.pi * self.D * k)
        return Wa_k

    def We(self, f):
        """temporal power spectrum"""
        """" M.Toyoshima P.4 eq. 10"""

        def integrand(k):
            combined_frequency = np.sqrt(k ** 2 + (f ** 2 / self.V ** 2))
            return self.Wa(combined_frequency) * np.conj(self.Ws(combined_frequency))

        integral, _ = quad(lambda k: (integrand(k)), 0, np.inf)
        We_f = (self.Tr / self.V) * integral
        return We_f
    def plot_normalized_intensity(self):
        if self.RescaledIntensity is None or self.time_values is None:
            print("Please run the 'calculate_normalized_intensity' method first.")
            return
        plt.figure(figsize=(10, 5))
        plt.plot(self.time_values, np.abs(self.RescaledIntensity))
        plt.title(f'Normalized Intensity vs Time [rms wind speed : {self.V}, Aperture Diameter: {self.D}]')
        plt.xlabel('Time (s)')
        plt.ylabel('Normalized Intensity')
        plt.grid(True)
        plt.show()
    def compute_intensity_time_series(self):

        self.We_f_values = np.array([self.We(f) for f in self.frequencies])
        # print(self.We_f_values)
        theta_rand_f = np.random.uniform(0, 2 * np.pi, len(self.frequencies))   #random phase perturbation from 0 to 2π

        """" M.Toyoshima P.5 eq. 17 """

        We_prime_f = self.We_f_values * np.exp(1j * theta_rand_f)                    #the power spectral density (PSD) for generating the random time-varying signals

        """" M.Toyoshima P.5 eq. 18"""

        e_t = ifft(We_prime_f)                                                       #the time-varying electrical signal

        self.RescaledIntensity = e_t / np.mean(e_t)* 1/self.scintillation_index
        time_spacing = 1 / (self.frequencies[-1] - self.frequencies[0])
        num_time_points = len(self.frequencies)
        self.time_values = np.linspace(0, num_time_points * time_spacing, num_time_points, endpoint=False)
        self.plot_normalized_intensity()





