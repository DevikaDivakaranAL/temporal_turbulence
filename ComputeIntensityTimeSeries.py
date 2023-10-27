import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j1
from scipy.integrate import quad
from scipy.fft import ifft
from LinkedAttribute import *
class ComputeIntensityTimeSeries:
    def __init__(self, l0, L0, Tr, V, D, C_n2, frequencies):
        self.l0 = l0
        self.L0 = L0
        self.Tr = Tr
        self.V = V
        self.D = D
        self.C_n2 = C_n2
        self.frequencies = frequencies
        self.k_m = 5.92 / l0
        self.k_0 = 2 * np.pi / L0
        self.We_f_values = None
        self.RescaledIntensity = None
        self.time_values = None

        self.scintillation_index =

    def Ws(self, k):
        return np.sqrt((0.033 * self.C_n2 * np.exp(-k ** 2 / self.k_m ** 2)) / (k ** 2 + self.k_0 ** 2) ** (11 / 6))

    def Wa(self, k):
        return (np.pi * self.D ** 2) / 4 * 2 * j1(np.pi * self.D * k) / (np.pi * self.D * k)

    def We(self, f):
        def integrand(k):
            combined_frequency = np.sqrt(k ** 2 + (f ** 2 / self.V ** 2))
            return self.Wa(combined_frequency) * np.conj(self.Ws(combined_frequency))

        integral, _ = quad(lambda k: integrand(k), 0, np.inf)
        return (self.Tr / self.V) * integral

    def calculate_normalized_intensity(self):
        self.We_f_values = np.array([self.We(f) for f in self.frequencies])
        theta_rand_f = np.random.uniform(0, 2 * np.pi, len(self.frequencies))
        We_prime_f = self.We_f_values * np.exp(1j * theta_rand_f)
        e_t = ifft(We_prime_f)
        self.RescaledIntensity = (e_t / np.mean(e_t)) * 0.093
        time_spacing = 1 / (self.frequencies[-1] - self.frequencies[0])
        num_time_points = len(self.frequencies)
        self.time_values = np.linspace(0, num_time_points * time_spacing, num_time_points, endpoint=False)

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

# Constants
l0 = 0.004
L0 = 1.6
Tr = -2
V = 76
D = 0.32
C_n2 = 2.68e-14
frequencies = np.linspace(1, 10e4, 100000)

# Initialize the simulation
simulation = IntensityFluctuationSimulation(l0, L0, Tr, V, D, C_n2, frequencies)

# Calculate the normalized intensity
simulation.calculate_normalized_intensity()

# Plot the results
simulation.plot_normalized_intensity()
