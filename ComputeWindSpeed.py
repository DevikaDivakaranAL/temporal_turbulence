import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

class WindSpeedModel(object):
    def __init__(self, Vg=8, Ws=0.5):
        """
        Initialize the wind speed model.

        :param Vg: Ground wind speed in m/s, default is 8 m/s
        :param Ws: Slew rate in deg/s, default is 0.5 deg/s
        """
        self.Vg = Vg
        self.Ws = Ws
        self.rms_wind_speed = None
        self.compute_rms_wind_speed()

    def Vb(self, h):
        """
        Bufton wind model function.

        :param h: Altitude in meters
        :return: Wind speed at altitude h
        """
        return (self.Ws * h + self.Vg + 30 * np.exp(-((h - 9400) / 4800) ** 2))

    def compute_rms_wind_speed(self):
        """
        Calculate the rms wind speed.
        """
        h_lower = 5e3  # Lower limit of integration in meters
        h_upper = 20e3  # Upper limit of integration in meters
        integral, _ = quad(lambda h: self.Vb(h) ** 2, h_lower, h_upper)
        self.rms_wind_speed = np.sqrt(integral / (h_upper - h_lower))

    def plot_wind_speeds(self, altitudes):
        """
        Plot wind speeds as a function of altitude.

        :param altitudes: Array of altitudes in meters
        """
        Vb_wind_speeds = [self.Vb(altitude) for altitude in altitudes]
        RMS_wind_speeds = [self.rms_wind_speed for _ in altitudes]

        plt.figure(figsize=(10, 6))
        plt.plot(altitudes / 1000, Vb_wind_speeds, label='Vb Wind Speed')
        plt.plot(altitudes / 1000, RMS_wind_speeds, label='RMS Wind Speed')
        plt.title('Wind Speed vs. Altitude')
        plt.xlabel('Altitude (km)')
        plt.ylabel('Wind Speed (m/s)')
        plt.legend()
        plt.grid(True)
        plt.show()


# Usage
wind_model = WindSpeedModel(Vg=8, Ws=0.5)
altitudes = np.linspace(0, 20000, 200)
wind_model.plot_wind_speeds(altitudes)
