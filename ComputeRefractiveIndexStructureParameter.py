import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from LinkedAttribute import *

class ComputeRefractiveIndexStructureParameter(object):
    rms_wind_speed = LinkedAttribute('windSpeedObject')  # Linking to rms_wind_speed attribute from windSpeedObject

    def __init__(self, windSpeedObject, daynight_model, wind_height_max, sunset, sunrise, time, hv_ground_cst=1.7e-14):
        """
        Initializes the ComputeRefractiveIndexStructureParameter class.

        Parameters:
        - windSpeedObject: An instance of a wind speed calculation class.
        - daynight_model (str): 'day' or 'night' to select the appropriate model.
        - wind_height_max (float): Maximum height for wind calculation in m.
        - sunset, sunrise (float): Sunset and sunrise times for nighttime model.
        - time (float): Current time for nighttime model.
        - hv_ground_cst (float): Hufnagel-Valley ground constant.
        """
        self.windSpeedObject = windSpeedObject
        self.daynight_model = daynight_model
        self.hv_ground_cst = hv_ground_cst
        self.h0 = 122  # Reference height parameter
        self.hg = 500  #  height from SL parameter
        self.c2n = None  # Refractive index structure parameter
        self.wind_height_max = wind_height_max
        self.wind_speed = windSpeedObject.rms_wind_speed  # RMS wind speed
        self.sunset = sunset
        self.sunrise = sunrise
        self.time = time

    def compute_HV_daytime(self, h):
        """
        Computes the refractive index structure parameter (Cn²) using the Hufnagel-Valley model for daytime.

        Parameters:
        - h (float): Height in m.

        Returns:
        - RI (float): Refractive index structure parameter at height h.
        """
        term1 = 0.00594 * (self.wind_speed / 27) ** 2 * (10 ** -5 * h) ** 10 * mp.exp(-h / 1000)
        term2 = 2.7 * 10 ** -16 * mp.exp(-h / 1500)
        term3 = self.hv_ground_cst * mp.exp(-h / 100)
        RI = term1 + term2 + term3
        return RI

    def calculate_power_law_parameter(self, TH):
        """
        Calculates a power law parameter based on the time after sunset (TH).

        Parameters:
        - TH (float): Time after sunset normalized to the length of the night.

        Returns:
        - p (float): Power law parameter.
        """
        if 0.75 < TH < 3.5:
            p = -0.11 * (12 - TH) ** 2 + 1.83 * (12 - TH) - 6.22
        elif 3.5 <= TH < 8.5:
            p = 1.45 - 0.02 * (TH - 6) ** 2
        elif 8.5 <= TH < 11.25:
            p = -0.048 * TH ** 2 + 0.68 * TH - 1.06
        else:
            p = np.nan  # Return NaN if TH is out of bounds
        return p

    def compute_HAP_nighttime(self, h, sunset, sunrise, time):
        """
        Computes the refractive index structure parameter (Cn²) using the Hufnagel-Andrews-Philips model for nighttime.

        Parameters:
        - h (float): Height in m.
        - sunset, sunrise, time (float): Environmental time parameters.

        Returns:
        - RI (float): Refractive index structure parameter at height h.
        """
        TP = (sunrise - sunset) / 12  # Normalizing night length to 12 hours
        TH = (time - sunrise) / TP   # Time after sunrise normalized to TP
        power_law_parameter = self.calculate_power_law_parameter(TH)
        RI = (1 * (0.00594 * ((self.wind_speed / 27) ** 2) * (((h + self.hg) * 10 ** -5) ** 10) * mp.exp(-(h + self.hg) / 1000))
              + 2.7e-16 * mp.exp(-(h + self.hg) / 1500)
              + 1e-16 * (self.h0 / h) ** power_law_parameter)
        return RI

    def compute_c2n_fct(self):
        """
        Computes the refractive index structure parameter (Cn²) based on the selected day or night model.

        Sets:
        - self.c2n: Function for Cn² computation.
        - self.c2n_value: Computed Cn² value.
        """
        if self.daynight_model == 'day':
            self.c2n = self.compute_HV_daytime
            self.c2n_value = float(self.compute_HV_daytime(self.wind_height_max))
        elif self.daynight_model == 'night':
            self.c2n = self.compute_HAP_nighttime
            self.c2n_value = float(self.compute_HAP_nighttime(self.wind_height_max, self.sunset, self.sunrise, self.time))
            print(self.c2n_value)
        else:
            raise ValueError("Invalid day/night model. Please choose 'day' or 'night'.")
