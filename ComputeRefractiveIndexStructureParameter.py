import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt

from LinkedAttribute import *
class ComputeRefractiveIndexStructureParameter(object):
    rms_wind_speed = LinkedAttribute('windSpeedObject')  # Modify this line to access the rms_wind_speed attribute
    def __init__(self, windSpeedObject, daynight_model, wind_height_max, sunset, sunrise, time, hv_ground_cst=1.7e-14):
        self.windSpeedObject = windSpeedObject
        self.daynight_model = daynight_model
        self.hv_ground_cst = hv_ground_cst
        self.h0 = 122
        self.hg = 500
        self.c2n = None
        self.wind_height_max = wind_height_max
        self.wind_speed = windSpeedObject.rms_wind_speed
        self.sunset = sunset
        self.sunrise = sunrise
        self.time = time

    def compute_HV_daytime(self, h):
        # daytime
        # Andrews p.481 eq 1
        # Computes c2n with Hufnagel-Valley model
        term1 = 0.00594 * (self.wind_speed / 27) ** 2 * (10 ** -5 * h) ** 10 * mp.exp(-h / 1000)
        term2 = 2.7 * 10 ** -16 * mp.exp(-h / 1500)
        term3 = self.hv_ground_cst * mp.exp(-h / 100)
        RI = term1 + term2 + term3
        return RI

    def calculate_power_law_parameter(self, TH):
        if 0.75 < TH < 3.5:
            p = -0.11 * (12 - TH) ** 2 + 1.83 * (12 - TH) - 6.22
        elif 3.5 <= TH < 8.5:
            p = 1.45 - 0.02 * (TH - 6) ** 2
        elif 8.5 <= TH < 11.25:
            p = -0.048 * TH ** 2 + 0.68 * TH - 1.06
        else:
            p = np.nan  # If TH is out of bounds, return NaN
        return p

    def compute_HAP_nighttime(self, h, sunset, sunrise, time):
        TP = (self.sunrise - self.sunset) / 12
        TH = (self.time - self.sunrise) / TP
        power_law_parameter = self.calculate_power_law_parameter(TH)
        RI = (1 * (0.00594 * ((self.wind_speed / 27) ** 2) * (((h + self.hg) * 10 ** -5) ** 10) * mp.exp(-(h + self.hg) / 1000))
              + 2.7e-16 * mp.exp(-(h + self.hg) / 1500)
              + 1e-16 * (self.h0 / h) ** power_law_parameter)

        return RI

    def compute_c2n_fct(self):
        if self.daynight_model == 'day':
            self.c2n = self.compute_HV_daytime
            self.c2n_value = float(self.compute_HV_daytime(self.wind_height_max))
        elif self.daynight_model == 'night':
            self.c2n = self.compute_HAP_nighttime
            self.c2n_value = float(self.compute_HAP_nighttime(self.wind_height_max,self.sunset, self.sunrise, self.time))
            print(self.c2n_value)
        else:
            raise ValueError("Invalid c2n_model. Please choose 'day' or 'night'.")
