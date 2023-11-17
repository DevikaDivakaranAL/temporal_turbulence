import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

class ComputeWindSpeed(object):
    def __init__(self, Vg, slew, wind_height_min, wind_height_max, geometry):
        self.wind_height_min = wind_height_min
        self.wind_height_max = wind_height_max
        self.Vg = Vg
        self.slew = slew
        self.geometry = geometry
        self.rms_wind_speed = self.compute_wind_speed()

    def Vb(self, height, Vg, slew):
        #bufton wind model with slew rate
        VB = slew * height + Vg + 30 * np.exp(-((height - 9400) / 4800) ** 2)
        return VB
    def compute_rms_wind_speed(self, h_lower, h_upper):
        integral, _ = quad(lambda h: self.Vb(h_upper - h, self.Vg, self.slew) ** 2, h_lower, h_upper)
        return np.sqrt(integral / (h_upper - h_lower))


    def compute_wind_speed(self):
        if self.geometry == 'uplink':
            h_lower = self.wind_height_min
            h_upper = self.wind_height_max
            return self.compute_rms_wind_speed(h_lower, h_upper)
        elif self.geometry == 'downlink':
            h_lower = self.wind_height_max
            h_upper = self.wind_height_min
            return self.compute_rms_wind_speed(h_lower, h_upper)
        else:
            raise ValueError("Invalid geometry. Please choose 'uplink' or 'downlink'.")
