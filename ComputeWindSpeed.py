import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

class ComputeWindSpeed(object):
    def __init__(self, Vg, Ws, height, geometry='uplink'):
        self.height = height
        self.Vg = Vg
        self.Ws = Ws
        self.geometry = geometry
        self.rms_wind_speed = self.compute_wind_speed()

    def Vb(self, height, Vg, Ws):
        VB = Ws * height + Vg + 30 * np.exp(-((height - 9400) / 4800) ** 2)
        return VB

    def compute_rms_wind_speed_downlink(self):
        h_lower = 1
        h_upper = self.height
        integral, _ = quad(lambda h: self.Vb(h, self.Vg, self.Ws) ** 2, h_lower, h_upper)
        return np.sqrt(integral / (h_upper - h_lower))

    def compute_rms_wind_speed_uplink(self):
        h_lower = 1
        h_upper = self.height
        integral, _ = quad(lambda h: self.Vb(h, self.Vg, self.Ws) ** 2, h_lower, h_upper)
        return np.sqrt(integral / (h_upper - h_lower))

    def compute_wind_speed(self):
        if self.geometry == 'uplink':
            return self.compute_rms_wind_speed_uplink()
        elif self.geometry == 'downlink':
            return self.compute_rms_wind_speed_downlink()
        else:
            raise ValueError("Invalid geometry. Please choose 'uplink' or 'downlink'.")
