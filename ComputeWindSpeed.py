import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

class ComputeWindSpeed(object):
    def __init__(self, Vg, Ws, geometry='uplink'):
        self.Vg = Vg
        self.Ws = Ws
        self.geometry = geometry
        self.rms_wind_speed = None
        self.compute_wind_speed()
    def Vb(self, h):
        return (self.Ws * h + self.Vg + 30 * np.exp(-((h - 9400) / 4800) ** 2))

    def compute_rms_wind_speed_downlink(self):
        h_lower = 5e3
        h_upper = 20e3
        integral, _ = quad(lambda h: self.Vb(h) ** 2, h_lower, h_upper)
        self.rms_wind_speed = np.sqrt(integral / (h_upper - h_lower))

    def compute_rms_wind_speed_uplink(self):
        h_lower = 5e3
        h_upper = 20e3
        integral, _ = quad(lambda h: self.Vb(h) ** 2, h_lower, h_upper)
        self.rms_wind_speed = np.sqrt(integral / (h_upper - h_lower))

    def compute_wind_speed(self):
        if self.geometry == 'uplink':
            self.compute_rms_wind_speed_uplink()
        elif self.geometry == 'downlink':
            self.compute_rms_wind_speed_downlink()
        else:
            raise ValueError("Invalid geometry. Please choose 'uplink' or 'downlink'.")
