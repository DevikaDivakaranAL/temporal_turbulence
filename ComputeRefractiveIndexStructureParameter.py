import numpy as np
import mpmath as mp
from LinkedAttribute import *
class ComputeRefractiveIndexStructureParameter(object):
    rms_wind_speed = LinkedAttribute('WindSpeedModel')
    def __init__(self, rms_wind_speed, c2n_model='day', hv_ground_cst=1.7e-14):
        self.c2n_model = c2n_model
        self.hv_ground_cst = hv_ground_cst
        self.wind_speed = rms_wind_speed
        self.h0 = mp.mpf(1000)
        self.hg = mp.mpf(100)
        self.c2n = None
        self.compute_c2n_fct()

    def compute_HV_daytime(self, h):
        h = mp.mpf(h)
        refractive_index_structure = (
            0.00594 * (self.wind_speed / 27) ** 2 * (h * 10 ** -5) ** 10 * mp.exp(-h / 1000) +
            2.7e-16 * mp.exp(-h / 1500) + self.hv_ground_cst * mp.exp(-h / 100)
        )
        return refractive_index_structure

    def compute_HV_nighttime(self, h):
        h = mp.mpf(h)
        refractive_index_structure = 0.1 * (
            0.00594 * (self.wind_speed / 27) ** 2 * ((h + self.hg) / 10 ** 5) ** 10 * mp.exp(-(h + self.hg) / 1000) +
            2.7e-16 * mp.exp(-(h + self.hg) / 1500) + self.hv_ground_cst * (self.h0 / h) ** (4 / 3)
        )
        return refractive_index_structure

    def compute_c2n_fct(self):
        if self.c2n_model == 'day':
            self.c2n = self.compute_HV_daytime
        elif self.c2n_model == 'night':
            self.c2n = self.compute_HV_nighttime
        else:
            raise ValueError("Invalid c2n_model. Please choose 'day' or 'night'.")
