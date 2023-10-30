import numpy as np
import mpmath as mp
from LinkedAttribute import *
class ComputeRefractiveIndexStructureParameter(object):
    rms_wind_speed = LinkedAttribute('windSpeedObject')  # Modify this line to access the rms_wind_speed attribute
    def __init__(self, windSpeedObject, c2n_model, hv_ground_cst=1.7e-14):
        self.c2n_model = c2n_model
        self.hv_ground_cst = hv_ground_cst
        self.h0 = 122
        self.hg =2000
        self.c2n = None
        self.wind_speed = windSpeedObject.rms_wind_speed

    def compute_HV_daytime(self, h):

        RI = (
            0.00594 * (self.wind_speed / 27) ** 2 * (h * 10 ** -5) ** 10 * mp.exp(-h / 1000) +
            2.7e-16 * mp.exp(-h / 1500) + self.hv_ground_cst * mp.exp(-h / 100)
        )
        return RI

    def compute_HAP_nighttime(self, h):

        RI = 1 * (
            0.00594 * (self.wind_speed / 27) ** 2 * ((h + self.hg) / 10 ** 5) ** 10 * mp.exp(-(h + self.hg) / 1000) +
            2.7e-16 * mp.exp(-(h + self.hg) / 1500) + self.hv_ground_cst * (self.h0 / h) ** (4 / 3)
        )
        return RI
        # RI = (0.00594 * (self.wind_speed / 27) ** 2 * (h ** 10) * mp.exp(-h/1000) +
        #       3.02 * 10 ** -17 * mp.exp(-h/1500) +
        #       1.9 * 10 ** -15 * mp.exp(-h/100))
        # return RI

    def compute_c2n_fct(self):
        if self.c2n_model == 'day':
            self.c2n = self.compute_HV_daytime
        elif self.c2n_model == 'night':
            self.c2n = self.compute_HAP_nighttime
        else:
            raise ValueError("Invalid c2n_model")
