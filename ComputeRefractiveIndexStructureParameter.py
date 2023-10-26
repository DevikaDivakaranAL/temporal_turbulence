import numpy as np
import mpmath as mp

class RefractiveIndexDayNight(object):
    """Class to compute the refractive index structure parameter (Cn^2) using different models."""
    def __init__(self, c2n_model='day', wind_speed=None, hv_ground_cst=1.7e-14):
        self.c2n_model = c2n_model
        self.hv_ground_cst = hv_ground_cst
        self.wind_speed = wind_speed
        self.h0 = 1000
        self.hg = 100
        self.c2n = None
        self.compute_c2n_fct()

    def compute_HV_daytime(self, h):
        """Computes Cn^2 using the Hufnagel-Valley model for daytime conditions."""
        refractive_index_structure = (
            0.00594 * (self.wind_speed / 27.) ** 2. * (h * 10 ** -5.) ** 10. * mp.exp(-h / 1000.) +
            2.7e-16 * mp.exp(-h / 1500.) + self.hv_ground_cst * mp.exp(-h / 100.)
        )
        return refractive_index_structure

    def compute_HV_nighttime(self, h):
        """Computes Cn^2 using the Hufnagel-Andrews-Phillips model for nighttime conditions."""
        refractive_index_structure = 0.1 * (
            0.00594 * (self.wind_speed / 27.) ** 2. * ((h + self.hg) / 10 ** 5.) ** 10. * mp.exp(-(h + self.hg) / 1000.) +
            2.7e-16 * mp.exp(-(h + self.hg) / 1500.) + self.hv_ground_cst * (self.h0 / h) ** 4 / 3
        )
        return refractive_index_structure

    def compute_c2n_fct(self):
        """Assigns the appropriate Cn^2 computation function based on the selected model."""
        if self.c2n_model == 'day':
            self.c2n = self.compute_HV_daytime
        elif self.c2n_model == 'night':
            self.c2n = self.compute_HV_nighttime
        else:
            raise ValueError("Invalid c2n_model. Please choose 'day' or 'night'.")
