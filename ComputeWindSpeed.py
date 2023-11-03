import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

class ComputeWindSpeed(object):
    def __init__(self, Vg, Ws, height, geometry):
        self.height = height
        self.Vg = Vg
        self.Ws = Ws
        self.geometry = geometry

        self.u_star = 0.069  # friction velocity, in m/s
        self.z_0 = 0.03  # roughness length, for open terrain, in meters
        self.kappa = 0.41  # von Karman constant, dimensionless
        #self.U_ref = 10  # reference wind speed at reference height, in m/s
        #self.z_ref = 10  # reference height, was 0.01 km, now in meters
        self.alpha = 1 / 7  # power-law exponent, dimensionless


        self.rms_wind_speed = self.compute_wind_speed()

    def Vb(self, height):
        """Computes the wind speed using the Bufton wind model"""
        if height <= 1000:  # logarithmic model for up to 1,000 meters
            VB = (self.u_star / self.kappa) * np.log(height / self.z_0)
        else:  # power law model for above 1,000 meters
            z_ref = height - 100
            U_ref = self.Ws * z_ref + self.Vg + 30 * np.exp(-((z_ref - 9400) / 4800) ** 2)
            VB = U_ref * (height / z_ref) ** self.alpha
        return VB


    def compute_rms_wind_speed_downlink(self):
        """Computes the rms wind speed for downlink"""
        """" Andrews, p.481 eq. 2 """
        h_lower, h_upper = self.height, 1                                # from height of consideration to the ground
        integral, _ = quad(lambda h: self.Vb(h_upper - h) ** 2, h_lower, h_upper)
        return np.sqrt(integral / (h_upper - h_lower))

    def compute_rms_wind_speed_uplink(self):
        """Computes the rms wind speed for uplink"""
        """" Andrews, p.481 eq. 2 """
        h_lower, h_upper = 1, self.height                                # from ground to the height of consideration
        integral, _ = quad(lambda h: self.Vb(h) ** 2, h_lower, h_upper)
        return np.sqrt(integral / (h_upper - h_lower))

    def compute_wind_speed(self):
        if self.geometry == 'uplink':
            return self.compute_rms_wind_speed_uplink()
        elif self.geometry == 'downlink':
            return self.compute_rms_wind_speed_downlink()
        else:
            raise ValueError("Invalid geometry. Please choose 'uplink' or 'downlink'.")
