import numpy as np
from scipy.integrate import quad
from math import sqrt, sin, cos, radians, asin, degrees

class ComputeWindSpeed(object):
    def __init__(self, Vg, wind_height_min, wind_height_max, elevation, SLT_altitude, geometry):
        """
        Initializes the ComputeWindSpeed class.

        Parameters:
        - Vg (float): Ground wind speed in m/s.
        - wind_height_min (float): Minimum height for wind calculation in m.
        - wind_height_max (float): Maximum height for wind calculation in m.
        - elevation (float): Elevation angle in degrees.
        - SLT_altitude (float): Satellite altitude in km.
        - geometry (str): 'uplink' or 'downlink'.
        """
        self.wind_height_min = wind_height_min
        self.wind_height_max = wind_height_max
        self.Vg = Vg
        self.geometry = geometry
        self.elevation = elevation
        self.SLT_altitude = SLT_altitude
        self.rms_wind_speed = self.compute_wind_speed()
        self.slew = None
    def calculate_slew(self):
        """
        Calculates the slew rate based on satellite parameters and elevation angle.

        Returns:
        - slew (float): Slew rate in degrees/s.
        """
        gravitational_constant = 0.0000000000667
        earth_mass = 5.972E+24  # Mass of the Earth in kg
        radius_earth = 6371  # Radius of Earth in km

        # Calculate orbital velocity
        r_1 = self.SLT_altitude + radius_earth  # Orbital radius in km
        orbital_velocity = sqrt((gravitational_constant * earth_mass) / (r_1 * 1000)) / 1000  # in km/s

        # Calculate angle ψ and Earth Central angle
        psi = degrees(asin(radius_earth / r_1 * sin(radians(self.elevation + 90))))
        EC = 90 - self.elevation - psi  # Earth Central angle

        # Determine arc length and tangential velocity
        L = radius_earth * sin(radians(EC)) / sin(radians(psi))  # Arc length in km
        v_t = orbital_velocity * cos(radians(psi))  # Tangential velocity in km/s

        # Calculate slew rate
        slew = degrees(v_t / L)  # Slew rate in degrees/s
        return slew
    def Vb(self, height, Vg):
        """
        Calculates the Bufton wind speed at a given height.

        Parameters:
        - height (float): Height in m.
        - Vg (float): Ground wind speed in m/s.

        Returns:
        - VB (float): Bulk wind speed at the given height.
        """
        self.slew = self.calculate_slew()
        VB = self.slew * height + Vg + 30 * np.exp(-((height - 9400) / 4800) ** 2)
        return VB

    def compute_rms_wind_speed(self, h_lower, h_upper):
        """
        Computes the RMS wind speed over a given height range.

        Parameters:
        - h_lower (float): Lower height in m.
        - h_upper (float): Upper height in m.

        Returns:
        - (float): RMS wind speed over the given height range.
        """
        integral, _ = quad(lambda h: self.Vb(h_upper - h, self.Vg) ** 2, h_lower, h_upper)
        return sqrt(integral / (h_upper - h_lower))

    def compute_wind_speed(self):
        """
        Determines the height range for wind speed computation based on geometry and computes RMS wind speed.

        Returns:
        - (float): Computed RMS wind speed.
        """
        if self.geometry == 'uplink':
            h_lower = self.wind_height_min
            h_upper = self.wind_height_max
        elif self.geometry == 'downlink':
            h_lower = self.wind_height_max
            h_upper = self.wind_height_min
        else:
            raise ValueError("Invalid geometry. Please choose 'uplink' or 'downlink'.")

        return self.compute_rms_wind_speed(h_lower, h_upper)
