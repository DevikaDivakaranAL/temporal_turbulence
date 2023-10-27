# This is a sample Python script.

import numpy as np
import matplotlib.pyplot as plt
from MainFadeSimulation import *

# Set whether to compute for a range of elevation angles
elvrange_case = 0  # 0 for a single elevation angle, 1 for a range

# Elevation angle settings
ElevationAngle = 15
ElevationAngleMin = 15
ElevationAngleMax = 35
ElevationAngleMStep = 5

# Height settings
MinHeight = 0
MaxHeight = 30000

# Adjust settings based on elvrange_case
if elvrange_case == 0:
    ElevationAngleMin = ElevationAngle
    ElevationAngleMax = ElevationAngle + ElevationAngleMStep
    printing = True
else:
    ElevationAngleMax = ElevationAngleMax + ElevationAngleMStep
    printing = False

# Initialize lists to store results
RMS_wind_speeds = []
heights = []

# Main loop
for heightloop in range(MinHeight, MaxHeight, 100):
    for elvloop in range(ElevationAngleMin, ElevationAngleMax, ElevationAngleMStep):
        height, rms_wind_speed, Cn2, RL, SI = main_fade_simulation(
            'uplink', 'day', heightloop, 1.82e-05, 1.2, 1550, 1E-4, 1200, 0,
            elvloop, 8, 0.1, hv_ground_cst=1.7e-14, altApertureDiameter=0.02,
            printResults=printing, C_r=0, integration_step_multiplier=0.1, nr_transmitters=1,
            compute_only_fades=False
        )
        RMS_wind_speeds.append(rms_wind_speed)
        heights.append(height)
        print('Height:', height)
        print('RMS Wind Speed:', rms_wind_speed)
        print('Rytov Variance:', RL)
        print('Scintillation Index:', SI)
        print('Cn2:', Cn2)

# Plotting results
plt.figure(figsize=(10, 5))
plt.plot(heights, RMS_wind_speeds, color='b')
plt.title('RMS Wind Speed vs Height')
plt.xlabel('Height (m)')
plt.ylabel('RMS Wind Speed (m/s)')
plt.grid(True)
plt.show()
