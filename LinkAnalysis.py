# This is a sample Python script.

import numpy as np
import matplotlib.pyplot as plt
from MainFadeSimulation import *

# Set whether to compute for a range of elevation angles
elvrange_case = 1  # 0 for a single elevation angle, 1 for a range

# Elevation angle settings
ElevationAngle = 15
ElevationAngleMin = 15
ElevationAngleMax = 35
ElevationAngleMStep = 5

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
Scintillation_index = []
Rytov_variance = []

# Height range for wind speed calcu
MinHeight = 122
MaxHeight = 30000
geometry ='uplink'
# Main loop

for elvloop in range(ElevationAngleMin, ElevationAngleMax, ElevationAngleMStep):
    #for height in range(MinHeight,MaxHeight, 1000):
    elevation, rms_wind_speed, Cn2, RL, SI = main_fade_simulation(
        geometry, 'day', 30e3, 1.82e-05, 1.2, 1550, 1E-4, 1200, 0,
        elvloop, 10, 0.1, 1.7e-14, 0.004, 1.6, transmission_losses=-2,
        altApertureDiameter=0.02, printResults=printing, C_r=0, integration_step_multiplier=0.1, nr_transmitters=1, compute_only_fades=False
    )

    print('Elevation:', elevation)
    print('RMS Wind Speed:', rms_wind_speed)
    # print('Rytov Variance:', RL)
    print('Scintillation Index:', SI)
    # print('Cn2:', Cn2)

#
# fig, axs = plt.subplots(3, 1, figsize=(10, 15))
# link_geometry = geometry
# HVmodel = 'night'
# # Plot RMS Wind Speed
# axs[0].plot(heights, RMS_wind_speeds, color='b')
# axs[0].set_title(f' ({link_geometry}, {HVmodel})')
# axs[0].set_xlabel('Height (m)')
# axs[0].set_ylabel('RMS Wind Speed (m/s)')
# axs[0].grid(True)
#
# # Plot Rytov Variance (RL)
# axs[1].plot(heights, Rytov_variance, color='g')
# axs[1].set_xlabel('Height (m)')
# axs[1].set_ylabel('Rytov Variance')
# axs[1].grid(True)
#
# # Plot Scintillation Index (SI)
# axs[2].plot(heights, Scintillation_index, color='r')
# axs[2].set_xlabel('Height (m)')
# axs[2].set_ylabel('Scintillation Index')
# axs[2].grid(True)
#
# # Adjust layout
# plt.tight_layout()
# plt.show()