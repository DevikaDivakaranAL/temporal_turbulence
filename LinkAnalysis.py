import numpy as np
import matplotlib.pyplot as plt
from MainFadeSimulation import *

# Define the range of elevation angles for the simulation
elvrange_case = 1  # 0 for a single elevation angle, 1 for a range

# Elevation angle settings
ElevationAngle = 15  # Default elevation angle
ElevationAngleMin = 15  # Minimum elevation angle for the range
ElevationAngleMax = 35  # Maximum elevation angle for the range
ElevationAngleStep = 5  # Step size for incrementing the elevation angle

# Adjust settings based on the chosen elevation angle case
if elvrange_case == 0:
    # For a single elevation angle
    ElevationAngleMin = ElevationAngle
    ElevationAngleMax = ElevationAngle + ElevationAngleStep
    printing = True
else:
    # For a range of elevation angles
    ElevationAngleMax = ElevationAngleMax + ElevationAngleStep
    printing = False

# Main simulation loop
geometry = 'uplink'  # Choose between 'uplink' and 'downlink'
for elvloop in range(ElevationAngleMin, ElevationAngleMax, ElevationAngleStep):
    # Run the main fade simulation with the specified parameters
    slew, mean_fade_loss, fadeLevel_dB, surgeLevel_dB, elevation, rms_wind_speed, Cn2, RL, SI = main_fade_simulation(
                                                                                                geometry, 'day', 1, 30000, 18.00, 6.00, 4.00,
                                                                                                1.82e-05, 1.2, 1550, 1E-4, 1200, 0, elvloop, 5,
                                                                                                1.7e-14, 0.004, 1.6, transmission_losses=-2, altApertureDiameter=0.02,
                                                                                                printResults=printing, C_r=0, integration_step_multiplier=0.1, nr_transmitters=1, compute_only_fades=False
                                                                                            )

    print('Elevation:', elevation)
    print('wind speed:', rms_wind_speed)
    print('Scintillation Index:', SI)
    print('Cn2:', Cn2)
    print('Mean fade loss:', mean_fade_loss)
    print('Fade Level in dB:', fadeLevel_dB)
    print('Surge Level in dB:', surgeLevel_dB)
