import numpy as np
import matplotlib.pyplot as plt
from MainFadeSimulation import *

# Elevation and Slew Rate Settings
elvrange_case = 1  # 0 for a single elevation angle, 1 for a range

ElevationAngle = 15
ElevationAngleMin = 15
ElevationAngleMax = 35
ElevationAngleMStep = 5

# Slew Rate Settings
slew_rate = 0.1


# Adjust settings based on elvrange_case
if elvrange_case == 0:
    ElevationAngleMin = ElevationAngle
    ElevationAngleMax = ElevationAngle + ElevationAngleMStep
    printing = True
else:
    ElevationAngleMax = ElevationAngleMax + ElevationAngleMStep
    printing = False

# Adjust settings based on slewrange_case


# Mapping of elevation angles to corresponding slew rates
slew_rate_mapping = {
    15: 0.1,
    20: 0.2,
    25: 0.3,
    30: 0.4,
    35: 0.5
}

# Main loop
geometry = 'uplink'  # 'uplink' or 'downlink'
for elvloop in range(ElevationAngleMin, ElevationAngleMax, ElevationAngleMStep):
    # Get the slew rate for the current elevation angle
    current_slew_rate = slew_rate_mapping.get(elvloop, slew_rate)

    # Run the main fade simulation with the current elevation and slew rate
    mean_fade_loss, fadeLevel_dB, surgeLevel_dB, elevation, rms_wind_speed, Cn2, RL, SI = main_fade_simulation(
        geometry, 'day', 30000, 1.82e-05, 1.2, 1550, 1E-4, 1200, 0,
        elvloop, 5, current_slew_rate, 1.7e-14, 0.004, 1.6, transmission_losses=-2,
        altApertureDiameter=0.02, printResults=printing, C_r=0, integration_step_multiplier=0.1, nr_transmitters=1, compute_only_fades=False
    )

    print('Elevation:', elevation)
    print('wind speed:', rms_wind_speed)
    print('Slew:', current_slew_rate)
    print('Scintillation Index:', SI)
    print('Cn2:', Cn2)
    print('Mean fade loss:', mean_fade_loss)
    print('Fade Level in dB:', fadeLevel_dB)
    print('Surge Level in dB:', surgeLevel_dB)
