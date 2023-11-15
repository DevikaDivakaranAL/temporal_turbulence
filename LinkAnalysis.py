import numpy as np
import matplotlib.pyplot as plt
from MainFadeSimulation import *

# Elevation and Slew Rate Settings
elvrange_case = 1  # 0 for a single elevation angle, 1 for a range
slewrange_case = 0  # 0 for a single slew rate, 1 for a range

ElevationAngle = 15
ElevationAngleMin = 15
ElevationAngleMax = 35
ElevationAngleMStep = 5

# Slew Rate Settings
slew_rate = 0.1
slewMin = 0.0
slewMax = 0.5
slewStep = 0.1

# Adjust settings based on elvrange_case
if elvrange_case == 0:
    ElevationAngleMin = ElevationAngle
    ElevationAngleMax = ElevationAngle + ElevationAngleMStep
    printing = True
else:
    ElevationAngleMax = ElevationAngleMax + ElevationAngleMStep
    printing = False

# Adjust settings based on slewrange_case
if slewrange_case == 0:
    slewMin = slew_rate
    slewMax = slew_rate + slewStep
    printing = True
else:
    slewMax = slewMax + slewStep
    printing = False

# Mapping of elevation angles to corresponding slew rates
slew_rate_mapping = {
    15: 0.1,
    20: 0.2,
    25: 0.3,
    30: 0.4,
    35: 0.5
}

# Main loop
geometry = 'downlink'  # 'uplink' or 'downlink'
for elvloop in range(ElevationAngleMin, ElevationAngleMax, ElevationAngleMStep):
    # Get the slew rate for the current elevation angle
    current_slew_rate = slew_rate_mapping.get(elvloop, slew_rate)

    # Run the main fade simulation with the current elevation and slew rate
    elevation, rms_wind_speed, Cn2, RL, SI = main_fade_simulation(
        geometry, 'day', 30000, 1.82e-05, 1.2, 1550, 1E-4, 1200, 0,
        elvloop, 5, current_slew_rate, 1.7e-14, 0.004, 1.6, transmission_losses=-2,
        altApertureDiameter=0.02, printResults=printing, C_r=0, integration_step_multiplier=0.1, nr_transmitters=1, compute_only_fades=False
    )

    print('Elevation:', elevation)
    print('wind speed:', rms_wind_speed)
    print('Slew:', current_slew_rate)
    print('Scintillation Index:', SI)
    print('Cn2:', Cn2)
