import numpy as np
import matplotlib.pyplot as plt
from MainFadeSimulation import *

# Elevation and Slew Rate Settings
elvrange_case = 1  # 0 for a single elevation angle, 1 for a range
slewrange_case = 0 # 0 for a single slew rate, 1 for a range

ElevationAngle = 15 # change elevation angle in case only one elevation angle is to be considered
ElevationAngleMin = 15
ElevationAngleMax = 35
ElevationAngleStep = 5

# Slew Rate Settings
Slew_rate = 0.2 #change slew rate based on user requirement in case  only one slew rate is to be considered
Slew_rateMin = 0.1
Slew_rateMax = 0.5
Slew_rateStep = 0.2

# Adjust settings based on elvrange_case
if elvrange_case == 0:
    ElevationAngleMin = ElevationAngle
    ElevationAngleMax = ElevationAngle + ElevationAngleStep
    printing = True
else:
    ElevationAngleMax = ElevationAngleMax + ElevationAngleStep
    printing = False

# Adjust settings based on slewrange_case
if slewrange_case == 0:
    Slew_rateMin = Slew_rate
    Slew_rateMax = Slew_rate + Slew_rateStep
    printing = True
else:
    Slew_rateMax = Slew_rateMax + Slew_rateStep
    printing = False


# Main loop
geometry = 'uplink'  # 'uplink' or 'downlink'
for elvloop in range(ElevationAngleMin, ElevationAngleMax, ElevationAngleStep):
    # Get the slew rate for the current elevation angle
    for slewloop in np.arange(Slew_rateMin, Slew_rateMax, Slew_rateStep):

        # Run the main fade simulation with the current elevation and slew rate
        mean_fade_loss, fadeLevel_dB, surgeLevel_dB, elevation, rms_wind_speed, Cn2, RL, SI = main_fade_simulation(
            geometry, 'day', 1, 30000, 18.00,6.00,4.00,
            1.82e-05, 1.2, 1550, 1E-4, 1200, 0, elvloop, 5,
            slewloop, 1.7e-14, 0.004, 1.6, transmission_losses=-2, altApertureDiameter=0.02,
            printResults=printing, C_r=0, integration_step_multiplier=0.1, nr_transmitters=1, compute_only_fades=False
        )

        print('Elevation:', elevation)
        print('wind speed:', rms_wind_speed)
        print('Slew:', Slew_rate)
        print('Scintillation Index:', SI)
        print('Cn2:', Cn2)
        print('Mean fade loss:', mean_fade_loss)
        print('Fade Level in dB:', fadeLevel_dB)
        print('Surge Level in dB:', surgeLevel_dB)
