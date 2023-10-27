# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import math
import csv
import copy
import sys

import numpy as np
import matplotlib.pyplot as plt
from pip._vendor.distlib.compat import raw_input
#
# sys.path.append('OpticalTurbulenceModel_v4_1/')

from MainFadeSimulation import *

# CHECK THE ELEVATION ANGLE VALUE HERE!
elvrange_case = 0               # set this to 1 or 0

# if computation needs to be done for a single elevation angle, set variable 'elvrange_case = 0'. Set the min elevation
# angle value in degrees as the required value. Leave the other parameters as it is.
ElevationAngle = 15             # in degrees

# if computation needs to be done on a range of elevation angles, set variable 'elvrange_case = 1'. Set the min and max
# elevation angle limits and the step size, all in degrees
ElevationAngleMin = 15          # min elevation angle in degrees
ElevationAngleMax = 35          # max elevation angle in degrees
ElevationAngleMStep = 5         # step size in degrees
MinHeight = 0
MaxHeight = 30000
if elvrange_case == 0:
    ElevationAngleMin = ElevationAngle
    ElevationAngleMax = ElevationAngleMin + ElevationAngleMStep
    printing = True
elif elvrange_case == 1:
    ElevationAngleMax = ElevationAngleMax + ElevationAngleMStep
    printing = False
for heightloop in range(MinHeight, MaxHeight, 100):
    for elvloop in range(ElevationAngleMin, ElevationAngleMax, ElevationAngleMStep):
         main_fade_simulation('uplink', 'day', heightloop,1.82e-05, 1.2, 1550,
                                                                           1E-4,  1200, 0,
                                                                           elvloop, 8, 0.1,
                                                                           hv_ground_cst=1.7e-14,
                                                                           altApertureDiameter=0.02,
                                                                           printResults=printing, C_r=0,
                                                                           integration_step_multiplier=0.1, nr_transmitters=1,
                                                                           compute_only_fades=False)

        #print(elvloop, mean_fade_loss, fadeLevel_dB, surgeLevel_dB)