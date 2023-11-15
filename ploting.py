import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt

def vb(height, Vg, Ws):
    # Bufton wind model with slew rate
    VB = Ws * height + Vg + 30 * np.exp(-((height - 9400) / 4800) ** 2)
    return VB

def compute_HV_daytime(h, Vb):
    # daytime
    # Andrews p.481 eq 1
    # Computes c2n with Hufnagel-Valley model
    term1 = 0.00594 * (Vb / 27) ** 2 * (10 ** -5 * h) ** 10 * np.exp(-h / 1000)
    term2 = 2.7 * 10 ** -16 * np.exp(-h / 1500)
    term3 = 1.7e-13  * np.exp(-h / 100)
    RI = term1 + term2 + term3
    # RI = (
    #         0.00594 * (Vb / 27) ** 2 * (10 ** -5 * h) ** 10 * np.exp(-h / 1000) +
    #         2.7 * 10 ** -16 * np.exp(-h / 1500) + 1.7e-13  * np.exp(-h / 100)
    # )
    return RI

def calculate_power_law_parameter(TH):
    if 0.75 < TH < 3.5:
        p = -0.11 * (12 - TH) ** 2 + 1.83 * (12 - TH) - 6.22
    elif 3.5 <= TH < 8.5:
        p = 1.45 - 0.02 * (TH - 6) ** 2
    elif 8.5 <= TH < 11.25:
        p = -0.048 * TH ** 2 + 0.68 * TH - 1.06
    else:
        p = np.nan  # If TH is out of bounds, return NaN
    return p

def compute_HAP_nighttime(h, Vb, sunset=18.00, sunrise=6.00, time=2.00, hg=500, h0=122):
    TP = (sunrise - sunset) / 12
    TH = (time - sunrise) / TP
    power_law_parameter = calculate_power_law_parameter(TH)
    print(TP,TH, power_law_parameter)
    term1 = 1 * (0.00594 * ((Vb / 27) ** 2) * (((h + hg) * 10 ** -5) ** 10) * np.exp(-(h + hg) / 1000))
    term2 = 2.7e-16 * np.exp(-(h + hg) / 1500)
    term3 = 1e-16 * ( h0 / h) ** power_law_parameter
    RI = term1 + term2 + term3
    return RI


MinHeight = 1
MaxHeight = 30000
stepHeight = 500
heights = np.arange(MinHeight, MaxHeight, stepHeight)   # Heights from 0 to 30000 meters
Vg = 5  # Ground speed in km/h
slewrange_case = 1
slew_rate = 0.3
slewMin = 0.0
slewMax = 0.5
slewStep = .1
geometry = 'uplink'  # or 'downlink'
model = 'night'
h0 = 1
hg = 50

# Adjust settings based on slewrange_case
if slewrange_case == 0:
    slewMin = slew_rate
    slewMax = slewMax + slewStep
    printing = True
else:
    slewMax = slewMax + slewStep
    printing = False

#plt.figure(figsize=(10, 6))
# Main loop
for Ws in np.arange(slewMin, slewMax, slewStep):
    VB_values = [vb(height, Vg, Ws) for height in heights]
#     print(f'Slew : {Ws}, RMS Wind Speed: {VB_values}')
#     plt.plot(heights, VB_values, label=f'Slew Rate = {Ws}°/s')
#
# plt.xlabel('Height (m)')
# plt.ylabel('Wind speed (m/s)')
# plt.title(f'Wind speed vs Height for Different Slew Rates [for {geometry}]')
# plt.legend()
# plt.grid(True)
# plt.show()

RI_values = []
Vbb_values = []
for h in heights:
    Vb = vb(h, Vg, slew_rate)  # Calculate VB for each height
    Vbb_values.append(Vb)
    #print(Vb, h)
    if model == 'day':
        RI = float(compute_HV_daytime(h, Vb))
        print(Vb, h, RI)
    elif model == 'night':
        RI = float(compute_HAP_nighttime(h, Vb))
    else:
        raise ValueError("Invalid model. Please choose 'day' or 'night'.")
    RI_values.append(RI)

RI_values_day = [compute_HV_daytime(h, vb(h, Vg, slew_rate)) for h in heights]
RI_values_night = [compute_HAP_nighttime(h, vb(h, Vg, slew_rate)) for h in heights]
plt.figure(figsize=(10, 6))
plt.plot(heights, Vbb_values, label=f'Slew Rate = {slew_rate}°/s')
plt.xlabel('Height (m)')
plt.ylabel('Wind speed (m/s)')
plt.title(f'Wind speed vs Height for Different Slew Rates [for {geometry}]')
plt.legend()
plt.grid(True)
plt.show()
# Plotting
plt.figure(figsize=(10, 6))
plt.semilogx(heights, RI_values_day, label='Day Model', color='red')
plt.semilogy(heights, RI_values_day)
plt.semilogx(heights, RI_values_night, label='Night Model', color='blue')
plt.semilogy(heights, RI_values_night)
plt.xlabel(r'$C_n^2 \ (m^{-2/3})$')
plt.ylabel('Height (m) - Log Scale')
plt.title('Cn2 vs Height')
plt.legend()
plt.grid(True)
plt.show()
