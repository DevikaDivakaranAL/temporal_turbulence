# Defining the functions for calculating RI for day and night models

def vb(height, Vg, Ws):
    # Bufton wind model with slew rate
    VB = Ws * height + Vg + 30 * np.exp(-((height - 9400) / 4800) ** 2)
    return VB

def compute_HV_daytime(h, Vb):
    # Computes c2n with Hufnagel-Valley model for daytime
    term1 = 0.00594 * (Vb / 27) ** 2 * (10 ** -5 * h) ** 10 * np.exp(-h / 1000)
    term2 = 2.7 * 10 ** -16 * np.exp(-h / 1500)
    term3 = 1.7e-13  * np.exp(-h / 100)
    RI = term1 + term2 + term3
    return RI

def compute_HAP_nighttime(h, Vb, power_law_parameter=0.67, hg=500, h0=122):
    # Computes c2n with HAP model for nighttime
    term1 = 1 * (0.00594 * ((Vb / 27) ** 2) * (((h + hg) * 10 ** -5) ** 10) * np.exp(-(h + hg) / 1000))
    term2 = 2.7e-16 * np.exp(-(h + hg) / 1500)
    term3 = 1.7e-14 * ( h0 / h) ** power_law_parameter
    RI = term1 + term2 + term3
    return RI

# Height parameters
MinHeight = 1
MaxHeight = 30000
stepHeight = 500
heights = np.arange(MinHeight, MaxHeight, stepHeight)  # Heights from 1 to 30000 meters

# Wind model parameters
Vg = 5  # Ground speed in m/s
slew_rate = 0.3  # Slew rate in degrees/s

# Calculating RI values for daytime and nighttime models
RI_values_day = [compute_HV_daytime(h, vb(h, Vg, slew_rate)) for h in heights]
RI_values_night = [compute_HAP_nighttime(h, vb(h, Vg, slew_rate)) for h in heights]

# Plotting
plt.figure(figsize=(10, 6))
plt.semilogx(RI_values_day, heights, label='Day Model')
plt.semilogy(RI_values_day, heights)
plt.semilogx(RI_values_night, heights, label='Night Model')
plt.semilogy(RI_values_night, heights)
plt.xlabel(r'$C_n^2 \ (m^{-2/3})$ - Log Scale')
plt.ylabel('Height (m) - Log Scale')
plt.title('Cn2 vs Height for Day and Night Models')
plt.legend()
plt.grid(True)
plt.show()
