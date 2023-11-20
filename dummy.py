# Create an instance of ComputeWindSpeed
compute_wind_speed = ComputeWindSpeed(
    Vg=10,                    # Ground Velocity in m/s
    wind_height_min=0,        # Minimum Wind Height in km
    wind_height_max=15,       # Maximum Wind Height in km
    elevation=15,             # Elevation in degrees
    SLT_altitude=300,         # Satellite LEO Altitude in km
    geometry='uplink'         # Geometry ('uplink' or 'downlink')
)

# Calculate slew rate
slew_rate = compute_wind_speed.calculate_slew()

# Print the calculated slew rate
print(f"Slew Rate: {slew_rate} degrees/s")

# Calculate RMS wind speed
rms_wind_speed = compute_wind_speed.rms_wind_speed

# Print the calculated RMS wind speed
print(f"RMS Wind Speed: {rms_wind_speed} m/s")
