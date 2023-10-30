MinHeight = 0
MaxHeight = 30000

# Uplink: Height increases from MinHeight to MaxHeight with a step of 500
for height in range(MinHeight, MaxHeight + 1, 500):
    print("Uplink:", height)

# Downlink: Height decreases from MaxHeight to MinHeight with a step of -500
for height in range(MaxHeight, MinHeight - 1, -500):
    print("Downlink:", height)
