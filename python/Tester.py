from python.Atmosphere import pressure, temperature, density

height = 40000  # meters
pressure = pressure(height)
temp = temperature(height)
dens = density(height)
print("Pressure:", pressure, "- Temperature:", temp, "- Density:", dens)
