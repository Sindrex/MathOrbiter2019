from src.Atmosphere import pressure, temperature, density
from src.Oppg4 import m, rocketThrust

# Atmosphere test
height = 40000  # meters
pressure = pressure(height)
temp = temperature(height)
dens = density(height)
print("Pressure:", pressure, "- Temperature:", temp, "- Density:", dens)

# Saturn-V test
for t in range(0, 1200, 100):
    print("t:", t, "Mass:", m(t), "Thrust:", rocketThrust(t))