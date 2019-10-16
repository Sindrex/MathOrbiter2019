"""
Created on Sep 18 2019
@author: sindrhpa
"""

from src.Atmosphere import pressure, temperature, density, airResistance

import matplotlib.pyplot as plt
import numpy as np

# Atmosphere test
x = list(np.linspace(0, 50000, 5000))
pres = [pressure(i) for i in x]
temp = [temperature(i) for i in x]
dens = [density(i) for i in x]
airres = [airResistance(i, 1, 1, 1) for i in x]

plt.subplot(4, 1, 1)
plt.plot(x, pres)
plt.title('Pressure, Temperature, Density, Air Resistance based on height')
plt.ylabel('Pressure (Pa)')

plt.subplot(4, 1, 2)
plt.plot(x, temp)
plt.ylabel('temperature (K)')

plt.subplot(4, 1, 3)
plt.plot(x, dens)
plt.ylabel('density (kg/m^3)')

plt.subplot(4, 1, 4)
plt.plot(x, airres)
plt.xlabel('height (m)')
plt.ylabel('Air Resistance (N)')

plt.show()
