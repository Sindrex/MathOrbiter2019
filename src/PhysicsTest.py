"""
Created on Oct 21 2019
@author: sindrhpa
"""

import matplotlib.pyplot as plt
from src.Physics import *
from src.EarthMoonData import *
from src.Oppg4 import *

m = m(1300)

hmax = 8 * 10**6

X = list(np.linspace(0, hmax, 5000))
Y = [forceGravity(m, earth_mass, i + earth_diameter_eqv/2) for i in X]

plt.plot(X, Y)
plt.show()