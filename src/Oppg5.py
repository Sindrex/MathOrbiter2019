"""
Created on Oct 4 2019
@author: sindrhpa
"""

import numpy as np
from src.Oppg4 import *
from src.Physics import *
from src.EarthMoonData import *
from src.Atmosphere import *
import matplotlib.pyplot as plt

# Equation from https://www.narom.no/undervisningsressurser/sarepta/rocket-theory/rocket-engines/the-rocket-equation/
# and http://www.me.umn.edu/courses/me4090_summer2017/RocketMotionSlides.pdf

g = 9.81
t = stage1_duration
ve = stage1_exhaustVel
dv1 = -ve * np.log(stage1_grossmass/stage1_drymass) - g * t

t = stage2_duration
ve = stage2_exhaustVel
dv2 = -ve * np.log(stage2_grossmass/stage2_drymass) - g * t

t = stage3_duration
ve = stage3_exhaustVel
dv3 = -ve * np.log(stage3_grossmass/stage3_drymass) - g * t

dv = dv1 + dv2 + dv3
print(dv)

# Euler
h_arr = []
dt = 1
h = 0
v = 0
for t in np.arange(0, total_duration, dt):
    r = h + earth_diameter_eqv/2
    Fg = - forceGravity(m(t), earth_mass, r)
    Fr = rocketThrust(t)
    Fd = -airResistance(h, rocket_Cd, rocket_cross_area, v) #-np.sign(v)
    force = Fr + Fg + Fd
    a = force/m(t)
    v = v + a * dt
    h = h + v * dt
    print("h:", h, "v:", v, "Fr:", Fr, "Fd:", Fd, "Fg:", Fg, "SumForce:", force)
    if(h <= 0):
        print("Crashed into earth")
        #print("Max height:", h, "m")
        break
    h_arr.append(h)

#print(h_arr)

X = list(np.linspace(0, total_duration/dt, len(h_arr)))
plt.plot(X, h_arr)
plt.show()