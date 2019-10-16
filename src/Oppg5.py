"""
Created on Oct 4 2019
@author: sindrhpa
"""

import numpy as np
from src.Oppg4 import *
from src.Physics import *
from src.EarthMoonData import *
from src.Atmosphere import *
from src.RK45 import *
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
print("dv:", dv)


def F(t, y):
    h, v = y
    r = h + earth_diameter_eqv/2
    Fg = -forceGravity(m(t), earth_mass, r)
    Fr = rocketThrust(t)
    Fd = -airResistance(h, rocket_Cd, rocket_cross_area, v) #-np.sign(v)
    force = Fr + Fg + Fd
    a = force/m(t)
    return np.array([v, a])


# RK45
y0 = np.array([0, 0])  # [v, h]
t0 = 0
step = 1200*1e-5
T = 1e-4
t, w = RK45(F, y0, t0, total_duration, step, T)
h, v = zip(*w)

X = list(np.linspace(0, total_duration, len(h)))

plt.plot(X, h)
plt.title('Height based on time')
plt.ylabel('height (m)')
plt.xlabel('time (s)')
plt.show()
