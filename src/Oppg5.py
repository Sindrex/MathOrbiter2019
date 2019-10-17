"""
Created on Oct 4 2019
@author: sindrhpa
"""

from src.Oppg4 import *
from src.Physics import *
from src.EarthMoonData import *
from src.Atmosphere import *
from src.RK45 import *
import matplotlib.pyplot as plt

def F(t, y):
    h, v = y
    r = h + earth_diameter_eqv/2
    Fg = -forceGravity(m(t), earth_mass, r)
    Fr = rocketThrust(t)
    Fd = -airResistance(h, rocket_Cd, rocket_cross_area, v) #-np.sign(v)
    force = Fr + Fg + Fd
    a = force/m(t)
    return np.array([v, a])

if __name__ == "__main__":
    # RK45
    y0 = np.array([0, 0])  # [v, h]
    t0 = 0
    step = 1200*1e-5
    T = 1e-4
    t, w = RK45(F, y0, t0, total_duration, step, T)
    h, v = zip(*w)

    X = list(np.linspace(0, total_duration, len(h)))

    plt.plot(t, h)
    plt.title('Height based on time')
    plt.ylabel('height (m)')
    plt.xlabel('time (s)')
    plt.show()
