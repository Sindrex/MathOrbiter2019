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
    T = 1e-10
    exTime = 5000
    t, w = RK45(F, y0, t0, total_duration + exTime, step, T)
    h, v = zip(*w)


    print("H end:", h[len(h) - 1])
    print("V_end:", v[len(v) - 1])
    #plt.plot(t, h)
    #plt.title('Height based on time')
    #plt.ylabel('height (m)')
    #plt.xlabel('time (s)')
    #plt.plot(t, v)
    #plt.show()

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('height (m)', color=color)
    ax1.plot(t, h, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    # ax1.grid(b=True)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('velocity (m/s)', color=color)  # we already handled the x-label with ax1
    ax2.plot(t, v, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    # ax2.grid(b=True)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()