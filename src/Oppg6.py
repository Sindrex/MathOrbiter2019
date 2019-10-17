"""
Created on Oct 17 2019
@author: sindrhpa
"""

from src.Oppg4 import *
from src.Physics import *
from src.EarthMoonData import *
from src.Atmosphere import *
from src.RK45 import *
import matplotlib.pyplot as plt

theta = 1/2 * np.pi

def F(t, y):
    x, y, vx, vy = y
    r = y + earth_diameter_eqv/2
    Fg = -forceGravity(m(t), earth_mass, r)
    Fr = rocketThrust(t)
    v_len = np.sqrt(vx ** 2 + vy ** 2)
    Fd = -airResistance(y, rocket_Cd, rocket_cross_area, v_len) #-np.sign(v)
    if v_len > 0:
        v_norm_x = vx / v_len
        v_norm_y = vy / v_len
    else:
        v_norm_x = 0
        v_norm_y = 0
    force_x = Fr*np.cos(theta) + Fd*v_norm_x
    force_y = Fr*np.sin(theta) + Fd*v_norm_y + Fg
    m_t = m(t)
    ax = force_x/m_t
    ay = force_y/m_t
    return np.array([vx, vy, ax, ay])


if __name__ == "__main__":
    # RK45
    y0 = np.array([0, 0, 0, 0])  # [x, y, vx, vy]
    t0 = 0
    step = 1200*1e-5
    T = 1e-4
    t, w = RK45(F, y0, t0, total_duration, step, T)
    x, y, vx, vy = zip(*w)
    print(x)
    print(y)
    print(vx)
    print(vy)
    X = list(np.linspace(0, total_duration, len(y)))

    plt.plot(t, y)
    plt.title('Height based on time')
    plt.ylabel('height (m)')
    plt.xlabel('time (s)')
    plt.show()
