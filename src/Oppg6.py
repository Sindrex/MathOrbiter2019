"""
Created on Oct 17 2019
@author: sindrhpa
"""
from matplotlib import animation

from src.Oppg4 import *
from src.Physics import *
from src.EarthMoonData import *
from src.Atmosphere import *
from src.RK45 import *
import matplotlib.pyplot as plt

theta = np.pi/4.48  # 1/2 * np.pi (90 grader)

def F(t, y):
    x, y, vx, vy = y
    y_true = y + earth_diameter_eqv/2
    r2 = y_true**2 + x**2
    r = np.sqrt(r2)
    Fg = -forceGravity(m(t), earth_mass, r)
    Fgx = x/r * Fg
    Fgy = y_true/r * Fg
    Fr = rocketThrust(t)
    v_len = np.sqrt(vx ** 2 + vy ** 2)
    Fd = -airResistance(y, rocket_Cd, rocket_cross_area, v_len) #-np.sign(v)
    if v_len > 0:
        v_norm_x = vx / v_len
        v_norm_y = vy / v_len
    else:
        v_norm_x = 0
        v_norm_y = 0
    force_x = Fr*np.cos(theta) + Fd*v_norm_x + Fgx
    force_y = Fr*np.sin(theta) + Fd*v_norm_y + Fgy
    m_t = m(t)
    ax = force_x/m_t
    ay = force_y/m_t
    return np.array([vx, vy, ax, ay])


if __name__ == "__main__":
    # RK45, med vinkel
    y0 = np.array([0, 0, 0, 0])  # [x, y, vx, vy]
    t0 = 0
    step = 1200*1e-5
    T = 1e-6
    extra = 10000
    t, w = RK45(F, y0, t0, total_duration + extra, step, T)
    x, y, vx, vy = zip(*w)
    #print(x)
    #print(y)
    print("Highest x:", max(x))
    print("Highext y:", max(y))
    #print(vx)
    print("Highext vx:", max(vx))
    print("Highext vy:", max(vy))
    #print(vy)

    plt.plot(t, y, color='y')
    plt.plot(t, [earth_diameter_eqv/2 for i in t], color='r')
    plt.title('Height based on time')
    plt.ylabel('height (m)')
    plt.xlabel('time (s)')
    plt.plot(t, x, color='m')
    plt.show()

    def init():
        line.set_data([0], [0])
        return line,

    def animate(i):
        line.set_data(x[i], y[i])
        return line,


    fig = plt.figure()
    axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                           xlim=(-40.e6, 40.e6),
                           ylim=(-40.e6, 40.e6)
                           )
    line, = axes.plot([], [], color='r', marker='o', lw=3)
    anim = animation.FuncAnimation(fig, animate, repeat=True, interval=16, blit=True, init_func=init)
    plt.scatter(0, 0, s=50, c='blue')
    plt.show()