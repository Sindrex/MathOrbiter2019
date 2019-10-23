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

theta = np.pi/3  # np.pi/3 works
maxr = 10
minr = np.inf

def F(t, y):
    x, y, vx, vy = y
    y_true = y + earth_diameter_eqv/2
    r2 = y_true**2 + x**2
    r = np.sqrt(r2)
    Fg = -forceGravity(m(t), earth_mass, r)
    Fgx = x/r * Fg
    Fgy = y_true/r * Fg
    print("Fgx:", Fgx, "Fgy:", Fgy)
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

def F2(t, y):
    global maxr, minr
    #y0 = earth diameter/2
    x, y, vx, vy = y
    r = (y**2 + x**2)**0.5
    #y/r = sin(tetha) => theta = sin-1(y/r)
    angle = theta
    Fg = -forceGravity(m(t), earth_mass, r)
    Fgx = x/r * Fg
    Fgy = y/r * Fg
    Fr = rocketThrust(t)
    v_len = np.sqrt(vx ** 2 + vy ** 2)
    Fd = -airResistance(r, rocket_Cd, rocket_cross_area, v_len) #-np.sign(v)
    if v_len > 0:
        v_norm_x = vx / v_len
        v_norm_y = vy / v_len
    else:
        v_norm_x = 1
        v_norm_y = 1
    if r<7.15e6: #7.2 works
        Frx = Fr*np.cos(angle)
        Fry = Fr*np.sin(angle)
    else:
        Frx = (Fr * -y/r)
        Fry = (Fr * x/r)

    force_x = Frx + Fd*v_norm_x + Fgx
    force_y = Fry + Fd*v_norm_y + Fgy
    m_t = m(t)
    ax = force_x/m_t
    ay = force_y/m_t
    if(t>10000):
        print("t:", t, "r:", r, "Frx:", Frx, "Fry:", Fry)
        if(r > maxr):
            maxr = r
        if(r < minr):
            minr = r
    #print("x:", x, "y:", y, "m:", m_t, "Fdx:", Fd*v_norm_x, "Fdy:", Fd*v_norm_y, "Fgx:", Fgx, "Fgy:", Fgy, "force_x:", force_x, "force_y:", force_y, "ax:", ax, "ay:", ay)
    return np.array([vx, vy, ax, ay])


if __name__ == "__main__":
    # RK45, med vinkel
    y0 = np.array([0, earth_diameter_eqv/2, 0, 0])  # [x, y, vx, vy]
    t0 = 0
    step = 1200*1e-5
    T = 1e-6
    extra = 20000
    t, w = RK45(F2, y0, t0, total_duration + extra, step, T)
    x, y, vx, vy = zip(*w)
    print("Max r:", maxr, "Min r:", minr)

    if(False):
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
        i *= 1
        line.set_data(x[i], y[i])
        return line,


    fig = plt.figure()
    axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                           xlim=(-0.15e8, 0.15e8),
                           ylim=(-0.15e8, 0.15e8)
                           )
    line, = axes.plot([], [], color='r', marker='o', lw=3)
    anim = animation.FuncAnimation(fig, animate, repeat=True, interval=0, blit=True, init_func=init)
    #x = np.linspace(-earth_diameter_eqv/2, earth_diameter_eqv, 500)
   # y = np.sqrt(-x**2 + earth_diameter_eqv/2)
    #axes.plot(x, y, c='blue')
    #axes.plot(x, -y, color='blue')
    axes.add_artist(plt.Circle((0,0), earth_diameter_eqv/2, color='blue'))
    axes.add_artist(plt.Circle((0, 0), earth_diameter_eqv / 2 + 250e3, color='deepskyblue', alpha=0.2))
    plt.show()