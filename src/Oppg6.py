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

maxr = 10
minr = np.inf
prevT = -1

Fr_angle = 0.8751520472211607 #10m 1.0889561834065051 #np.pi/2.8  # np.pi/3 works
v_angle = 0 # 0 er +x, np.pi = -x
time_up = 8*60

# DEPRECATED
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
    force_x = Fr * np.cos(Fr_angle) + Fd * v_norm_x + Fgx
    force_y = Fr * np.sin(Fr_angle) + Fd * v_norm_y + Fgy
    m_t = m(t)
    ax = force_x/m_t
    ay = force_y/m_t
    return np.array([vx, vy, ax, ay])

def F2(t, y):
    global maxr, minr, prevT, v_angle, Fr_angle
    #y0 = earth diameter/2
    x, y, vx, vy = y
    r = np.sqrt(y**2 + x**2)
    Fg = -forceGravity(m(t), earth_mass, r)
    Fgx = x/r * Fg
    Fgy = y/r * Fg
    Fr = rocketThrust(t)
    v_len = np.sqrt(vx ** 2 + vy ** 2)
    Fd = -airResistance(r-earth_diameter_eqv/2, rocket_Cd, rocket_cross_area, v_len)

    if v_len > 0:
        v_norm_x = vx / v_len
        v_norm_y = vy / v_len
    else:
        v_norm_x = 0
        v_norm_y = 1

    #print("Angle:", angle, "v_angle:", v_angle, "Fr_angle:", Fr_angle)

    if t < time_up:
        Frx = 0
        Fry = Fr
    else:
        Frx = Fr * (v_norm_x * np.cos(Fr_angle) - v_norm_y * np.sin(Fr_angle))
        Fry = Fr * (v_norm_x * np.sin(Fr_angle) + v_norm_y * np.cos(Fr_angle))

    force_x = Frx + Fd*v_norm_x + Fgx
    force_y = Fry + Fd*v_norm_y + Fgy
    m_t = m(t)
    ax = force_x/m_t
    ay = force_y/m_t

    if prevT == t:
        print("same t:", t)
    prevT = t
    #print("x:", x, "y:", y, "m:", m_t, "Fdx:", Fd*v_norm_x, "Fdy:", Fd*v_norm_y, "Fgx:", Fgx, "Fgy:", Fgy, "force_x:", force_x, "force_y:", force_y, "ax:", ax, "ay:", ay)
    return np.array([vx, vy, ax, ay])


if __name__ == "__main__":
    # RK45, med vinkel theta
    y0 = np.array([0, earth_diameter_eqv/2, 0, 0])  # [x, y, vx, vy]
    t0 = 0
    step = 1200*1e-5
    T = 1e-8
    extra = 80000
    t, w = RK45(F2, y0, t0, total_duration + extra, step, T)
    x, y, vx, vy = zip(*w)

    print("Max r:", maxr, "Min r:", minr)

    # plotting
    if True:
        plt.plot(t, y, color='y')
        plt.plot(t, [earth_diameter_eqv/2 for i in t], color='r')
        plt.plot(t, [-earth_diameter_eqv / 2 for i in t], color='r')
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

    # Animering
    fig = plt.figure()
    axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                           xlim=(-0.2e8, 0.2e8),
                           ylim=(-0.2e8, 0.2e8)
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

    Fr_angleA = 0
    Fr_angleB = 0.45*np.pi
    Fr_angleM = (Fr_angleA + Fr_angleB) / 2

    while True:
        y0 = np.array([0, earth_diameter_eqv / 2, 0, 0])  # [x, y, vx, vy]
        t0 = 0
        step = 1200 * 1e-5
        T = 1e-2
        extra = 20000

        #A
        Fr_angle = Fr_angleM
        t, w = RK45(F2, y0, t0, total_duration + extra, step, T)
        x, y, vx, vy = zip(*w)
        bonustime_i = -1
        for i, e in enumerate(t):
            if e > 1200:
                bonustime_i = i
                break

        avg_radius = sum([np.sqrt(x_i**2 + y_i**2) for x_i, y_i in zip(x[bonustime_i:], y[bonustime_i:])])/(len(t)-bonustime_i)
        metric = avg_radius - 7e6
        if metric > 0:
            Fr_angleA = Fr_angleM
        else:
            Fr_angleB = Fr_angleM
        Fr_angleM = (Fr_angleA + Fr_angleB) / 2
        print(f'M: {Fr_angleM}\tmetric: {metric}')

