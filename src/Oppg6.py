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
import matplotlib
matplotlib.use('TkAgg')

prevT = -1

Fr_angle = 1.2658864862178714 #better than 10min: 1.147902139420116
# gives apoapsis 1.017294e+07, periapsis 9.000022e+06
v_angle = 0 # 0 er +x, np.pi = -x
time_up = 12*60

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
    T = 1e-10
    step = T/10
    extra = 3*60*60
    t, w = RK45(F2, y0, t0, total_duration + extra, step, T)
    x, y, vx, vy = zip(*w)

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
        dot.set_data([0], [0])
        return line, dot

    def animate(i):
        t_to_meet = 5*60*i/60
        correct_i = -1
        for i, e in enumerate(t):
            if e > t_to_meet:
                correct_i = i
                break

        print(t)
        line.set_data(x[:correct_i], y[:correct_i])
        dot.set_data(x[correct_i], y[correct_i])
        return line, dot

    # Animering
    fig = plt.figure()
    axes = fig.add_subplot(111, aspect='equal', autoscale_on=True,
                           xlim=(-1.5e7, 1.5e7),
                           ylim=(-1.5e7, 1.5e7)
                           )
    axes.grid(True, 'major')
    line, = axes.plot([], [], color='r', lw=1)
    dot, = axes.plot([], [], color='r', marker='o', lw=3)
    anim = animation.FuncAnimation(fig, animate, repeat=True, interval=(1/60), blit=True, init_func=init)
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

    prev_metric = 0
    T = 1e-12
    while True:
        y0 = np.array([0, earth_diameter_eqv / 2, 0, 0])  # [x, y, vx, vy]
        t0 = 0

        step = 1

        #A
        Fr_angle = Fr_angleM
        t, w = RK45(F2, y0, t0, total_duration + extra, step, T)
        x, y, vx, vy = zip(*w)
        extra_time_i = -1
        for i, e in enumerate(t):
            if e > total_duration:
                extra_time_i = i
                break

        radii = [np.sqrt(x_i ** 2 + y_i ** 2) for x_i, y_i in zip(x[extra_time_i:], y[extra_time_i:])]
        periapsis = min(radii)
        apoapsis = max(radii)
        metric = periapsis - 9e6

        if metric > 0:
            Fr_angleA = Fr_angleM
        else:
            Fr_angleB = Fr_angleM
        Fr_angleM = (Fr_angleA + Fr_angleB) / 2
        print(f'M: {Fr_angleM}\tmetric: {metric}\tapo: {apoapsis:5e}\tperi: {periapsis:5e}')
        if metric - prev_metric == 0:
            break
        prev_metric = metric

