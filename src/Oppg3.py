import src.RK45 as RK45
import src.EarthMoonData as EarthMoonData
import src.Physics as Physics
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation
matplotlib.use('TkAgg')

def F(t, y):
    x, y, vx, vy = y
    dist = (x**2 + y**2)**0.5
    f = Physics.forceGravity(EarthMoonData.moon_mass, EarthMoonData.earth_mass, dist)
    a = f / EarthMoonData.moon_mass
    ax = - a * x / dist
    ay = - a * y / dist
    return np.array([
        vx, vy, ax, ay
    ])


def init():
    """initialize animation"""
    line.set_data([], [])
    return line,


def animate(i):
    """perform animation step"""
    global t, w
    while t < i*60*60*24 * 16.666/1000 or w is None:
        t, w = next(it)
    x, y, vx, vy = w
    line.set_data(x, y)
    return line,


if __name__ == '__main__':
    y0 = np.array([
        (EarthMoonData.moon_apoapsis + EarthMoonData.moon_periapsis) / 2,
        0,
        0,
        -EarthMoonData.moon_velocity])  # [x_moon, y_moon, vx_moon, vy_moon]
    t0 = t = 0
    w = None
    step = 1e-8
    T = 1e-6
    total_duration = 60*60*24*30 * 10
    it = RK45.RK45Iterator(F, y0, t0, total_duration, step, T)
    fig = plt.figure()
    axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                           xlim=(-EarthMoonData.moon_apoapsis*1.1, EarthMoonData.moon_apoapsis*1.1),
                           ylim=(-EarthMoonData.moon_apoapsis*1.1, EarthMoonData.moon_apoapsis*1.1))

    line, = axes.plot([], [], color='grey', marker='o', lw=2)
    axes.scatter([0], [0], c='blue')
    anim = animation.FuncAnimation(fig,  # figure to plot in
                                   animate,  # function that is called on each frame
                                   repeat=True,
                                   interval=16.66666,
                                   blit=True,
                                   init_func=init  # initialization
                                   )

    plt.show()
