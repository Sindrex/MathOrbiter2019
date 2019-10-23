import src.RK45 as RK45
import src.EarthMoonData as EarthMoonData
import src.Physics as Physics
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation
matplotlib.use('TkAgg')
from timeit import default_timer

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
    time_text.set_text('')
    return line, time_text


def animate(i):
    """perform animation step"""
    global t, w, time_passed
    time_passed = default_timer() - start_anim_time
    while t < time_passed*60*60*24 or w is None:
        t, w = next(it)
        time_passed = default_timer() - start_anim_time
    x, y, vx, vy = w
    line.set_data(x, y)

    time_text.set_text(f'T = {time_passed:.2f} days')
    return line, time_text


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

    time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)

    start_anim_time = default_timer()
    time_passed = 0
    anim = animation.FuncAnimation(fig,  # figure to plot in
                                   animate,  # function that is called on each frame
                                   repeat=True,
                                   interval=0,
                                   blit=True,
                                   init_func=init  # initialization
                                   )

    plt.show()
