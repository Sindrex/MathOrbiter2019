import numpy as np
from timeit import default_timer
from matplotlib import pyplot as plt
from src import Oppg6


def RK45(F, y_0, t_0, t_1, h_0, T):
    """Runge Kutta Fehlberg, a Runge Kutta 4/5 embedded pair method for iterative solving of IVPs

    Parameters:
    F (function(t, y)): the derivative of y(t)
    y_0 (number(s)): the initial values for y
    t_0: start time
    t_1: end time
    h: initial step size
    T: error tolerance
    """
    h = h_0
    t = [t_0]
    w = [y_0]
    while t[-1] + h < t_1:
        h_new, t_new, w_new = RK45step(F, t[-1], w[-1], h, T)
        t.append(t_new)
        w.append(w_new)
        h = h_new

    while t[-1] < t_1:
        h = t_1 - t[-1]
        _, t_new, w_new = RK45step(F, t[-1], w[-1], h, T)  # get end value
        t.append(t_new)
        w.append(w_new)
    return t, w


def RK45Iterator(F, y_0, t_0, t_1, h_0, T):
    """Runge Kutta Fehlberg, a Runge Kutta 4/5 embedded pair method for iterative solving of IVPs

    Parameters:
    F (function(t, y)): the derivative of y(t)
    y_0 (number(s)): the initial values for y
    t_0: start time
    t_1: end time
    h: initial step size
    T: error tolerance
    """
    h = h_0
    t = t_0
    w = y_0
    while t + h < t_1:
        h, t, w = RK45step(F, t, w, h, T)
        yield t, w

    while t < t_1:
        h = t_1 - t
        _, t, w = RK45step(F, t, w, h, T)  # get end value
        yield t, w
    return t, w


def RK45step(F, t, w, h, T, retry=False):
    """A single step of RK45

    Parameters:
    F (function(t, y)): the derivative of y(t)
    t: number(s), the current x value
    w: the current approximation of y
    h: the step size
    T: error tolerance
    """

    while True:  # loops until return
        s1 = F(t, w)
        s2 = F(t + 1/4 * h, w + 1/4 * h*s1)
        s3 = F(t + 3/8 * h, w + 3/32 * h*s1 + 9/32 * h*s2)
        s4 = F(t + 12/13 * h, w + 1932/2197 * h*s1 - 7200/2197 * h*s2 + 7296/2197 * h*s3)
        s5 = F(t + h, w + 439/216 * h*s1 - 8*h*s2 + 3680/513 * h*s3 - 845/4104*h*s4)
        s6 = F(t + 1/2 * h, w - 8/27 * h*s1 + 2 * h*s2 - 3544/2565 * h*s3 + 1859/4104 * h*s4 - 11/40 * h*s5)
        w_new = w + h*(25/216 * s1 + 1408/2565 * s3 + 2197/4104 * s4 - 1/5 * s5)
        z_new = w + h*(16/135 * s1 + 6656/12825 * s3 + 28561/56430 * s4 - 9/50 * s5 + 2/55 * s6)
        t_new = t + h
        e = abs(np.linalg.norm(z_new - w_new))
        e = max(e, T**2) # protect against zero error making infinite jumps
        if np.isinf(e):
            e = T*1000
        w_norm = np.linalg.norm(w_new)
        if abs(w_norm) < 1e-6:
            w_norm = 1e-5

        if retry:
            h_new = h/2
        else:
            h_new = 0.8 * h * np.power(T / e, 0.2)

        if h_new == 0:
            h_new = T**5
        if e/w_norm > T:
            h = h_new
            retry = True
        else:
            return h_new, t_new, w_new


def _test_y1(t):
    return np.e**t


def _test_F1(t, y):
    """y(t) = e^t"""
    return y


def _test_y2(t):
    return np.e ** (t**3 / 3)


def _test_F2(t, y):
    return t*t*y


def _test_y3(t):
    return np.array([t, 2*t*t, np.e**(t*t)])


def _test_F3(t, y):
    return np.array([1, 4*t, 2*t*np.e**(t*t)])


def _test_errors():
    T = 1e-12
    h_0 = T/2
    t_0 = -1
    t_1 = 3

    y_0_1 = _test_y1(t_0)
    y_0_2 = _test_y2(t_0)
    y_0_3 = _test_y3(t_0)

    tees_1, w_1 = RK45(_test_F1, y_0_1, t_0, t_1, h_0, T)
    iterator2 = RK45Iterator(_test_F2, y_0_2, t_0, t_1, h_0, T)
    tees_3, w_3 = RK45(_test_F3, y_0_3, t_0, t_1, h_0, T)
    y_1 = [_test_y1(t) for t in tees_1]
    y_3 = [_test_y3(t) for t in tees_3]
    g_1 = [(y - w) for y, w in zip(y_1, w_1)]
    g_3 = [np.linalg.norm((y - w)) for y, w in zip(y_3, w_3)]
    r_1 = [(y - w)/w for y, w in zip(y_1, w_1)]
    r_3 = [np.linalg.norm((y - w))/np.linalg.norm(w) for y, w in zip(y_3, w_3)]

    print("RK45")
    print("t:            | w:            | y:            | error:        | rel. error")
    for i in range(len(tees_1)):
        print(f'{tees_1[i]:+5e} | {w_1[i]:+5e} | {y_1[i]:+5e} | {g_1[i]:+5e} | {r_1[i]:+5e}')

    print("RK45Iterator")
    print("t:            | w:            | y:            | error:        | rel. error")
    for t, w in iterator2:
        print(f'{t:+5e} | {w:+5e} | {_test_y2(t):+5e} | {(_test_y2(t)-w):+5e} | {(_test_y2(t)-w) / w:+5e}')

    print("RK45 with multivariable")
    print("t:            | error:         | rel. error")
    for i in range(len(tees_3)):
        print(f'{tees_3[i]:+5e} | {g_3[i]:+5e} | {r_3[i]:+5e}')


def _test_runtime():
    t_0 = 0
    t_1 = 1200

    y_0 = [0, 0, 0, 0]

    times = []
    tolerances = np.logspace(-1, -16, num=20)
    print("starting time tests:")
    for T in tolerances:
        start = default_timer()
        t, w = RK45(Oppg6.F, y_0, t_0, t_1, T, T)
        end = default_timer()
        print("solution w tolerance",T,"has length",len(t))
        times.append(end - start)
        print(end-start)

    ax = plt.gca()
    ax.plot(tolerances, times)
    print(*zip(tolerances, times), sep=', ')
    ax.set_xlabel("Tolerance")
    ax.set_ylabel("Time (seconds)")

    plt.show()


def _test_realtime():
    t_0 = 0
    t_1 = 30

    y_0 = [0, 0, 0, 0]

    times = []
    T = 1e-25
    it = RK45Iterator(Oppg6.F, y_0, t_0, t_1, T, T)
    print("starting time tests:")
    start = default_timer()
    for t, w in it:
        pass
    end = default_timer()
    print(end-start)


if __name__ == '__main__':
    #_test_errors()
    #_test_runtime()
    _test_realtime()

    t_0 = 0
    t_1 = 1200

    y_0 = [0, 0, 0, 0]
