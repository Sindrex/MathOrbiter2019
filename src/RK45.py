import numpy as np
import sys


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
    w_norm = np.linalg.norm(w_new)

    if retry:
        h_new = h/2
    else:
        h_new = 0.8 * h * np.power(T / e, 0.2)

    if e/w_norm > T:
        return RK45step(F, t, w, h_new, T, retry=True)
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


if __name__ == '__main__':
    T = 1e-6
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

    print("global feil 1:", g_1[-1], '\nakk. lokal feil 1:', sum(g_1))
    #print("t:      ", '\t|\t', "w:    ", '\t|\t', "y:      ", '\t|\t', "error:")
    #for i in range(len(tees_1)):
    #    print(f'{tees_1[i]:+5f}\t|\t{w_1[i]:+5f}\t|\t{y_1[i]:+5f}\t|\t{e_1[i]:+5f}')

    print("t:      ", '\t|\t', "w:    ", '\t|\t', "y:      ", '\t|\t', "error:")
    for t, w in iterator2:
        print(f'{t:+5f}\t|\t{w:+5f}\t|\t{_test_y2(t):+5f}\t|\t{(_test_y2(t)-w):+5f}')

    #print("t:      ", '\t|\t', "w:    ", '\t|\t', "y:      ", '\t|\t', "error:")
    #for t, w in iterator2):
    #    print(f'{tees_3[i]:+5f}\t|\t{e_3[i]:+5f}')
