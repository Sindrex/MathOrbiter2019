"""
Created on Sep 18 2019
@author: sindrhpa
"""

import numpy as np

_troposLim = 11000  # m
_lstratLim = 25000  # m

def airResistance(h, Cd, A, v):
    return 1/2 * Cd * density(h) * A * v*v  # N

def density(h):
    k = 3.4855  # K*s^2/m^2
    return pressure(h)/temperature(h)*k  # kg/m^3

def pressure(h):  # Pa
    if h <= _troposLim:
        return _troposP(h)
    elif h <= _lstratLim:
        return _lstratP(h)
    else:
        return _ustratP(h)

def temperature(h):  # K
    if h <= _troposLim:
        return _troposT(h)
    elif h <= _lstratLim:
        return _lstratT(h)
    else:
        return _ustratT(h)

def _troposP(h):
    k1 = 101.29*1000  # Pa
    k2 = 288.08  # Kelvin
    k3 = 5.256
    return k1 * (_troposT(h)/k2)**k3

def _troposT(h):
    k1 = 288.19  # Kelvin
    k2 = 0.00649  # K/m
    return k1 - k2*h

def _lstratP(h):
    k1 = 127.76*1000  # Pa
    k2 = 0.000157  # m^-1
    return k1 * np.exp(-(k2)*h)

def _lstratT(h):
    return 216.69  # Kelvin

def _ustratP(h):
    k1 = 2.488*1000  # Pa
    k2 = 216.6  # Kelvin
    k3 = -11.388
    return k1*(_ustratT(h)/k2)**k3

def _ustratT(h):
    k1 = 141.94  # Kelvin
    k2 = 0.00299  # K/m
    return k1 + k2*h
