"""
Created on Sep 30 2019
@author: sindrhpa
"""

def forceN2(m, a):
    return m*a

def forceGravity(m1, m2, r):
    g = 6.67430*10**-11  # m^3*kg^-1*s^-2
    return g*(m1*m2)/(r*r)

def forceRocketEngine(dm, ve):
    # F = -dm*ve
    # dm = derivert masse av tid
    # ve = relativ fart tid eksos i forhold til rakett
    return -dm*ve
