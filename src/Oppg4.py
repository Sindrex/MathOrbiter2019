"""
Created on Sep 30 2019
@author: sindrhpa
"""

import numpy as np

# SATURN-V

# N
stage1_thrust = 35100*1000  # N
stage2_thrust = 5141*1000  # N
stage3_thrust = 1033.1*1000  # N

# kg
stage1_drymass = 130000
stage1_grossmass = 2290000
stage1_fuelmass = stage1_grossmass - stage1_drymass
stage2_drymass = 40100
stage2_grossmass = 496200
stage2_fuelmass = stage2_grossmass - stage2_drymass
stage3_drymass = 13500
stage3_grossmass = 123000
stage3_fuelmass = stage3_grossmass - stage3_drymass
total_grossmass = stage1_grossmass + stage2_grossmass + stage3_grossmass
final_grossmass = stage3_grossmass - stage3_fuelmass

# burn duration, s
stage1_duration = 168
stage2_duration = 360
stage3_duration = 165 + 335
total_duration = stage1_duration + stage2_duration + stage3_duration

# delta mass, kg/s
stage1_dmass = stage1_fuelmass/stage1_duration
stage2_dmass = stage2_fuelmass/stage2_duration
stage3_dmass = stage3_fuelmass/stage3_duration

# exhaust velocity, ve = -F/dm
stage1_exhaustVel = -stage1_thrust/stage1_dmass
stage2_exhaustVel = -stage2_thrust/stage2_dmass
stage3_exhaustVel = -stage3_thrust/stage3_dmass

# other values
rocket_diameter = 10.1  # m
rocket_cross_area = np.pi * (rocket_diameter/2)**2  # m^2
rocket_Cd = 0.515

# mass after t sec
def m(t):
    delta = 0
    if t <= stage1_duration:
        delta = stage1_dmass * t
    elif t <= (stage1_duration + stage2_duration):
        t1 = t - stage1_duration
        delta = stage1_grossmass + stage2_dmass * t1
    elif t <= (stage1_duration + stage2_duration + stage3_duration):
        t1 = t - stage1_duration - stage2_duration
        delta = stage1_grossmass + stage2_grossmass + stage3_dmass * t1
    else:
        delta = stage1_grossmass + stage2_grossmass + stage3_fuelmass
    return total_grossmass - delta

# thrust after t sec
def rocketThrust(t):
    if t < stage1_duration:
        return stage1_thrust
    elif t < (stage1_duration + stage2_duration):
        return stage2_thrust
    elif t < (stage1_duration + stage2_duration + stage3_duration):
        return stage3_thrust
    else:
        return 0
