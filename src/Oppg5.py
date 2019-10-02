import numpy as np

from src.Oppg4 import *

# Equation from https://www.narom.no/undervisningsressurser/sarepta/rocket-theory/rocket-engines/the-rocket-equation/
# and http://www.me.umn.edu/courses/me4090_summer2017/RocketMotionSlides.pdf

g = 9.81
t = stage1_duration
ve = stage1_exhaustVel
dv1 = -ve * np.log(stage1_grossmass/stage1_drymass) - g * t

t = stage2_duration
ve = stage2_exhaustVel
dv2 = -ve * np.log(stage2_grossmass/stage2_drymass) - g * t

t = stage3_duration
ve = stage3_exhaustVel
dv3 = -ve * np.log(stage3_grossmass/stage3_drymass) - g * t

dv = dv1 + dv2 + dv3
print(dv)

print(dv*t)