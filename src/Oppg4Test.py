"""
Created on Sep 30 2019
@author: sindrhpa
"""

from src.Oppg4 import *

import matplotlib.pyplot as plt
import numpy as np

print(stage1_exhaustVel, stage2_exhaustVel, stage3_exhaustVel)

# Saturn-V test
X = list(np.linspace(0, 1200, 5000))
Y = [m(i) for i in X]
Z = [rocketThrust(i) for i in X]

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time (s)')
ax1.set_ylabel('mass (kg)', color=color)
ax1.plot(X, Y, color=color)
ax1.tick_params(axis='y', labelcolor=color)
# ax1.grid(b=True)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Thrust (N)', color=color)  # we already handled the x-label with ax1
ax2.plot(X, Z, color=color)
ax2.tick_params(axis='y', labelcolor=color)
# ax2.grid(b=True)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
