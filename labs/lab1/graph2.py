import numpy as np
import matplotlib.pyplot as plt

 # Plotting tutorials in Python
# Adding Multiple plots by superimposition
# Good for plots sharing similar x, y limits
# Using multiple plot commands
# Much better and preferred than previous

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-5, 5, 1000)
y = np.sin(x)
z = np.cos(x)
z1 = -1 * np.sin(x)

# values for making ticks in x and y axis
xnumbers = np.linspace(-5, 5, 15)
ynumbers = np.linspace(-1, 1, 11)

plt.plot(x, y, color='r', label='sin') # r - red colour
plt.plot(x, z, color='g', label='cos') # g - green colour
plt.plot(x, z1, color='b', label='-sin')
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Original functions")
plt.xticks(xnumbers)
plt.yticks(ynumbers)
plt.legend()
plt.grid()
plt.axis([0, 6.5, -1.1, 1.1]) # [xstart, xend, ystart, yend]
plt.show()
