import numpy as np
import matplotlib.pyplot as plt

_x, _y, _y1 = 'dx.txt', 'error_2.txt', 'error_2_2.txt'

x = np.genfromtxt(_x, dtype = 'double')
y = np.log1p(np.genfromtxt(_y, dtype = 'double'))
y1 = np.log1p(np.genfromtxt(_y1, dtype = 'double'))


xnumbers = np.linspace(0, 0.2, 5)
ynumbers = np.linspace(0, 0.2, 5)

plt.plot(x, y, color = 'b', label = "table sin'' O(h)", linestyle = '--')
plt.plot(x, y1, color='r', label="table sin'' O(h^2)", linestyle = '--')
plt.xlabel("X")
plt.ylabel("ln(E)")
plt.title("Approximate")
plt.xticks(xnumbers)
plt.yticks(ynumbers)
plt.legend()
plt.grid()
plt.axis([0, 0.2, 0, 0.2]) # [xstart, xend, ystart, yend]
plt.show()

