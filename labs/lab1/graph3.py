import numpy as np
import matplotlib.pyplot as plt

_x, _y, _y1 = 'dx.txt', 'error_1.txt', 'error_2.txt'

x = np.genfromtxt(_x, dtype = 'double')
y = np.log1p(np.genfromtxt(_y, dtype = 'double'))
y1 = np.log1p(np.genfromtxt(_y1, dtype = 'double'))


xnumbers = np.linspace(0, 1, 5)
ynumbers = np.linspace(0, 1, 10)

plt.plot(x, y, color = 'b', label = "table sin'", linestyle = '--')
plt.plot(x, y1, color='r', label="table sin''", linestyle = '--')
plt.xlabel("X")
plt.ylabel("ln(E+1)")
plt.title("Approximate")
plt.xticks(xnumbers)
plt.yticks(ynumbers)
plt.legend()
plt.grid()
plt.axis([0, 6.5, -1.1, 1.1]) # [xstart, xend, ystart, yend]
plt.show()


