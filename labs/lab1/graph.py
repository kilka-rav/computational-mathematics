import numpy as np
import matplotlib.pyplot as plt

_x, _y, _y1, _y2 = 'data_x.txt', 'data1.txt', 'data2.txt', 'data3.txt'

x = np.genfromtxt(_x, dtype = 'double')
y = np.genfromtxt(_y, dtype = 'double')
y1 = np.genfromtxt(_y1, dtype = 'double')
y2 = np.genfromtxt(_y2, dtype = 'double')


# values for making ticks in x and y axis
xnumbers = np.linspace(-5, 5, 15)
ynumbers = np.linspace(-1, 1, 11)

plt.plot(x, y, color = 'b', label = 'table sin', linestyle = '-')
plt.plot(x, y1, color='r', label="table sin'", linestyle = '-') # r - red colour
plt.plot(x, y2, color='g', label="table sin''", linestyle = '-') # g - green colour
plt.xlabel("X")
plt.ylabel("Y")
plt.title("table-grid function with step = 0.25")
plt.xticks(xnumbers)
plt.yticks(ynumbers)
plt.legend()
plt.grid()
plt.axis([0, 6.5, -1.1, 1.1]) # [xstart, xend, ystart, yend]
plt.show()
