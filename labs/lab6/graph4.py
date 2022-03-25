import numpy as np
import matplotlib.pyplot as plt

_x, _y = "x_runge_kutta.txt", "y_runge_kutta.txt"

x = np.genfromtxt(_x, dtype = 'double')
y = np.genfromtxt(_y, dtype = 'double')

y_orig = np.exp(x ** 2)

plt.scatter(x, y, color = 'r', label = 'solution', linestyle = ':')
plt.plot(x, y_orig, color = 'black', linestyle = ':')
plt.xlabel("X")
plt.ylabel("Y")
plt.title("runge_kutta")
plt.grid()
plt.show()
