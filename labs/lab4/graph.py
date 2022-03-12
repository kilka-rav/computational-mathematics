import numpy as np
import matplotlib.pyplot as plt

_x, _y = "x_newton.txt", "y_newton.txt"

_x_orig, _y_orig = "x_cub.txt", "y_cub.txt"

x_orig = np.genfromtxt(_x_orig, dtype = 'double')
y_orig = np.genfromtxt(_y_orig, dtype = 'double')




x = np.genfromtxt(_x, dtype = 'double')
y = np.genfromtxt(_y, dtype = 'double')

y1 = np.log(x) + np.sqrt(1 + x)

plt.scatter(x, y, color = 'r', label = 'newton', linestyle = ':')
plt.plot(x, y1, color = 'b', label = 'newton', linestyle = '-')
plt.scatter(x_orig, y_orig, color = 'black', linestyle = ':')
plt.xlabel("X")
plt.ylabel("Y")
plt.title("newton - r, original - b")
plt.grid()
plt.show()
