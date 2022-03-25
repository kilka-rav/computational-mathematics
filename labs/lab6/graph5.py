import numpy as np
import matplotlib.pyplot as plt

_x, _y , _y1, _y2, _y3 = "step.txt", "ex_err.txt", "im_err.txt", "recalc_err.txt", "runge_err.txt"

x = np.genfromtxt(_x, dtype = 'double')
y = np.genfromtxt(_y, dtype = 'double')
y1 = np.genfromtxt(_y1, dtype = 'double')
y2 = np.genfromtxt(_y2, dtype = 'double')
y3 = np.genfromtxt(_y3, dtype = 'double')

#x = 1 / x

x = np.log(x)
y = np.log(y)
y1 = np.log(y1)
y2 = np.log(y2)
y3 = np.log(y3)

k = y / x
k1 = y1 / x
k2 = y2 / x
k3 = y3 / x

print(np.mean(k))
print(np.mean(k1))
print(np.mean(k2))
print(np.mean(k3))

plt.plot(x, y, color = 'y', label = 'explicit', linestyle = ':')
plt.plot(x, y1, color = 'black', label = 'implcit', linestyle = ':')
plt.plot(x, y2, color = 'b', label = 'recalc', linestyle = ':')
plt.plot(x, y3, color = 'r', label = 'runge_kutta', linestyle = ':')
plt.xlabel("X")
plt.ylabel("Y")
plt.legend()
plt.title("runge_kutta")
plt.grid()
plt.show()


