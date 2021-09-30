import numpy as np
import matplotlib.pyplot as plt

def middle(k):
    return sum(k) / len(k)


_x, _y, _y1 = 'dx.txt', 'error_1.txt', 'error_1_1.txt'

x = np.log(np.genfromtxt(_x, dtype = 'double'))
y = np.log(np.genfromtxt(_y, dtype = 'double'))
y1 = np.log(np.genfromtxt(_y1, dtype = 'double'))

k1 = []
k2 = []

for i in range(0, len(x) - 1):
    k1.append((y[i+1] - y[i]) / (x[i+1] - x[i])) 
    k2.append((y1[i+1] - y1[i]) / (x[i+1] - x[i]))

print("k1 = ", middle(k1))
print("k2 = ", middle(k2))

xnumbers = np.linspace(0, 0.2, 5)
ynumbers = np.linspace(0, 0.2, 5)

plt.plot(x, y, color = 'b', label = "table sin' O(h^2)", marker = 'o', linestyle = '--')
plt.plot(x, y1, color='r', label="table sin' O(h)", marker = 'o', linestyle = '--')
plt.xlabel("ln(X)")
plt.ylabel("ln(E)")
plt.title("Approximate")
plt.xticks(xnumbers)
plt.yticks(ynumbers)
plt.legend()
plt.grid()
plt.axis([-8, 0, -20, 0]) # [xstart, xend, ystart, yend]
plt.show()


