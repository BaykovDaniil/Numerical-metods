import matplotlib.pyplot as plt
import numpy as np
import math
import sys


with open("data") as file:
	inf     = file.readlines()
first 		= inf[0].split()
second		= inf[1].split()
a 			= float(first[0])
b 			= float(first[1])
M 			= int(first[2])
N 			= int(first[3])
dots 		= second[:2*M]  # узлы аппроксимации
data 		= second[2*M:]  # данный о x, многочлене и функции в этих точках

x 			= np.array([float(data[i]) for i in range(0, len(data), 3)])
y_lagrange 	= np.array([float(data[i]) for i in range(1, len(data), 3)])
y_true 		= np.array([float(data[i]) for i in range(2, len(data), 3)])
x_dots 		= np.array([float(dots[i]) for i in range(0, len(dots), 2)])
y_dots 		= np.array([float(dots[i]) for i in range(1, len(dots), 2)])
x_special	= np.array([x_dots[i] for i in range(0, len(x_dots), N-1)])
y_special	= np.array([y_dots[i] for i in range(0, len(y_dots), N-1)])

# дальше только графика
fig, ax 	= plt.subplots()

# y_max 		= math.ceil(max(y_lagrange.max(), y_true.max()))
# y_min 		= min(y_lagrange.min(), y_true.min()) // 1

ax.plot(x, y_lagrange, "k:", linewidth=2.0, label="L(x)")
ax.plot(x, y_true, "k", linewidth=2.0, label="f(x)")
ax.plot(x_dots, y_dots, "o")
ax.plot(x_special, y_special, "o")
ax.legend()
ax.set(xlim=(a//1, math.ceil(b)))
ax.grid()
ax.minorticks_on()
plt.show()
