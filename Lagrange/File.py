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
N 			= int(first[2])
dots 		= second[:2*N]  # узлы аппроксимации
data 		= second[2*N:]  # данный о x, многочлене и функции в этих точках

x 			= np.array([float(data[i]) for i in range(0, len(data), 3)])
y_lagrange 	= np.array([float(data[i]) for i in range(1, len(data), 3)])
y_true 		= np.array([float(data[i]) for i in range(2, len(data), 3)])
x_dots 		= np.array([float(dots[i]) for i in range(0, len(dots), 2)])
y_dots 		= np.array([float(dots[i]) for i in range(1, len(dots), 2)])

# дальше только графика
fig, ax 	= plt.subplots()

ax.plot(x, y_lagrange, "k:", linewidth=2.0, label="L(x)")
ax.plot(x, y_true, "k", linewidth=2.0, label="f(x)")
ax.plot(x_dots, y_dots, "o")
ax.legend()
ax.set(xlim=(a, b))
ax.grid()
ax.minorticks_on()

plt.show()