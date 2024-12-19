import matplotlib.pyplot as plt
import numpy as np
import math
import sys


with open("data") as file:
	inf		= file.readlines()
first 		= inf[0].split()
second		= inf[1].split()
third		= inf[2].split()
fourth		= inf[3].split()
fifth		= inf[4].split()

a 			= float(first[0])
b 			= float(first[1])
M 			= int(first[2])

x1 			= np.array([float(second[i]) for i in range(0, len(second), 2)])
y1		 	= np.array([float(second[i]) for i in range(1, len(second), 2)])
x2 			= np.array([float(third[i]) for i in range(0, len(third), 2)])
y2		 	= np.array([float(third[i]) for i in range(1, len(third), 2)])
x3 			= np.array([float(fourth[i]) for i in range(0, len(fourth), 2)])
y3		 	= np.array([float(fourth[i]) for i in range(1, len(fourth), 2)])
x4 			= np.array([float(fifth[i]) for i in range(0, len(fifth), 2)])
y4		 	= np.array([float(fifth[i]) for i in range(1, len(fifth), 2)])



# дальше только графика
fig, ax 	= plt.subplots()

ax.plot(x1, y1, "o", label="h")
ax.plot(x2, y2, "o", label="h/2")
ax.plot(x3, y3, "o", label="Runge")
ax.plot(x4, y4, "k", linewidth=2.0, label="D_f(x)")
ax.legend()
ax.set(xlim=(a, b))
ax.grid()
ax.minorticks_on()
plt.show()

