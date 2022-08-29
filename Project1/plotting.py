import numpy as np
import matplotlib.pyplot as plt

x_values = []
u_values = []

infile = open("resultater.txt","r")
n = int(infile.readline())
for line in infile:
    words = line.split(",")
    x_values.append(float(words[0]))
    u_values.append(float(words[1]))
infile.close()

def analytic_answer(x):
    return 1-(1-np.exp(-10))*x - np.exp(-10*x)

x_values = np.array(x_values)
u_values = np.array(u_values)

plt.plot(x_values,u_values, label="Numerical solution")
plt.plot(x_values,analytic_answer(x_values), label="Exact solution")
plt.title("Numerical solution with n=%d points"%(n))
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()
plt.show()
