import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from common import*

def analytical_2D(x, y, t):
    """This is the analytical solution to the 2D diffusion equation"""
    return np.sin(np.pi*x)*np.sin(np.pi*y)*np.exp(-np.pi**2*t)

def explicit_2D(v_new, v_old, r, meshpoints, dt, total_T, t):
    """Implementing the explicit scheme in 2D to solve the diffusion equation
    v_new = old solution matrix
    v_old = new solution matrix
    r = r is dt/(h**2)
    meshpoints = number of points in the grid
    dt = spacing for the time
    total_T = total time the simulation runs for
    t = time variable, stop simulating when we reach total_T
     """
    while (t < total_T):
        for i in range(1, meshpoints-1):
            for j in range(1, meshpoints-1):
                v_new[i,j] = (1-4*r)*v_old[i, j] + r*(v_old[i+1, j] + v_old[i-1, j] + v_old[i, j+1] + v_old[i, j-1]);
        v_old = v_new
        t = t + dt

    return v_new, t

#How to use the program

#defining the parameters
r = 0.15 #in two dimensions r=dt/(h**2)
h = 0.1
dt = r*(h**2)
t = 0
T = 10 #total time simulation
xinitial = 0 # starting value for x
xfinal = 1 # final value for x

#equal number of points in x and y direction
points = int(((xfinal - xinitial)/h) + 1)
x = np.zeros(points)
y = np.zeros(points)
# setting up points in x and y directions
for i in range(points):
    x[i] = i*h
    y[i] = i*h

#create empty matrix for the solutions
u_new = np.zeros((points, points))
u_old = np.ones((points, points))

#setting the initial condition
u_old = np.ones((points, points))*( np.sin(np.pi*x)*np.sin(np.pi*y))
#u_old = np.random.uniform(low=28.5, high=55.5, size=(points, points))

#setting boundary conditions
u_old[0, :] = 0
u_old[-1, :] = 0
u_old[:, 0] = 0
u_old[:, -1] = 0

#explicit_2D(v_new, v_old, r, meshpoints, dt, total_T, t)
u_exp2D, time = explicit_2D(u_new, u_old, r, points, dt, T, t)
X,Y = np.meshgrid(x,y)
A = analytical_2D(X, Y, t=time)

print("Shape of matrices")
print("Analytical  ", A.shape)
print("Computed    ", u_exp2D.shape)
print()
print("Information about parameters")
print("The value for r   : ",r)
print("The value for h   : ",h)
print("The value for dt  : ",dt)
print("The time is       : ",time)
print("Time horizon is T : ", T)


plt.figure()
plt.contourf(x, y, u_exp2D, levels = 40, cmap = 'plasma')
plt.xlabel("x", fontsize = 15)
plt.title("Numerical solution in 2D, r = {}".format(r), fontsize = 16)
plt.ylabel("y", fontsize = 15)
cbar = plt.colorbar()
cbar.set_label("u(x,y,t = {})".format(np.round(time)), size = 12)
cbar.ax.tick_params(labelsize = 15)

plt.figure()
plt.contourf(X, Y, A, levels = 40, cmap = 'plasma')
plt.xlabel("x", fontsize = 15)
plt.title("Analytical solution", fontsize = 16)
plt.ylabel("y", fontsize = 15)
cbar = plt.colorbar()
cbar.set_label("u(x,y,t={})".format(np.round(time)), size = 12)
cbar.ax.tick_params(labelsize = 15)

fig = plt.figure()
ax = fig.gca(projection="3d")
surf = ax.plot_surface(x, y, u_exp2D, cmap=cm.Spectral, antialiased=False)
plt.title("Numerical solution of diffusion equation in 2D")
fig.colorbar(surf, shrink=0.5, aspect=5)

fig = plt.figure()
ax = fig.gca(projection="3d")
surf = ax.plot_surface(X,Y,A, cmap=cm.Spectral, antialiased=False)
plt.title("Analytical solution of diffusion equation in 2D")
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
