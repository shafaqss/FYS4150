import numpy as np
import matplotlib.pyplot as plt
from Tridiagonal import*
from common import*

def explicit1D(Nx, Nt, r, IC=0, Lb=0, Rb=1):
    """Implementing the explicit scheme in 1D to solve the diffusion equation
    Nx = number of points in x direction
    Nt = number of points in time
    r = r is dt/(dx**2)
    IC = initial condition u(x,0), defaults to 0
    Lb = left boundary value, u(0,t), defaults to 0
    Rb = rightboundary value, u(L,t) (here we have u(1,t)), defualts to 1
     """
    u   = np.zeros(Nx+1)
    u_new = np.zeros(Nx+1)

    #set initial conditions
    for i in range(0, Nx+1):
        u_new[i] = IC

    for n in range(0, Nt):
        for i in range(1, Nx):
            u[i] = (1-2*r)*u_new[i] + r*u_new[i-1] + r*u_new[i+1]
        #setting boundary conditions
        u[0] = Lb;  u[Nx] = Rb
        #updating the solution before next step
        u_new, u = u, u_new

    return u_new

def implicit1D(Nx, Nt, r, IC=0, Lb=0, Rb=1):
    """Implemetning the implcit scheme in 1D to solve the diffusion equation
    Note = The arguments are exactly the same as the previous function.
    """
    u   = np.zeros(Nx+1)  #solution at current time level
    u_new = np.zeros(Nx+1)  #solution at previous time level
    b = np.zeros(Nx+1)
    A = create_matrix(Nx, r)

    #set initial conditions
    for i in range(0, Nx+1):
        u_new[i] = IC

    for n in range(0, Nt):
        for i in range(1, Nx):
            b[i] = u_new[i]

        b[0] = Lb;  b[Nx] = Rb #setting boundary conditions
        #rr = forward_substitution(A, b)
        #solved = backward_substitution(A, rr)
        u[:] = tridiag(r, b, Nx)
        #uncomment line below, and comment line above,to see correct calculations
        #u[:] = np.linalg.solve(A, b)

        #updating the solution before next step
        u_new, u = u, u_new

    return u_new

def crank_nicolson1D(Nx, Nt, r, IC=0, Lb=0, Rb=1):
    """Implemetning the Crank Nicholson scheme in 1D to solve the diffusion equation
    Note = The arguments are exactly the same as the previous function.
    """
    u   = np.zeros(Nx+1)
    u_new = np.zeros(Nx+1)
    #A = np.zeros((Nx+1, Nx+1))
    rho = r/2
    b = np.zeros(Nx+1)
    A = create_matrix(Nx, rho)

    #set initial conditions
    for i in range(0, Nx+1):
        u_new[i] = IC

    #calculating the right hand side
    for n in range(0, Nt):
        b[1:-1] = u_new[1:-1] + rho*(u_new[:-2] - 2*u_new[1:-1] + u_new[2:])

        #setting boundary conditions
        b[0] = Lb;  b[-1] = Rb

        #solving the linear system of equations now
        u[:] = tridiag(rho, b, Nx)
        #uncomment line below, and comment line above,to see correct calculations
        #u[:] = np.linalg.solve(A, b)

        #updating the solution before next step
        u_new, u = u, u_new

    return u_new

#How to use the program

#defining the parameters
Nx = 10 #position
Nt = 100 #time
T = 0.1 #Time horizon
x = np.linspace(0, 1, Nx+1)
dx = x[1]-x[0]
t = np.linspace(0, T, Nt)
dt = t[1]-t[0]
#v = np.zeros([M, N+1]) #Explicit
r = dt/(dx**2)

print("Information about the scheme")
print("Nx is:", Nx, " and Nt is : ", Nt)
print("Time horizon is T : ", T)
print("The value for dt  : ", dt)
print("The vlaue for dx  : ", dx)
print("The value for r   : ", r)

u_explicit = explicit1D(Nx, Nt, r)
u_implicit = implicit1D(Nx, Nt, r)
u_crank = crank_nicolson1D(Nx, Nt, r)

#error(analytical, computed)
error_ex = error(u_analytical(x,T), u_explicit)
error_im = error(u_analytical(x,T), u_implicit)
error_cr = error(u_analytical(x,T), u_crank)

#relative error plot
plt.figure()
plt.plot(x, error_ex, label='Explicit')
plt.plot(x, error_im, label='Implicit')
plt.plot(x, error_cr, label='Crank Nicolson')
plt.title("Error of the numerical solution in 1D")
#plt.xlim([0, 1.1])
#plt.ylim([0, 1.1])
plt.legend()

plt.figure()
plt.plot(x, u_explicit, label='Explicit')
plt.plot(x, u_implicit, label='Implicit')
plt.plot(x, u_crank, label='Crank Nicolson')
plt.plot(x, u_analytical(x,T), label='Analytic')
#plt.xlim([0, 1.1])
plt.ylim([0, 1.1])
plt.title("Solution of diffusion equation in 1D")
plt.legend()

plt.show()
