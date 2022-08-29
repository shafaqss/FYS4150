import numpy as np
import matplotlib.pyplot as plt

def create_matrix(Nx, r):
    """Creating the matrix for the implicit schemes, however this is
    not needed to solve the equations using the Thomas algorithm for
    tridiaginal matrices
    """
    A = np.zeros((Nx+1, Nx+1))

    for i in range(1, Nx):
        A[i,i-1] = -r
        A[i,i+1] = -r
        A[i,i] = 1 + 2*r
    A[0,0] = A[Nx,Nx] = 1
    return A

def u_analytical(x,t, K=1000):
    """This is the analytical solution to the 1D diffusion equation"""
    k = np.arange(1, K+1, 1)
    X, K = np.meshgrid(x, k)
    ss = np.sum((2/np.pi) * (((-1)**K)/K)*np.exp(-t*(K*np.pi)**2)*np.sin(K*np.pi*X),0) + x
    return ss

def error(analytical, computed):
    """To calculate the error between the actual and computed solution"""
    err = np.abs(analytical - computed)/np.abs(analytical)
    return err
