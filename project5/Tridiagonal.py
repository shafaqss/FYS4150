"""This script contains the funstion to perform the thomas alorithm to
do solve the tridiagonal system for the implicit schemes. The only function
that gives me decent results is the function tridiag. I have tried to get the
function to work perfectly, but the best results i could get was from
tridiag"""
import numpy as np

def tridiag(alpha, u, N): #gives best results
    """
    Performs tridiagonal gaussian elimination, here the diagonal element is
    1+2*alpha, and the upper and lower elements are - alpha.
    alpha = value for the elements, 1+2*alpha and -alpha.
    u = right hand side
    N = length of the solution
    """
    d = np.zeros(N) + (1+2*alpha)
    b = np.zeros(N-1) - alpha

    # perform forward elimination
    for i in range(1,N):
        #normalizing row i
        b[i-1] /= d[i-1];
        u[i] /= d[i-1]
        d[i-1] = 1
        #elimination
        u[i+1] += u[i]*alpha
        d[i] += b[i-1]*alpha
    #Normalize bottom row
    u[N] /= d[N-1]
    d[N-1] = 1

    #backward substitution
    for i in range(N, 1, -1):
        u[i-1] -= u[i]*b[i-2]

    return u


def backward_substitution(matrix, b):
    """Implements backward substitution for a general matrix"""
    n,m  = matrix.shape
    x = np.zeros(n)
    x[n-1] = b[n-1] / matrix[n-1,n-1]

    for i in range(n-2,-1,-1):
        sum = 0
        for j in range(i+1,n):
            sum += matrix[i,j]*x[j]
        x[i] = (b[i]-sum)/matrix[i,i]

    return x

def forward_substitution(matrix, b):
    """Implements forward substitution for a general matrix"""
    n,m  = matrix.shape
    x = np.zeros(n)
    x[0] = b[0] / matrix[0,0]
    for i in range(1,n):
        sum = 0
        for j in range(0,i):
            sum += matrix[i,j]*x[j]
        x[i] = (b[i]-sum)/matrix[i,i]

    return x

def solver_tri(b, Nx, r):
    """Solving a linear system where the matrix is tridiagonal"""
    n = len(b)
    v = np.zeros(n)
    alpha = np.ones(Nx+1)*(1+2*r)
    beta = np.ones(Nx+1)*(-r)
    gamma = np.ones(Nx)*(-r)
    alpha[0] = alpha[-1] = 1
    beta[0] = beta[-1] = 0
    gamma[0] = 0

    delta = np.zeros(n)
    c = np.zeros(n)
    m = np.zeros(n)

    for k in range(1,n):
        m[k] = beta[k]/delta[k-1]
        delta[k] = alpha[k] - m[k]*gamma[k-1]
        c[k] = b[k] - m[k]*c[k-1]

    v[n-1] = c[n-1]/delta[n-1]
    for k in range(n-2, -1, -1):
        v[k] = ( c[k] - gamma[k]*v[k+1] )/ delta[k]
    return v


def solving_tri(matrix, b):
    """Solving a system of equations where the matrix is tridiagonal"""
    #forward substitution
    n,m  = matrix.shape
    x = np.zeros(n)
    x[0] = b[0] / matrix[0,0]
    for i in range(1,n):
        sum = 0
        for j in range(0,i):
            sum += matrix[i,j]*x[j]
        x[i] = (b[i]-sum)/matrix[i,i]

    #backward substitute
    x[n-1] = b[n-1] / matrix[n-1,n-1]
    for i in range(n-2,-1,-1):
        sum = 0
        for j in range(i+1,n):
            sum += matrix[i,j]*x[j]
        x[i] = (b[i]-sum)/matrix[i,i]

    return x

def thomas_tridiag(a,b,c,d):
    """Solution of a linear system of equations with tridiagonal matrix.
    a is the lower diagonal, b is th main diagonal, c is the upper
    diagonal, an d is the right hand side.
    """
    n = len(d)
    w= np.zeros(n-1,float)
    g= np.zeros(n, float)
    p = np.zeros(n,float)

    w[0] = c[0]/b[0]
    g[0] = d[0]/b[0]

    for i in range(1,n-1):
        w[i] = c[i]/(b[i] - a[i-1]*w[i-1])
    for i in range(1,n):
        g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1])
    p[n-1] = g[n-1]
    for i in range(n-1,0,-1):
        p[i-1] = g[i-1] - w[i-1]*p[i]
    return p
