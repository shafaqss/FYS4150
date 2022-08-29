## Project 5 - Diffusion equation

Information about scripts:

* diffusion1D.py - contains the explicit scheme, implicit scheme and the crank nicholson scheme for the diffusion equation in 1D. My implicit scheme, and my crank nicholson scheme do not work so perfectly. When I am trying to solve the triadiagonal matrix equation, my function that implements the thomas algorithm for the tridiagonal matrix, seems to have something wrong. The script ```Tridiagonal.py``` containd my function for perfoming solving the matrix system. I have really tried to find the error, I have tried different implementations, however I was unable to completly get a correct version. Therefore to get a correct answer from the implicit and crank nicholson scheme, I have used the numpy function ```np.linalg.solve(A, b)``` to solve the tridiagonal matrix system. Using this I get a correct answer. Just for the sake of doing th stability analysis, I used this function.

* diffusion2D.py - contains the explicit scheme for the diffusion equation in 2D.

* common.py - contains some functions needed in the script ```diffusion1D.py```.

* Tridiagonal.py - contains the function for the solving the tridiagonal matrix system for the 1D implicit and crank nicholson scheme. However I was not able to get the function to work properly. The function that gives me the most accurate solution is the ```tridiag(rho, b, Nx)``` function. You can see that I have really tried to get the tridiagonal function to work properly. This function is needed in the script ```diffusion1D.py```.

* ```FYSproject5shafaq.pdf``` is the final report for project 5.
