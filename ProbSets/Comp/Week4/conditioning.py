# conditioning.py
"""Volume 1B: Conditioning.
<John Wilson>
<Math 347>
<2/7/17>
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg as la
from sympy import subfactorial as sub
from sympy import integrate as integ
from sympy import exp
from sympy.abc import x
import math
from astropy.table import Table



#Problem 1
def prob1(A):
    """Take in a matrix A and compute its condition number using singular
    value decomposition. The largest singular value is divided by the smallest.
    """
    svds = la.svdvals(A)
    k = np.argsort(svds)
    if svds[k[0]] == 0:
        return np.inf
    else:
        return svds[k[-1]]/svds[k[0]]


# Problem 2
def prob2():
    """Randomly perturb w_coeff by replacing each coefficient a_i with
    a_i*r_i, where r_i is drawn from a normal distribution centered at 1 with
    standard deviation 1e-10.

    Plot the roots of 100 such experiments in a single graphic, along with the
    roots of the unperturbed polynomial w(x).

    Using the final experiment only, estimate the relative and absolute
    condition number (in any norm you prefer).

    Returns:
        Display a graph of all 100 perturbations.
        Print the values of relative and absolute condition numbers.
    """
    w_roots = np.arange(1, 21)
    w_coeffs = np.array([1, -210, 20615, -1256850, 53327946, -1672280820,
                        40171771630, -756111184500, 11310276995381,
                        -135585182899530, 1307535010540395,
                        -10142299865511450, 63030812099294896,
                        -311333643161390640, 1206647803780373360,
                        -3599979517947607200, 8037811822645051776,
                        -12870931245150988800, 13803759753640704000,
                        -8752948036761600000, 2432902008176640000])
    reals=np.arange(1,21)
    imags=np.zeros(20)
    for i in range(200):
        tempc=np.copy(w_coeffs)
        r=np.random.normal(1,1e-10,21)
        tempc*=r
        reals=np.append(reals,np.roots(np.poly1d(tempc)).real)
        imags=np.append(imags,np.roots(np.poly1d(tempc)).imag)
    plt.plot(reals[:20],imags[:20],'*')
    plt.plot(reals[20:],imags[20:],',')
    plt.show()
    pt_roots=np.roots(np.poly1d(tempc))
    k=la.norm(pt_roots-w_roots,np.inf)/la.norm(r,np.inf)
    k_r=k*la.norm(w_coeffs,np.inf)/la.norm(w_roots,np.inf)
    return k,k_r


# Problem 3
def eig_condit(M):
    """Approximate the condition number of the eigenvalue problem at M.

    Inputs:
        M ((n,n) ndarray): A square matrix.

    Returns:
        (float) absolute condition number of the eigenvalue problem at M.
        (float) relative condition number of the eigenvalue problem at M.
    """
    eigs=la.eig(M)[0]
    perturb=np.random.normal(0,1e-10,M.shape)+np.random.normal(0,1e-10,M.shape)*1j
    eigsp = la.eig(M+perturb)[0]
    k=la.norm(eigs-eigsp)/la.norm(perturb)
    return k,k*la.norm(M)/la.norm(eigs)

# Problem 4
def plot_eig_condit(x0=-100, x1=100, y0=-100, y1=100, res=200):
    """Create a grid [x0, x1] x [y0, y1] with the given resolution. For each
    entry (x,y) in the grid, find the relative condition number of the
    eigenvalue problem, using the matrix   [[1 x]
                                            [y 1]]  as the input.
    Use plt.pcolormesh() to plot the condition number over the entire grid.

    Inputs:
        x0 (float): min x-value.
        x1 (float): max x-value.
        y0 (float): min y-value.
        y1 (float): max y-value.
        res (int): number of points along each edge of the grid.
    """
    x = np.linspace(x0,x1,res)
    y = np.linspace(y0,y1,res)
    Z = np.empty((res,res))
    for i in range(res):
        for j in range(res):
            A=np.array([[1,x[i]],[y[j],1]])
            Z[i,j] = eig_condit(A)[1]
    X,Y = np.meshgrid(x,y)
    plt.pcolormesh(X,Y,Z,cmap='gray_r')
    plt.colorbar()
    plt.show()
    
# Problem 5            
def prob5(n):
    """
    Accepts an interger n and finds the least square solution to the stability data problem.
    Plots the forward errors.
    """

# Problem 6
def integral(n):
    """Calculate the integral from 0 to 1 of x^n e^{x-1} dx using the closed
    form solution (-1)^n !n + (-1)^{n+1} n!/e.
    """
    return (-1)**n*sub(n)+(-1)**(n+1)*(math.factorial(n)/np.exp(1))

def prob6():
    """For the values of n in the problem, compute integral(n). Compare
    the values to the actual values, and print your explanation of what
    is happening.
    """
    # Actual values of the integral at specified n.
    ns=[1,5,10,15,20,25,30,35,40,45,50]
    actual_values = []
    vals=[]
    rfe=[]
    for n in ns:
        act=float(integ(x**n * exp(x-1),(x,0,1)))
        comp=integral(n)
        vals.append(comp)
        actual_values.append(act)
        rfen=np.abs(act-comp)/np.abs(act)
        rfe.append(rfen)
    plt.plot(ns,rfe)
    plt.yscale('log')
    plt.show()
        
    #rfe=la.norm(actual-computed)/la.norm(actual)
    
    t=Table()
    t['ns']=ns
    t['Actual values']=actual_values
    t['Computed values']=vals
    #print(t)
    #print("There is error coming from the fact that values close together are being subtracted, so the computer can't store them properly because the difference is so small.")
    