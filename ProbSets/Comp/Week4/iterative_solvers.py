# iterative_solvers.py
"""Volume 1A: Iterative Solvers.
<John Wilson>
<Math 345>
<10/11/2016>
"""
import numpy as np
from scipy import linalg as la
from scipy import sparse
from matplotlib import pyplot as plt
import time as time
from scipy.sparse import linalg as spla

# Helper function
def diag_dom(n, num_entries=None):
    """Generate a strictly diagonally dominant nxn matrix.

    Inputs:
        n (int): the dimension of the system.
        num_entries (int): the number of nonzero values. Defaults to n^(3/2)-n.

    Returns:
        A ((n,n) ndarray): An nxn strictly diagonally dominant matrix.
    """
    if num_entries is None:
        num_entries = int(n**1.5) - n
    A = np.zeros((n,n))
    rows = np.random.choice(np.arange(0,n), size=num_entries)
    cols = np.random.choice(np.arange(0,n), size=num_entries)
    data = np.random.randint(-4, 4, size=num_entries)
    for i in range(num_entries):
        A[rows[i], cols[i]] = data[i]
    for i in range(n):
        A[i,i] = np.sum(np.abs(A[i])) + 1
    return A


# Problems 1 and 2
def jacobi_method(A, b, tol=1e-8, maxiters=100, plot=False):
    """Calculate the solution to the system Ax = b via the Jacobi Method.

    Inputs:
        A ((n,n) ndarray): A square matrix.
        b ((n,) ndarray): A vector of length n.
        tol (float, opt): the convergence tolerance.
        maxiters (int, opt): the maximum number of iterations to perform.
        plot (bool, opt): if True, plot the convergence rate of the algorithm.
            (this is for Problem 2).

    Returns:
        x ((n,) ndarray): the solution to system Ax = b.
    """
    d=np.diag(A)
    x0=np.zeros_like(d)
    i=0
    diff=1
    e = []
    while (diff>tol and i<maxiters):
        x=x0+(b-np.dot(A,x0))/d
        diff = la.norm(x0-x, ord=np.inf)
        e.append(la.norm(np.dot(A,x)-b, ord=np.inf))
        x0=x
        i+=1
    if plot == True:
        plt.semilogy(range(len(e)),e)
        plt.show()
    return x

def test_prob1(n, plot):
    print("Testing problem 1...")
    A=diag_dom(n)
    b=np.random.rand(n)
    x = jacobi_method(A,b,plot=plot)
    print(np.allclose(A@x,b))
    
    
    
    
# Problem 3
def gauss_seidel(A, b, tol=1e-8, maxiters=100, plot=False):
    """Calculate the solution to the system Ax = b via the Gauss-Seidel Method.

    Inputs:
        A ((n,n) ndarray): A square matrix.
        b ((n,) ndarray): A vector of length n.
        tol (float, opt): the convergence tolerance.
        maxiters (int, opt): the maximum number of iterations to perform.
        plot (bool, opt): if True, plot the convergence rate of the algorithm.

    Returns:
        x ((n,) ndarray): the solution to system Ax = b.
    """
    d=np.diag(A)
    x0=np.zeros_like(d)
    n=0
    diff=1
    e = []
    while (diff>tol and n<maxiters):
        x=[]
        for i in range(d.shape[0]):
            x.append(x0[i]+(1./d[i])*(b[i]-np.dot(np.transpose(A[i]),x0)))
        x=np.array(x)
        diff = la.norm(x0-x, ord=np.inf)
        e.append(la.norm(np.dot(A,x)-b, ord=np.inf))
        x0=x
        n+=1
    if plot == True:
        plt.semilogy(range(len(e)),e)
        plt.show()
    return x

def test_prob3(n,plot):
    print("Testing problem 3...")
    A=diag_dom(n)
    b=np.random.rand(n)
    x = gauss_seidel(A,b,plot=plot)
    print(np.allclose(A@x,b))
    
 
# Problem 4
def sparse_gauss_seidel(A, b, tol=1e-8, maxiters=100):
    """Calculate the solution to the sparse system Ax = b via the Gauss-Seidel
    Method.

    Inputs:
        A ((n,n) csr_matrix): An nxn sparse CSR matrix.
        b ((n,) ndarray): A vector of length n.
        tol (float, opt): the convergence tolerance.
        maxiters (int, opt): the maximum number of iterations to perform.

    Returns:
        x ((n,) ndarray): the solution to system Ax = b.
    """
    if type(A) != sparse.csr_matrix:   
        A=sparse.csr_matrix(A)
    x0=np.zeros_like(b)
    x=np.ones_like(b)
    n=0
    diff=1
    while (diff>tol and n<maxiters):
        for i in range(len(b)):
            rowstart = A.indptr[i]
            rowend = A.indptr[i+1]
            Aix=np.dot(A.data[rowstart:rowend], x[A.indices[rowstart:rowend]])
            x[i] = x0[i] + 1.0 / A[i, i] * (b[i] - Aix)
        diff = la.norm(x0-x, ord=np.inf)
        x0=np.copy(x)
        n+=1
    return x

def test_prob4(n):
    print("Testing problem 4...")
    A = sparse.csr_matrix(diag_dom(n))
    b = np.random.rand(n)
    x = sparse_gauss_seidel(A,b)
    print(np.allclose(A@x,b))
    
# Problem 5
def sparse_sor(A, b, omega, tol=1e-8, maxiters=100):
    """Calculate the solution to the system Ax = b via Successive Over-
    Relaxation.

    Inputs:
        A ((n,n) csr_matrix): An nxn sparse matrix.
        b ((n,) ndarray): A vector of length n.
        omega (float in [0,1]): The relaxation factor.
        tol (float, opt): the convergence tolerance.
        maxiters (int, opt): the maximum number of iterations to perform.

    Returns:
        x ((n,) ndarray): the solution to system Ax = b.
    """
    if type(A) != sparse.csr_matrix:   
        A=sparse.csr_matrix(A)
    x0=np.zeros_like(b)
    x=np.ones_like(b)
    n=0
    diff=1
    while (diff>tol and n<maxiters):
        for i in range(len(b)):
            rowstart = A.indptr[i]
            rowend = A.indptr[i+1]
            Aix=np.dot(A.data[rowstart:rowend], x[A.indices[rowstart:rowend]])
            x[i] = x0[i] + omega / A[i, i] * (b[i] - Aix)
        diff = la.norm(x0-x, ord=np.inf)
        x0=np.copy(x)
        n+=1
    return x, n

def test_prob5(n):
    print("Testing problem 5...")
    A = sparse.csr_matrix(diag_dom(n))
    b = np.random.rand(n)
    x = sparse_sor(A,b,1)[0]
    print(np.allclose(A@x,b))
    
def finite_difference(n):
    """Return the A and b described in the finite difference problem that
    solves Laplace's equation.
    """
    B=np.diag([-4]*n)+np.diag([1]*(n-1),-1) + np.diag([1]*(n-1),1)
    A=sparse.block_diag([B]*n)
    A.setdiag(1,-n)
    A.setdiag(1,n)
    
    b=np.zeros(n**2)
    c=np.arange(0,n**2,n)
    c0=np.arange(n-1,n**2,n)
    
    b[c]=-100
    b[c0]=-100
    
    return A.tocsr(),b

# Problem 6
def prob6(n, omega, tol = 1e-8, maxiters = 100, plot = False):
    A, b = finite_difference(n)
    u, its = sparse_sor(A, b, 1.0, tol=tol, maxiters = maxiters)
    if plot:
        n = int(np.sqrt(len(u)))
        x = np.linspace(0,10,n)
        y = np.linspace(0,10,n)
        X,Y=np.meshgrid(x,y)
        U = u.reshape((n,n))
        plt.pcolormesh(X,Y,U, cmap = "coolwarm")
        plt.colorbar()
        plt.show()
    return u, its
    
def test_prob6(plot):
    print("Testing problem 6...")
    u, n = prob6(20, 1.75, tol=1e-2, maxiters = 1000, plot=plot)
    A, b = finite_difference(20)
    print(np.amax(A@u - b))
    

# Problem 7
def compare_omega():
    """Time sparse_sor() with omega = 1, 1.05, 1.1, ..., 1.9, 1.95, tol=1e-2,
    and maxiters = 1000 using the A and b generated by finite_difference()
    with n = 20. Plot the times as a function of omega.
    """
    A,b=finite_difference(20)
    dom = np.linspace(1,1.95,20)
    iters=[]
    for omega in dom:
        u, n = sparse_sor(A,b,omega,tol=1e-2,maxiters=1000)
        iters.append(n)
    plt.plot(dom,iters)
    plt.show()
    return dom[np.argmin(iters)]

if __name__ == "__main__":
    plot=True
    test_prob1(100,plot)
    test_prob3(100,plot)
    test_prob4(500)
    test_prob5(500)
    test_prob6(plot)
    print("Testing problem 7...")
    o = compare_omega()
    print("Optimal omega is {}".format(o))