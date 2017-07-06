# qr_lstsq_eigs.py
"""Volume 1A: QR 2 (Least Squares and Computing Eigenvalues).
<Name>
<Class>
<Date>"""

import numpy as np
from cmath import sqrt
from scipy import linalg as la
from matplotlib import pyplot as plt


# Problem 1
def least_squares(A, b):
    """Calculate the least squares solutions to Ax = b using QR decomposition.

    Inputs:
        A ((m,n) ndarray): A matrix of rank n <= m.
        b ((m, ) ndarray): A vector of length m.

    Returns:
        x ((n, ) ndarray): The solution to the normal equation.
    """
    Q,R=la.qr(A,mode="economic")
    return la.solve_triangular(R,np.dot(Q.T,b))
    

# Problem 2
def line_fit():
    """Load the data from housing.npy. Use least squares to calculate the line
    that best relates height to weight.

    Plot the original data points and the least squares line together.
    """
    data=np.load("housing.npy")
    year=data[:,0]
    value=data[:,1]
    one=np.ones_like(year)
    X=np.column_stack((year,one))
    x=least_squares(X,value)
    pred=(x[0]*year+x[1])
    
    plt.plot(year,value,'k*')
    plt.plot(year,pred)
        
    
    plt.show()
    


# Problem 3
def polynomial_fit():
    """Load the data from housing.npy. Use least squares to calculate
    the polynomials of degree 3, 6, 9, and 12 that best fit the data.

    Plot the original data points and each least squares polynomial together
    in individual subplots.
    """
    data=np.load("housing.npy")
    
    value=data[:,1]
    yrs=data[:,0]
    
    p3=np.vander(data[:,0],3)
    p6=np.vander(data[:,0],6)
    p9=np.vander(data[:,0],9)
    p12=np.vander(data[:,0],12)
    
    ls3=la.lstsq(p3,value)[0]
    ls6=la.lstsq(p6,value)[0]
    ls9=la.lstsq(p9,value)[0]
    ls12=la.lstsq(p12,value)[0]
    
    x=np.linspace(0,16,100)
    
    plt.subplot(221)
    plt.plot(yrs,np.dot(p3,ls3))
    plt.plot(yrs,value,'k*')

    plt.subplot(222)
    plt.plot(yrs,np.dot(p6,ls6))
    plt.plot(yrs,value,'k*')
    
    plt.subplot(223)
    plt.plot(yrs,np.dot(p9,ls9))
    plt.plot(yrs,value,'k*')

    plt.subplot(224)
    plt.plot(yrs,np.dot(p12,ls12))
    plt.plot(yrs,value,'k*')    
    
    plt.show()
    

def plot_ellipse(a, b, c, d, e):
    """Plot an ellipse of the form ax^2 + bx + cxy + dy + ey^2 = 1."""
    theta = np.linspace(0, 2*np.pi, 200)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    A = a*(cos_t**2) + c*cos_t*sin_t + e*(sin_t**2)
    B = b*cos_t + d*sin_t
    r = (-B + np.sqrt(B**2 + 4*A))/(2*A)

    plt.plot(r*cos_t, r*sin_t, lw=2)
    plt.gca().set_aspect("equal", "datalim")

# Problem 4
def ellipse_fit():
    """Load the data from ellipse.npy. Use least squares to calculate the
    ellipse that best fits the data.

    Plot the original data points and the least squares ellipse together.
    """
    data=np.load("ellipse.npy")
    x=data[:,0]
    y=data[:,1]
    A=np.column_stack((x**2,x,x*y,y,y**2))
    
    ls=la.lstsq(A,np.ones_like(x))[0]
    plot_ellipse(ls[0],ls[1],ls[2],ls[3],ls[4])
    plt.plot(x,y,'k*')
    plt.show()



# Problem 5
def power_method(A, N=20, tol=1e-12):
    """Compute the dominant eigenvalue of A and a corresponding eigenvector
    via the power method.

    Inputs:
        A ((n,n) ndarray): A square matrix.
        N (int): The maximum number of iterations.
        tol (float): The stopping tolerance.

    Returns:
        (foat): The dominant eigenvalue of A.
        ((n, ) ndarray): An eigenvector corresponding to the dominant
            eigenvalue of A.
    """
    m,n=np.shape(A)
    x0=np.random.random(n)
    x0/=la.norm(x0)
    i=0
    diff=tol+1
    while (i<N and diff>tol):
        x1=np.dot(A,x0)
        x1/=la.norm(x1)
        diff=la.norm(x1-x0, ord=np.inf)
        i+=1
        x0=np.copy(x1)
    return np.dot(x0.T,np.dot(A,x0)), x0

# Problem 6
def qr_algorithm(A, N=50, tol=1e-12):
    """Compute the eigenvalues of A via the QR algorithm.

    Inputs:
        A ((n,n) ndarray): A square matrix.
        N (int): The number of iterations to run the QR algorithm.
        tol (float): The threshold value for determining if a diagonal block
            is 1x1 or 2x2.

    Returns:
        ((n, ) ndarray): The eigenvalues of A.
    """
    m,n=A.shape
    S = la.hessenberg(A)
    for k in xrange(N):
        Q,R=la.qr(S)
        S=np.dot(R,Q)
    eigs=[]
    i=0
    while i < n:
        if i == n-1:
            eigs.append(S[i,i])
        elif abs(S[i+1,i])<tol:
            eigs.append(S[i,i])
         
        else:
            a,b,c,d=S[i:i+2,i:i+2].ravel()
            B=-1*(a+d)
            C=a*d-b*c
            D=np.sqrt(B**2-4*C)
            eigs+=[(-B+D)/2.,(-B-D)/2.]
            i+=1
        i+=1
            
    return np.sort(eigs)
    
if __name__=="__main__":
    A=np.random.random((10,10))
    A=A+A.T
    e=qr_algorithm(A)
    