# qr_decomposition.py
"""Volume 1A: QR 1 (Decomposition).
<John Wilson>
<Math 345>
<10/18/2016>
"""

import numpy as np
from scipy import linalg as la


# Problem 1
def qr_gram_schmidt(A):
    """Compute the reduced QR decomposition of A via Modified Gram-Schmidt.

    Inputs:
        A ((m,n) ndarray): A matrix of rank n.

    Returns:
        Q ((m,n) ndarray): An orthonormal matrix.
        R ((n,n) ndarray): An upper triangular matrix.
    """
    m,n = A.shape
    Q = np.copy(A)
    R = np.zeros((n,n))
    for i in xrange(n):
        R[i,i] = la.norm(Q[:,i])
        Q[:,i] /= R[i,i]
        for j in xrange(i+1,n):
            R[i,j] = np.dot(Q[:,j].T,Q[:,i])
            Q[:,j] -= R[i,j]*Q[:,i]
    return Q, R


# Problem 2
def abs_det(A):
    """Use the QR decomposition to efficiently compute the absolute value of
    the determinant of A.

    Inputs:
        A ((n,n) ndarray): A square matrix.

    Returns:
        (float) the absolute value of the detetminant of A.
    """
    Q,R=qr_gram_schmidt(A)
    determinant=R[0,0]
    for i in xrange(1,R.shape[0]):
        determinant *= R[i,i]
    return abs(determinant)


# Problem 3
def solve(A, b):
    """Use the QR decomposition to efficiently solve the system Ax = b.

    Inputs:
        A ((n,n) ndarray): An invertible matrix.
        b ((n, ) ndarray): A vector of length n.

    Returns:
        x ((n, ) ndarray): The solution to the system Ax = b.
    """
    Q,R=qr_gram_schmidt(A)
    y=np.dot(Q.T,b)
    n=A.shape[0]
    x=np.zeros(n)
    for k in reversed(xrange(n)):
        x[k]=(1./R[k,k])*(y[k]-np.dot(R[k,k:],x[k:]))
    return x
    

# Problem 4
def qr_householder(A):
    """Compute the full QR decomposition of A via Householder reflections.

    Inputs:
        A ((m,n) ndarray): A matrix of rank n.

    Returns:
        Q ((m,m) ndarray): An orthonormal matrix.
        R ((m,n) ndarray): An upper triangular matrix.
    """
    sign = lambda x: 1 if x >= 0 else -1
    m,n=A.shape
    R=np.copy(A)
    Q=np.eye(m)
    for k in xrange(n):
        u = np.copy(R[k:,k])
        u[0] += sign(u[0])*la.norm(u)
        u /= la.norm(u)
        R[k:,k:] -= np.outer(2*u, np.dot(u.T,R[k:,k:]))
        Q[k:,:] -= np.outer(2*u, np.dot(u.T,Q[k:,:]))
        
    return Q.T, R

# Problem 5
def hessenberg(A):
    """Compute the Hessenberg form H of A, along with the orthonormal matrix Q
    such that A = QHQ^T.

    Inputs:
        A ((n,n) ndarray): An invertible matrix.

    Returns:
        H ((n,n) ndarray): The upper Hessenberg form of A.
        Q ((n,n) ndarray): An orthonormal matrix.
    """
    sign = lambda x: 1 if x >= 0 else -1
    m,n=A.shape
    H=np.copy(A).astype(np.float)
    Q=np.eye(m)
    for k in xrange(n-2):
        u=np.copy(H[k+1:,k])
        u[0] += sign(u[0])*la.norm(u)
        u /= la.norm(u)
        H[k+1:,k:]-= 2*np.outer(u, np.dot(u.T,H[k+1:,k:]))
        H[:,k+1:] -= 2*np.outer(np.dot(H[:,k+1:],u),u.T)
        Q[k+1:,:] -= np.outer(2*u, np.dot(u.T,Q[k+1:,:]))
    return H, Q.T
        
if __name__=="__main__":
    A=np.random.random((8,8))
    H,Q=hessenberg(A)
    np.allclose(np.triu(H,-1),H)
    np.allclose(np.dot(np.dot(Q,H),Q.T),A)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
