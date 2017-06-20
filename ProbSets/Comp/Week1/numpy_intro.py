# numpy_intro.py
"""Introductory Labs: Intro to NumPy.
<John Wilson>
<MATH 345>
<9/6/2016>
"""
import numpy as np

def prob1():
	"""Define the matrices A and B as arrays. Return the matrix product AB."""
	A = np.array( [ [3, -1, 4],[1, 5, -9] ] )
	B = np.array( [ [2, 6, -5, 3],[5, -8, 9, 7],[9, -3, -2, -3] ] )
	C = np.dot(A,B)
	return C


def prob2():
    """Define the matrix A as an array. Return the matrix -A^3 + 9A^2 - 15A."""
    A = np.array( [ [3, 1, 4],[1, 5, 9],[-5, 3, 1] ] )
    A2 = np.dot(A,A)
    A3 = np.dot(A2,A)
    return (-1 * A3) + (9 * A2) - (15 * A)



def prob3():
    """Define the matrices A and B as arrays. Calculate the matrix product ABA,
    change its data type to np.int64, and return it.
    """
    A = np.triu(np.ones(7))
    B = np.tril(np.full_like(A, -6))+5
    C = np.dot(np.dot(A,B),A)
    C = C.astype(np.int64)
    return C
    

def prob4(A):
    """Make a copy of 'A' and set all negative entries of the copy to 0.
    Return the copy.

    Example:
        >>> A = np.array([-3,-1,3])
        >>> prob4(A)
        array([0, 0, 3])
    """
    B = np.copy(A)
    mask = B <= 0
    B[mask] = 0
    return B
    


def prob5():
    """Define the matrices A, B, and C as arrays. Return the block matrix
                                | 0 A^T I |
                                | A  0  0 |,
                                | B  0  C |
    where I is the identity matrix of appropriate size and each 0 is a matrix
    of all zeros, also of appropriate sizes.
    """
    A = np.array( [ [0,2,4],[1,3,5] ] )
    B = np.tril(np.full((3,3),3,dtype=np.int))
    C = np.triu(np.tril(np.full((3,3),-2,np.int)))
    D = np.hstack((np.zeros_like(B),A.T,np.eye(3)))
    E = np.hstack((A,np.zeros((2,2),dtype=np.int),np.zeros((2,3),dtype=np.int)))
    F = np.hstack((B,np.zeros((3,2),dtype=np.int),C))
    return np.vstack((D,E,F))


def prob6(A):
    """Divide each row of 'A' by the row sum and return the resulting array.

    Example:
        >>> A = np.array([[1,1,0],[0,1,0],[1,1,1]])
        >>> prob6(A)
        array([[ 0.5       ,  0.5       ,  0.        ],
               [ 0.        ,  1.        ,  0.        ],
               [ 0.33333333,  0.33333333,  0.33333333]])
    """
    B = np.sum(A,axis=1,dtype=np.float)
    return np.divide(A,B.reshape((len(A),1)))


def prob7():
    """Given the array stored in grid.npy, return the greatest product of four
    adjacent numbers in the same direction (up, down, left, right, or
    diagonally) in the grid.
    """
    grid = np.load("grid.npy")
    horz_grid = grid[:,:-3] * grid[:,1:-2] * grid[:,2:-1] * grid[:,3:]
    vert_grid = grid[:-3,:] * grid[1:-2,:] * grid[2:-1,:] * grid[3:,:]
    right_diag = grid[:-3,:-3] * grid[1:-2,1:-2] * grid[2:-1,2:-1] * grid[3:,3:]
    left_diag = grid[:-3,3:] * grid[1:-2,2:-1] * grid[2:-1,1:-2] * grid[3:,:-3]
    return max(np.max(horz_grid),np.max(vert_grid),np.max(right_diag),np.max(left_diag))
	
