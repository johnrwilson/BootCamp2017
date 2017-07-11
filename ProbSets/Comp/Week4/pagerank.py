# pagerank.py
"""Volume 1: The Page Rank Algorithm.
<John Wilson>
<Math 347>
<3/27/2017>
"""
import scipy.sparse as sp
import numpy as np
from scipy import linalg as la

# Problem 1
def to_matrix(filename, n):
    """Return the nxn adjacency matrix described by datafile.

    Parameters:
        datafile (str): The name of a .txt file describing a directed graph.
        Lines describing edges should have the form '<from node>\t<to node>\n'.
        The file may also include comments.
    n (int): The number of nodes in the graph described by datafile

    Returns:
        A SciPy sparse dok_matrix.
    """
    A = sp.dok_matrix((n,n))
    with open(filename,'r') as myfile:
        for line in myfile:
            line = line.strip().split()
            try:
                if  A[int(line[0]),int(line[1])] == 0:
                    A[int(line[0]),int(line[1])]=1
            except:
                pass
    return A


# Problem 2
def calculateK(A,N):
    """Compute the matrix K as described in the lab.

    Parameters:
        A (ndarray): adjacency matrix of an array
        N (int): the datasize of the array

    Returns:
        K (ndarray)
    """
    Am = np.copy(A)
    for i in xrange(len(A)):
        if not Am[i].any():
            for j in xrange(len(Am[i])):
                A[i,j]+=1
    D = np.zeros_like(Am)
    for i in xrange(len(Am)):
        D[i,i] = np.sum(A[i])
    return la.solve(D,A).T

def sparse_gauss_seidel(A, b, stol=1e-8, maxiters=1000):
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
    if type(A) != sp.csr_matrix:   
        A=sp.csr_matrix(A)
    x0=np.random.rand(len(b))
    x=np.ones_like(b)
    n=0
    diff=1
    while (diff>stol and n<maxiters):
        for i in xrange(len(b)):
            rowstart = A.indptr[i]
            rowend = A.indptr[i+1]
            Aix=np.dot(A.data[rowstart:rowend], x[A.indices[rowstart:rowend]])
            x[i] = x0[i] + 1.0 / A[i, i] * (b[i] - Aix)
        diff = la.norm(x0-x)
        x0=np.copy(x)
        n+=1
    return x
    
# Problem 3
def iter_solve(adj, N=None, d=.85, tol=1E-5):
    """Return the page ranks of the network described by 'adj'.
    Iterate through the PageRank algorithm until the error is less than 'tol'.

    Parameters:
        adj (ndarray): The adjacency matrix of a directed graph.
        N (int): Restrict the computation to the first 'N' nodes of the graph.
            If N is None (default), use the entire matrix.
        d (float): The damping factor, a float between 0 and 1.
        tol (float): Stop iterating when the change in approximations to the
            solution is less than 'tol'.

    Returns:
        The approximation to the steady state.
    """
    if N:
        K = calculateK(adj[:N-1,:N-1],N)
    else:
        A = np.copy(adj)
        N = len(A)
        K = calculateK(A,N)
    b = ((1-d)/N)*np.ones(N)
    M = np.eye(N)-d*K
    return sparse_gauss_seidel(M, b, stol=tol)


# Problem 4
def eig_solve(adj, N=None, d=.85):
    """Return the page ranks of the network described by 'adj'. Use SciPy's
    eigenvalue solver to calculate the steady state of the PageRank algorithm

    Parameters:
        adj (ndarray): The adjacency matrix of a directed graph.
        N (int): Restrict the computation to the first 'N' nodes of the graph.
            If N is None (default), use the entire matrix.
        d (float): The damping factor, a float between 0 and 1.
        tol (float): Stop iterating when the change in approximations to the
            solution is less than 'tol'.

    Returns:
        The approximation to the steady state.
    """
    if N:
        K = calculateK(adj[:N-1,:N-1],N)
    else:
        A = np.copy(adj)
        N = len(A)
        K = calculateK(A,N)
    B = d*K+float((1-d)/N)*np.ones((N,N))
    E = la.eig(B)
    i = np.argmax(E[0])
    return 1./np.sum(E[1][:,i])*E[1][:,i]


# Problem 5
def team_rank(filename='ncaa2013.csv'):
    """Use iter_solve() to predict the rankings of the teams in the given
    dataset of games. The dataset should have two columns, representing
    winning and losing teams. Each row represents a game, with the winner on
    the left, loser on the right. Parse this data to create the adjacency
    matrix, and feed this into the solver to predict the team ranks.

    Parameters:
        filename (str): The name of the data file.
    Returns:
        ranks (list): The ranks of the teams from best to worst.
        teams (list): The names of the teams, also from best to worst.
    """
    K = np.loadtxt(filename,delimiter=',',skiprows=1,dtype =str)
    d={}
    d_inv={}
    i=0
    for line in K:
        if line[0] not in d.keys():
            d[line[0]]=i
            i+=1
        if line[1] not in d.keys():
            d[line[1]]=i
            i+=1
    d_inv = {v: k for k, v in d.iteritems()}
    n = len(d.keys())
    A = np.zeros((n,n))
    for line in K:
        try:
            A[d[line[1]],d[line[0]]]=1
        except:
            pass
    ranks = eig_solve(A,d=0.7)
    r = np.argsort(ranks)
    best=[]
    for ind in r:
        best.append(d_inv[ind])
    return np.sort(ranks)[::-1].astype(float),best[::-1]
    

if __name__ == "__main__":
    A = to_matrix("matrix.txt",8).toarray()
    #print calculateK(A,8).toarray()
    #x= iter_solve(A)
    #print eig_solve(A)
    print team_rank()