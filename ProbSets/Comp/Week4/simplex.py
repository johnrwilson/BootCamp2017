# simplex.py
"""Volume 2B: Simplex.
<John Wilson>
<Math 323>
<1/18/2017>

Problems 1-6 give instructions on how to build the SimplexSolver class.
The grader will test your class by solving various linear optimization
problems and will only call the constructor and the solve() methods directly.
Write good docstrings for each of your class methods and comment your code.

prob7() will also be tested directly.
"""

import numpy as np

# Problems 1-6
class SimplexSolver(object):
    """Class for solving the standard linear optimization problem

                        maximize        c^Tx
                        subject to      Ax <= b
                                         x >= 0
    via the Simplex algorithm.
    """

    def __init__(self, c, A, b):
        """

        Parameters:
            c (1xn ndarray): The coefficients of the linear objective function.
            A (mxn ndarray): The constraint coefficients matrix.
            b (1xm ndarray): The constraint vector.

        Raises:
            ValueError: if the given system is infeasible at the origin.
        """
        d = b < 0
        if d.any():
            raise ValueError("The system is infeasible at the origin!")
        self.A = np.array(A)
        self.c = c
        self.b = b
        
        m,n = A.shape
        self.m,self.n=m,n
        self.L = list(range(n,n+m)+range(n))
        
        self.initialize_tableau()
        
    def initialize_tableau(self):
        self.Abar = np.bmat([self.A,np.eye(self.m)])
        self.cbar = np.hstack((np.hstack(self.c),np.hstack(np.zeros((self.m,1)))))
        #np.bmat([[np.reshape(self.cbar,(1,len(self.cbar)))],[self.Abar]]) 
        self.T = np.bmat([[np.reshape(np.array([0]),(1,1)),-1.*np.reshape(self.cbar,(1,len(self.cbar))),np.array([[1]])],[np.reshape(self.b,(len(self.b),1)),self.Abar,np.zeros((len(self.b),1))]])
        
    def find_pivot(self):
        col = 1
        row = 1
        while self.T[0,col] >= 0:
            col+=1
        best = np.inf
        mask = self.T[:,col]>0
        if not mask.any():
            raise ValueError("The problem is unbounded")
        for num, thing in enumerate(self.T[1:,col]):
            if mask[num+1]:
                if self.T[num+1,0]/float(thing)<best:
                    best = self.T[num+1,0]/float(thing)
                    row = num+1
        self.col = col
        self.row = row
    
    def pivot(self):
        ind1 = self.row-1
        ind2 = self.L.index(self.col-1)
        self.L[ind1],self.L[ind2] = self.L[ind2],self.L[ind1]
        self.T[self.row,:]/=float(self.T[self.row,self.col])
        for i in xrange(len(self.T)):
            if i == self.row:
                pass
            else:
                factor = self.T[i,self.col]/float(self.T[self.row,self.col])
                self.T[i,:] -= factor * self.T[self.row]
    
    def solve(self):
        """Solve the linear optimization problem.

        Returns:
            (float) The maximum value of the objective function.
            (dict): The basic variables and their values.
            (dict): The nonbasic variables and their values.
        """
        finished = False
        while not finished:
            mask = self.T[0,1:]>=0
            if mask.all():
                finished = True
                break
            self.find_pivot()
            self.pivot()
        optimum = self.T[0,0]
        basic = {}
        for i in xrange(self.m):
            basic[self.L[i]] = self.T[i+1,0]
        non_basic={}
        for i in xrange(self.m,len(self.L)):
            non_basic[self.L[i]] = 0
        return optimum,basic,non_basic


# Problem 7
def prob7(filename='productMix.npz'):
    """Solve the product mix problem for the data in 'productMix.npz'.

    Parameters:
        filename (str): the path to the data file.

    Returns:
        The minimizer of the problem (as an array).
    """
    data = np.load(filename)
    A = data['A']
    p = data['p']
    m = data['m']
    d = data['d']
    A = np.bmat([[A],[np.eye(A.shape[1])]])
    b = np.concatenate((m,d))
    sol = SimplexSolver(p,A,b)
    basic = sol.solve()[1]
    answer =[]
    for i in xrange(4):
        answer.append(basic[i])
    return answer
    
if __name__ == "__main__":
    """
    c = np.array([3,2])
    A = np.array([[1,-1],[3,1],[4,3]])
    b = np.array([2,5,7])
    test = SimplexSolver(c,A,b)
    """