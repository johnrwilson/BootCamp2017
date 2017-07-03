import numpy as np
from matplotlib import pyplot as plt
import scipy.linalg as la

#Problem 1

def prob1():
    A = np.array([[0.6,0.1,-0.3],[0.5,-0.4,0.2],[1.0,-0.2,1.1]])
    b = np.array([12.,10.,-1.])

    x0 = np.random.rand(3)

    x1 = np.dot(A,x0)+b

    while np.linalg.norm(x1-x0) > 1e-13:
        x0 = np.copy(x1)
        x1 = np.dot(A,x0)+b
    
    return x1
    
def prob1b():
    A = np.array([[0.6,0.1,-0.3],[0.5,-0.4,0.2],[1.0,-0.2,1.1]])
    b = np.array([12.,10.,-1.])
    return la.solve((A-np.eye(3)),-b)
    
def prob3():
    beta = 0.96
    w = np.array([0.5,1.0,1.5])
    p = np.array([0.2,0.4,0.4])
    w_bars=[]
    cs = np.linspace(1,2,100)
    for c in cs:
        x0=2.0
        x1=1.0
        while abs(x0 - x1) > 1e-13:
            x0 = x1
            x1 = c * (1 - beta) + beta * np.sum(np.maximum(w,x0)*p)
        w_bars.append(x1)
    plt.plot(cs,w_bars)
    plt.xlabel("Compensation value")
    plt.ylabel("Reservation wage")
    plt.savefig('res_wage.png')
    plt.show()