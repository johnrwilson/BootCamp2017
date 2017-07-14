# newtons_method.py
"""Volume 1B: Newton's Method.
<John Wilson>
<Math 347>
<1/30/2017>
"""
import numpy as np
from scipy import linalg as la
from matplotlib import pyplot as plt
from scipy import optimize as opt

# Problem 1
def Newtons_method(f, x0, Df, iters=15, tol=1e-5, alpha=1):
    """Use Newton's method to approximate a zero of a function.

    Inputs:
        f (function): A function handle. Should represent a function
            from R to R.
        x0 (float): Initial guess.
        Df (function): A function handle. Should represent the derivative
             of f.
        iters (int): Maximum number of iterations before the function
            returns. Defaults to 15.
        tol (float): The function returns when the difference between
            successive approximations is less than tol.

    Returns:
        A tuple (x, converged, numiters) with
            x (float): the approximation for a zero of f.
            converged (bool): a Boolean telling whether Newton's method
                converged.
            numiters (int): the number of iterations computed.
    """
    n = 0
    converged = False
    if np.isscalar(x0):
        xold = x0+2
        xnew = x0
        while n < iters and np.abs(xnew - xold) > tol:
            xold = xnew
            xnew = xold - alpha*float(f(xold))/Df(xold)
            n += 1
        if abs(xnew - xold) < tol:
            converged=True
    else:
        n = len(x0)
        xold = np.copy(x0) + 1
        xnew = np.copy(x0)
        while la.norm(f(xold)-f(xnew))>tol and n < iters:
            xold = np.copy(xnew)
            xnew -= alpha * la.solve(Df(xnew),f(xnew))
            n+=1
        if la.norm(f(xold)-f(xnew))<tol:
            converged = True
    return xnew, converged, n

def test_prob1(f,df):
    print("Testing problem 1...")
    x, conv, n = Newtons_method(f, 1, df)
    print("The algorithm converged to x={} in {} iterations.".format(x, n))


# Problem 2
def prob2():
    """Given P1[(1+r)**N1-1] = P2[1-(1+r)**(-N2)], if N1 = 30, N2 = 20,
    P1 = 2000, and P2 = 8000, use Newton's method to determine r.
    Return r.
    """
    f = lambda x: (1+x)**30 + 4*(1+x)**-20 - 5
    Df = lambda x: 30*(1+x)**29 - 80*(1+x)**-21
    return Newtons_method(f,.05,Df)[0]
    
def test_prob2():
    print("Testing problem 2...")
    x = prob2()
    print("Problem 2 calculates r = {}".format(x))
    


# Problem 4: Modify Newtons_method and implement this function.
def prob3(alpha):
    """Find an alpha < 1 so that running Newtons_method() on f(x) = x**(1/3)
    with x0 = .01 converges. Return the complete results of Newtons_method().
    """
    f = lambda x: np.sign(x)*np.power(np.abs(x), 1./3)
    Df = lambda x: np.power(np.abs(x), 2./3)**-1*1./3
    return Newtons_method(f,.01,Df,alpha=alpha)
    
def test_prob3():
    print("Testing problem 3...")
    x, conv, n = prob3(1)
    print("With alpha = 1, the method converged: {}".format(conv))
    x, conv, n = prob3(.4)
    print("With alpha = 0.4, the method converged to: {}".format(x))



def prob4(plot = False):
    """Plot effectiveness of newtons_method for different values of alpha.
    Displays a plot of number if iterations v the alpha used.
    """
    f = lambda x: np.sign(x)*np.power(np.abs(x), 1./3)
    Df = lambda x: np.power(np.abs(x), 2./3)**-1*1./3
    avals = []
    nvals = []
    for a in np.linspace(0.01,1,100):
        x, conv, n = Newtons_method(f,0.01,Df,alpha=a,iters = 1000)
        if conv:
            avals.append(a)
            nvals.append(n)
    if plot:
        plt.plot(avals,nvals,label='x^(1/3)')
        plt.legend()
        plt.title("Convergence times for alpha values which actually converged")
        plt.xlabel("Alpha values")
        plt.ylabel("Number of iterations to complete")
        plt.show()
    
# Problem 5
# This problem is implemented above, see Newtons_method function

# Problem 6

def prob6(x0):
    """Solve the system using Newton's method and Newton's method with
    backtracking
    """
    f = lambda x: np.array([4.*x[0]*x[1]-x[0],-1.*x[0]*x[1]+(1-x[1])*(1+x[1])])
    Df = lambda x: np.array([[ 4.*x[1]-1,4.*x[0] ],[ -1.*x[1],-x[0]-2.*x[1] ]])
    a=Newtons_method(f,x0,Df)
    b=Newtons_method(f,x0,Df,alpha=.25)
    return a, b
    
def test_prob6():
    print("Testing problem 6...")
    converged = False
    x0=-.2
    y0=.2
    x, y = prob6(np.array([x0,y0]))
    print("Using ({},{}) as an initial point converged to {} using alpha = 1 and to {}  when alpha = 0.25.".format(x0,y0,x[0],y[0]))
    
    
    
    
# Problem 7
def plot_basins(f, Df, roots, xmin, xmax, ymin, ymax, numpoints=1000, iters=15, colormap='brg'):
    """Plot the basins of attraction of f.

    INPUTS:
        f (function): Should represent a function from C to C.
        Df (function): Should be the derivative of f.
        roots (array): An array of the zeros of f.
        xmin, xmax, ymin, ymax (float,float,float,float): Scalars that define the domain
            for the plot.
        numpoints (int): A scalar that determines the resolution of the plot. Defaults to 100.
        iters (int): Number of times to iterate Newton's method. Defaults to 15.
        colormap (str): A colormap to use in the plot. Defaults to 'brg'.
    """
    xreal=np.linspace(xmin,xmax,numpoints)
    ximag=np.linspace(xmin,xmax,numpoints)
    Xreal,Ximag=np.meshgrid(xreal,ximag)
    Xold=Xreal+1j*Ximag
    Xnew=Xold-f(Xold)/Df(Xold)
    
    for i in range(iters-1):
        Xold=Xnew
        Xnew=Xold-f(Xold)/Df(Xold)
    final=np.zeros_like(Xnew).astype(np.float)
    for i in range(numpoints):
        for j in range(numpoints):
            final[i,j] = np.argmin(np.abs(roots-Xnew[i,j]))
    plt.pcolormesh(Xreal,Ximag,final,cmap=colormap)
    plt.show()
    

def test_prob7():
    """Run plot_basins() on the function f(x) = x^3 - 1 on the domain
    [-1.5,1.5]x[-1.5,1.5].
    """
    print("Testing problem 7...")
    f = lambda x: x**3-1
    Df = lambda x: 3*x**2
    roots=np.array([1,-.5+(3.**.5)*1j/2,-.5-(3.**.5)*1j/2])
    xmax=1.5
    xmin=-1.5
    ymax=1.5
    ymin=-1.5
    plot_basins(f,Df,roots,xmin,xmax,ymin,ymax)
    
    f = lambda x: x**3 - x
    Df = lambda x: 3 * x ** 2 - 1
    roots = np.array([-1., 0, 1.])
    plot_basins(f,Df,roots,xmin,xmax,ymin,ymax)
    
    
    
    
    
if __name__ == "__main__":
    plot = True
    f = lambda x: np.exp(x) - 2
    df = lambda x: np.exp(x)
    test_prob1(f, df)
    test_prob2()
    test_prob3()
    print("Testing problem 4...")
    prob4(plot)
    test_prob6()
    test_prob7()