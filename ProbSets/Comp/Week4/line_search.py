# one_dimensional_optimization.py
"""Volume 2B: 1-D Optimization.
<John Wilson>
<Math 323>
<2/7/2017>
"""
import math
import numpy as np
import time

# Problem 1
def golden_section(f, a, b, niter=100, tol=1e-3):
    """Find the minimizer of the unimodal function f on the interval [a,b]
    using the golden section search method.

    Inputs:
        f (function): unimodal scalar-valued function on R.
        a (float): left bound of the interval of interest.
        b (float): right bound of the interval of interest.
        niter (int): number of iterations to compute.

    Returns:
        the approximated minimizer (the midpoint of the final interval).
    """
    rho=.5*(3.-math.sqrt(5))
    n=0
    while n<niter and abs(a - b)>tol*2:
        aprime=a+rho*(b-a)
        bprime=a+(1-rho)*(b-a)
        if f(aprime)<=f(bprime):
            b=bprime
        else:
            a=aprime
        n+=1
    if abs(a - b)>tol*2:
        print("Did not converge.")
    else:
        return .5*(a+b), n

def test_prob1(f):
    start = time.clock()
    v, n = golden_section(f, 0, 3)
    print("Golden Section took {} steps to converge in {} seconds.".format(n, time.clock()-start))
    
# Problem 2
def bisection(df, a, b, niter=100, tol=1e-3):
    """Find the minimizer of the unimodal function with derivative df on the
    interval [a,b] using the bisection algorithm.

    Inputs:
        df (function): derivative of a unimodal scalar-valued function on R.
        a (float): left bound of the interval of interest.
        b (float): right bound of the interval of interest.
        niter (int): number of iterations to compute.
    """
    n=0
    while n<niter and np.abs(a - b) > 2*tol:
        mid=.5*(a+b)
        if df(mid)<0:
            a=mid
        else:
            b=mid
        n+=1
    if np.abs(a - b) > 2*tol:
        print("Did not converge.")
    else:
        return .5*(a+b), n

def test_prob2(df):
    start = time.clock()
    v, n = bisection(df, 0, 3)
    print("Bisection took {} steps to converge in {} seconds.".format(n, time.clock()-start))
    
# Problem 3
def newton1d(f, df, ddf, x, niter=100, tol = 1e-3):
    """Minimize the scalar function f with derivative df and second derivative
    df using Newton's method.

    Parameters
        f (function): A twice-differentiable scalar-valued function on R.
        df (function): The first derivative of f.
        ddf (function): The second derivative of f.
        x (float): The initial guess.
        niter (int): number of iterations to compute.

    Returns:
        The approximated minimizer.
    """
    xnew = x
    n = 0
    diff = 1
    while n < niter and diff > tol:
        xold = xnew
        xnew = xold - float(df(xold))/ddf(xold)
        n += 1
        diff = np.abs(xnew - xold)
    if np.abs(xnew - xold) > tol:
        print("Did not converge to minimizer.")
    else:
        return xnew

def test_prob3(f, df, ddf):
    x = newton1d(f, df, ddf, 0)
    print("Newton's method converged to the minimizer {} with value {}.".format(x, f(x)))
    x = newton1d(f, df, ddf, 10)
    print("Newton's mothod did not work for initial guess {}".format(10))
    
# Problem 4
def secant1d(f, df, x0, x1, niter=10, tol = 1e-3):
    """Minimize the scalar function f using the secant method.

    Inputs:
        f (function): A differentiable scalar-valued function on R.
        df (function): The first derivative of f.
        x0 (float): A first initial guess.
        x1 (float): A second initial guess.
        niter (int): number of iterations to compute.

    Returns:
        The approximated minimizer.
    """
    xold = x0
    xnew = x1
    n = 0
    while n < niter and abs(xnew - xold) > tol:
        temp = xnew
        xnew = xnew - df(xnew)*(xnew-xold)/(df(xnew)-df(xold))
        xold = temp
        n += 1
    return xnew

def test_prob4(f, df):
    x = secant1d(f, df, 0, -1)
    print("Newton's method converged to the minimizer {} with value {}.".format(x, f(x)))
    x = secant1d(f, df, 1, 0)
    print("Newton's mothod did not work for initial guess {}, where it converged to {} with value {}".format((1,0),x, f(x)))
    
    
# Problem 5
def backtracking(f, slope, x, p, a=1, rho=.9, c=10e-4):
    """Do a backtracking line search to satisfy the Wolfe Conditions.
    Return the step length.

    Inputs:
        f (function): A scalar-valued function on R.
        slope (float): The derivative of f at x.
        x (float): The current approximation to the minimizer.
        p (float): The current search direction.
        a (float): Initial step length (set to 1 in Newton and quasi-Newton
            methods).
        rho (float): Parameter in (0,1).
        c (float): Parameter in (0,1).

    Returns:
        The computed step size.
    """
    while f(x+a*p)>f(x)+c*a*slope*p:
        a *= rho
        
    return a
    
if __name__ =="__main__":
    f = lambda x: np.exp(x)-4.*x
    df = lambda x: np.exp(x)-4.
    ddf = lambda x: np.exp(x)
    g=lambda x: x**2+np.sin(x)+np.sin(10.*x)
    dg=lambda x: 2.*x + np.cos(x) + 10.*np.cos(10.*x)
    h = lambda x: -x**2 + np.sin(5*x)
    dh = lambda x: -2*x + 5. * np.cos(5*x)
    ddh = lambda x: -2 - 25. * np.sin(5*x)
    print("Testing problem 1")
    test_prob1(f)
    print("Testing problem 2")
    test_prob2(df)
    print("Testing problem 3")
    test_prob3(h,dh,ddh)
    print("Testing problem 4")
    test_prob4(g,dg)