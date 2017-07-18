import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from scipy import integrate
from scipy import linalg as la
from sympy import prime


"""
Grading may be most easily accomplished by running the piece of code
which can be found under the if __name__ == "__main__" block.
First change test to True and then run the file in ipython, if you want.
"""



# Problem 1
def prob1(g,a,b,N,method):
    if a == b:
        return 0
    if a > b:
        a, b = b, a
    if method == "midpoint":
        points = np.linspace(0,N-1,N)
        points = (points * 2 + 1) * float(b - a) / (2 * N) + a
        return np.sum(g(points)) * float(b - a) / N
    elif method == "trapezoid":
        points = np.linspace(0,N,N+1)
        points = points * float(b - a) / N + a
        inside = g(points[0]) + 2. * np.sum(g(points[1:-1])) + g(points[-1])
        return float(b - a) * inside / (2. * N)
    elif method == "Simpsons":
        points = np.linspace(0,2*N,2*N+1)
        points = points * float(b - a) / (2 * N) + a
        inside = g(points[0]) + 4. * np.sum(g(points[1::2])) + 2. * np.sum(g(points[2:-2:2])) + g(points[-1])
        return inside / (6. * (N )) * (b-a)
    else:
        raise ValueError("Integration method must equal 'midpoint', 'trapezoid', or 'Simpsons'")
        
def test_prob1():
    g = lambda x: 0.1*x**4 - 1.5*x**3 + 0.53*x**2 + 2.*x + 1.
    true = 4373.33
    x = prob1(g,-10,10,1000,'midpoint')
    print("Midpoint error: {}".format(abs(x-true)))
    x = prob1(g,-10,10,1000,'trapezoid')
    print("Trapezoid error: {}".format(abs(x-true)))
    x = prob1(g,-10,10,300,'Simpsons')#[0]
    print("Simpsons error: {}".format(abs(x-true)))
    
    
# Problem 2
def normal_dist(mu,sigma,N,k):
    int_len = 2 * k * sigma
    node_len = int_len / (N - 1)
    nodes = np.zeros(N, dtype=float)
    low_node = mu - k * sigma
    for i in range(N):
        nodes[i] = low_node + i * node_len
    cdf = stats.norm.cdf
    weights = np.zeros(N)
    for i in range(N):
        if i == 0:
            weights[i] = cdf(.5 * (nodes[i] + nodes[i+1]), loc = mu, scale = sigma)
        elif i == N-1:
            weights[i] = 1 - cdf(.5 * (nodes[i] + nodes[i-1]), loc = mu, scale = sigma)
        else:
            weights[i] = cdf(.5 * (nodes[i] + nodes[i+1]), loc = mu, scale = sigma) -\
                                cdf(.5 * (nodes[i] + nodes[i-1]), loc = mu, scale = sigma)
    return nodes, weights
    

def test_prob2(printing = False):
    if printing:
        print(normal_dist(0, 1, 100, 3))
    else:
        normal_dist(0, 1, 100, 3)
        print("No errors!")
    
    
# Problem 3
def lognorm(mu, sigma, N, k):
    nodes, weights = normal_dist(mu, sigma, N, k)
    nodes = np.exp(nodes)
    return nodes, weights
    

def test_prob3(printing = False):
    if printing:
        print(lognorm(0, 1, 100, 3))
    else:
        lognorm(0, 1, 100, 3)
        print("No errors!")
    
    
# Problem 4
def prob4():
    n, w = lognorm(10.5, 0.8, 100, 4)
    expected_value = np.dot(n, w)
    actual_value = np.exp(10.5 + (0.8)**2 / 2)
    print("The computed expected value is {}.".format(expected_value))
    print("The actual expected value is {}.".format(actual_value))
    print("There was a relative error of {}.".format(abs(expected_value - actual_value) / actual_value))
    
def test_prob4():
    prob4()

# Problem 5
# First several helper functions are defined
# For the theory behind this method, see the ACME page on Gaussian Quadrature
# Note that I wrote this code myself :)
def shift(f, a, b, plot=False):
    """Shift the function f on [a, b] to a new function g on [-1, 1] such that
    the integral of f from a to b is equal to the integral of g from -1 to 1.

    Inputs:
        f (function): a scalar-valued function on the reals.
        a (int): the left endpoint of the interval of integration.
        b (int): the right endpoint of the interval of integration.
        plot (bool): if True, plot f over [a,b] and g over [-1,1] in separate
            subplots.

    Returns:
        The new, shifted function.
    """
    g = lambda x: f((b-a)/2.*x+(a+b)/2.)
    
    if plot == True:
        domg = np.linspace(-1,1,200)
        domf = np.linspace(a,b,200)
        
        plt.subplot(211)
        plt.plot(domf,f(domf))
        plt.title('f(x)')
        
        plt.subplot(212)
        plt.plot(domg,g(domg))
        plt.title('g(x)')
            
        plt.show()
        
    return g
    
def estimate_integral(f, a, b, points, weights):
    """Estimate the value of the integral of the function f over [a,b].

    Inputs:
        f (function): a scalar-valued function on the reals.
        a (int): the left endpoint of the interval of integration.
        b (int): the right endpoint of the interval of integration.
        points ((n,) ndarray): an array of n sample points.
        weights ((n,) ndarray): an array of n weights.

    Returns:
        The approximate integral of f over [a,b].
    """
    s = (b-a)/2.
    g = shift(f,a,b)
    return s*np.dot(weights.T,g(points))

def construct_jacobi(gamma, alpha, beta):
    """Construct the Jacobi matrix."""
    a = -1. * beta / alpha
    l = len(beta)
    b = np.zeros(l-1)
    for i in range(l-1):
        b[i] = ((gamma[i+1])/(alpha[i]*alpha[i+1]))**.5
    J = np.diag(a) + np.diag(b,-1) + np.diag(b,1)
    return J

def points_and_weights(n):
    """Calculate the points and weights for a quadrature over [a,b] with n
    points.

    Returns:
        points ((n,) ndarray): an array of n sample points.
        weights ((n,) ndarray): an array of n weights.
    """
    g = np.zeros(n)
    a = np.zeros(n)
    b = np.zeros(n)
    for k in range(n):
        g[k] = float(k)/(k+1)
        a[k] = (2.*(k+1)-1)/(k+1)
    J = construct_jacobi(g,a,b)
    #print J
    x,v = la.eig(J)
    w=2.*(v[0]**2)
    return x,w
    
def gaussian_quadrature(f, a, b, n):
    """Using the functions from the previous problems, integrate the function
    'f' over the domain [a,b] using 'n' points in the quadrature.
    """
    points, weights = points_and_weights(n)
    return estimate_integral(f,a,b,points,weights)    
    
def prob5():
    g = lambda x: 0.1*x**4 - 1.5*x**3 + 0.53*x**2 + 2.*x + 1.
    return gaussian_quadrature(g, -10, 10, 3)
    
def test_prob5():
    x = prob5()
    print("Gaussian quadrature yielded an approximation of {}.".format(x))
    true = 4373.33
    print("The error is {}.".format(abs(true - x)))
    
def prob6():
    g = lambda x: 0.1*x**4 - 1.5*x**3 + 0.53*x**2 + 2.*x + 1.
    return integrate.quad(g, -10, 10)
    
def test_prob6():
    x = prob6()
    print("Scipy's Gaussian quadrature yielded an approximation of {}.".format(x))
    true = 4373.33
    print("The error is {}.".format(abs(true - x[0])))
    
def prob7(g, xbounds, ybounds, N):
    xmin, xmax = xbounds
    ymin, ymax = ybounds
    vol = (xmax - xmin) * (ymax - ymin)
    xs = np.random.uniform(xmin, xmax, N)
    ys = np.random.uniform(ymin, ymax, N)
    val = 0
    for i in range(N):
        val += g([xs[i],ys[i]])
    return vol * val / N
    
def h(pts):
    x = pts[0]
    y = pts[1]
    if x ** 2 + y ** 2 <= 1:
        return 1
    else:
        return 0

def test_prob7():
    diff = 1
    N=1
    while diff > 1e-5:
        approx = prob7(h, [-1,1], [-1,1], N)
        diff = abs(approx - 3.1415)
        N += 1
    print("It took N={} before the approximation was accurate to 4 decimals.".format(N))
    print("The approximated value is {}.".format(approx))
    
def prob8(n, d, method):
    points = np.zeros((n, d))
    if method == "Weyl":
        for i in range(n):
            for j in range(d):
                point = (i + 1) * np.sqrt(prime(j + 1))
                points[i,j] = point - int(point)
    elif method == "Haber":
        for i in range(n):
            for j in range(d):
                point = (i + 1) * (i + 2) / 2 * np.sqrt(prime(j + 1))
                points[i,j] = point - int(point)
    elif method == "Niederreiter":
        for i in range(n):
            for j in range(d):
                point = (i + 1) * 2 ** (float(j + 1) / (j + 2))
                points[i,j] = point - int(point)
    elif method == "Baker":
        for i in range(n):
            for j in range(d):
                point = (i + 1) * np.exp(prime(j + 1))
                points[i,j] = point - int(point)
    return points
    
def test_prob8(plot):
    if plot:
        w = prob8(1000, 2, "Weyl")
        wx = w[:,0]
        wy = w[:,1]
        h = prob8(1000, 2, "Haber")
        hx = h[:,0]
        hy = h[:,1]
        n = prob8(1000, 2, "Niederreiter")
        nx = n[:,0]
        ny = n[:,1]
        b = prob8(1000, 2, "Baker")
        bx = b[:,0]
        by = b[:,1]
        plt.subplot(221)
        plt.plot(wx, wy, ',')
        plt.title("Weyl")
        plt.subplot(222)
        plt.plot(hx, hy, ',')
        plt.title("Haber")
        plt.subplot(223)
        plt.plot(nx, ny, ',')
        plt.title("Niederreiter")
        plt.subplot(224)
        plt.plot(bx, by, ',')
        plt.title("Baker")
        plt.show()
    
    
def prob9(g, N):
    x = prob8(N, 2, "Weyl")
    x = x * 2 - 1
    xs = x[:,0]
    ys = x[:,1]
    val = 0
    vol = 4
    for i in range(N):
        val += g([xs[i],ys[i]])
    return vol * val / N


def test_prob9():
    diff = 1
    N=1
    while diff > 1e-5:
        approx = prob9(h, N)
        diff = abs(approx - 3.1415)
        N += 1
    print("It took N={} before the approximation was accurate to 4 decimals.".format(N))
    print("The approximated value is {}.".format(approx))
    
    
    
if __name__ == "__main__":
    test = False
    if test:
        print("Testing Problem 1")
        test_prob1()
        print("Testing Problem 2")
        test_prob2()
        print("Testing Problem 3")
        test_prob3()
        print("Testing Problem 4")
        test_prob4()
        print("Testing Problem 5")
        test_prob5()
        print("Testing Problem 6")
        test_prob6()
        print("Testing Problem 7")
        test_prob7()
        print("Testing Problem 8")
        test_prob8(True)
        print("Testing Problem 9")
        test_prob9()
    
    pass