import numpy as np
from matplotlib import pyplot as plt


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
        inside = g(points[0]) + 2. * np.sum(g(points[1:-2])) + g(points[-1])
        return float(b - a) * inside / (2 * N)
    elif method == "Simpsons":
        points = np.linspace(0,2*N,2*N+1)
        points = points * float(b - a) / (2 * N) + a
        inside = g(points[0]) + 4. * np.sum(g(points[1::2])) + 2. * np.sum(g(points[2:-2:2])) + g(points[-1])
        return float(b - a) * inside / (3. * (N + 1))
    else:
        raise ValueError("Integration method must equal 'midpoint', 'trapezoid', or 'Simpsons'")
        
def test_prob1():
    g = lambda x: 0.1*x**4 - 1.5*x**3 + 0.53*x**2 + 2.*x + 1.
    true = 4373.33
    x = prob1(g,-10,10,500,'midpoint')
    print("Midpoint error: {}".format(abs(x-true)))
    x = prob1(g,-10,10,500,'trapezoid')
    print("Trapezoid error: {}".format(abs(x-true)))
    x = prob1(g,-10,10,500,'Simpsons')
    print("Simpsons error: {}".format(abs(x-true)))
        
if __name__ == "__main__":
    pass
    