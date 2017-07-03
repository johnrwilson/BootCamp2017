# matplotlib_intro.py
"""Introductory Labs: Intro to Matplotlib.
<John Wilson>
<Math 345>
<9/2/16>
"""

import numpy as np
from matplotlib import pyplot as plt

def var_of_means(n):
    """Construct a random matrix A with values drawn from the standard normal
    distribution. Calculate the mean value of each row, then calculate the
    variance of these means. Return the variance.

    Inputs:
        n (int): The number of rows and columns in the matrix A.

    Returns:
        (float) The variance of the means of each row.
    """
    y = np.random.randn(n,n)
    x = np.mean(y, axis=1)
    return np.var(x)


def prob1():
    """Create an array of the results of var_of_means() with inputs
    n = 100, 200, ..., 1000. Plot and show the resulting array.
    """
    x = list()
    #generate a list of the means
    for i in xrange(100, 1001, 100):
        x.append(var_of_means(i))
    y = np.array(x)
    plt.plot(y)
    plt.show()


def prob2():
    """Plot the functions sin(x), cos(x), and arctan(x) on the domain
    [-2pi, 2pi]. Make sure the domain is refined enough to produce a figure
    with good resolution.
    """
    x = np.linspace(-2*np.pi,2*np.pi,50)
    y1 = np.sin(x)
    y2 = np.cos(x)
    y3 = np.arctan(x)
    plt.plot(x,y1)
    plt.plot(x,y2)
    plt.plot(x,y3)
    plt.show()


def prob3():
    """Plot the curve f(x) = 1/(x-1) on the domain [-2,6].
        1. Split the domain so that the curve looks discontinuous.
        2. Plot both curves with a thick, dashed magenta line.
        3. Change the range of the y-axis to [-6,6].
    """
    #split the domain across the discontinuity
    x1 = np.linspace(-2,.9999,20)
    x2 = np.linspace(1.0001,6,30)
    plt.plot(x1, (1/(x1-1)), "m--", lw = 6)
    plt.plot(x2, (1/(x2-1)), "m--", lw = 6)
    plt.ylim(-6,6)
    plt.show()


def prob4():
    """Plot the functions sin(x), sin(2x), 2sin(x), and 2sin(2x) on the
    domain [0, 2pi].
        1. Arrange the plots in a square grid of four subplots.
        2. Set the limits of each subplot to [0, 2pi]x[-2, 2].
        3. Give each subplot an appropriate title.
        4. Give the overall figure a title.
        5. Use the following line colors and styles.
              sin(x): green solid line.
             sin(2x): red dashed line.
             2sin(x): blue dashed line.
            2sin(2x): magenta dotted line.
    """
    x = np.linspace(0,2*np.pi,50)
    x2 = 2*x
    
    plt.subplot(221)
    plt.plot(x,np.sin(x),"g-")
    plt.title("y=sin(x)")
    plt.axis([0,2*np.pi,-2,2])
    
    plt.subplot(222)
    plt.plot(x,np.sin(x2),"r--")
    plt.title("y=sin(2x)")
    plt.axis([0,2*np.pi,-2,2])
    
    plt.subplot(223)
    plt.plot(x,2*np.sin(x),"b--")
    plt.title("y=2sin(x)")
    plt.axis([0,2*np.pi,-2,2])
    
    plt.subplot(224)
    plt.plot(x,2*np.sin(x2),"m:")
    plt.title("y=2sin(2x)")
    plt.axis([0,2*np.pi,-2,2])
    
    plt.suptitle("Variations on sin(x)")
    plt.show()

def prob5():
    """Visualize the data in FARS.npy. Use np.load() to load the data, then
    create a single figure with two subplots:
        1. A scatter plot of longitudes against latitudes. Because of the
            large number of data points, use black pixel markers (use "k,"
            as the third argument to plt.plot()). Label both axes.
        2. A histogram of the hours of the day, with one bin per hour.
            Label and set the limits of the x-axis.
    """
    df = np.load('FARS.npy')
    
    plt.subplot(121)
    plt.plot(df[:,1],df[:,2],"k,")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.gca().set_aspect("equal")
    
    plt.subplot(122)
    plt.hist(df[:,0], bins=np.arange(0,25))
    plt.xlim(0,24)
    plt.xlabel("Time of the incident")
    
    plt.show()
    
    


def prob6():
    """Plot the function f(x,y) = sin(x)sin(y)/xy on the domain
    [-2pi, 2pi]x[-2pi, 2pi].
        1. Create 2 subplots: one with a heat map of f, and one with a contour
            map of f. Choose an appropriate number of level curves, or specify
            the curves yourself.
        3. Choose a non-default color scheme.
        4. Add a colorbar to each subplot.
    """
    x = np.linspace(-2*np.pi,2*np.pi)
    y = x.copy()
    X,Y = np.meshgrid(x,y)
    Z = np.sin(X) * np.sin(Y) / ( X * Y)
    
    plt.subplot(121)
    plt.pcolormesh(X,Y,Z,cmap="magma")
    plt.axis([-2*np.pi,2*np.pi,-2*np.pi,2*np.pi])
    plt.colorbar()
    
    plt.subplot(122)
    plt.contour(X,Y,Z,10,cmap="viridis")
    plt.colorbar()
    
    plt.show()
    