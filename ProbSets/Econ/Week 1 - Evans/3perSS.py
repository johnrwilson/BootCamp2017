#Solve for the steady-state in 3-per OG model

#Import libraries
import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize as opt

#Household parameters
nvec = np.array([1.0,1.0,0.2])
yrs_live = 60
S = 3.0
beta_annual = 0.96
beta = beta_annual ** (yrs_live / S)
sigma = 3.

#Firm's parameters
alpha = 0.35
A = 1.0
delta_annual = 0.05
delta = 1 - (1.-delta_annual)**(yrs_live/S)
 

def rfunc(K, L, params):
    A, alpha, delta = params
    r = alpha * A * ((L / K) ** (1-alpha)) - delta
    return r

def wfunc(K, L, params):
    A, alpha = params
    w = (1 - alpha) * A * ((K / L) ** alpha)
    return w

def muc1(n, w, r, b2, sigma):
    c1 = n * w - b2
    return c1 ** (-1. * sigma)

def muc2(n, w, r, b2, b3, sigma):
    c2 = n * w + (1. + r) * b2 - b3
    return c2 ** (-1. * sigma)
    
def muc3(n, w, r, b2, b3, sigma):
    c3 = n * w + (1. + r) * b3
    return c3 ** (-1. * sigma)
    
def steady_state(b2b3, args):
    beta, sigma, alpha, A, delta, nvec = args
    b2, b3 = b2b3
    K = np.sum(b2b3)
    L = np.sum(nvec)
    w = wfunc(K, L, (A, alpha))
    r = rfunc(K, L, (A, alpha, delta))
    MUc1 = muc1(nvec[0], w, r, b2, sigma)
    MUc2 = muc2(nvec[1], w, r, b2, b3, sigma)
    MUc3 = muc3(nvec[2], w, r, b2, b3, sigma)
    error1 = MUc1 - beta * (1.+r) * MUc2
    error2 = MUc2 - beta * (1.+r) * MUc3
    return np.array([error1, error2])
    

sol = opt.root(steady_state, np.array([.1,.1]), \
               args = [beta, sigma, alpha, A, delta, nvec])
               
b2, b3 = sol['x']
K = np.sum(sol['x'])
L = np.sum(nvec)
w = (1 - alpha) * A * ((K / L) ** alpha)
r = alpha * A * ((L / K) ** (1-alpha)) - delta
c1 = nvec[0] * w - b2
c2 = nvec[1] * w + (1. + r) * b2 - b3
c3 = nvec[2] * w + (1. + r) * b3
print "c1 = {}".format(c1)
print "c2 = {}".format(c2)
print "c3 = {}".format(c3)
print "b2 = {}".format(b2)
print "b3 = {}".format(b3)
print "w = {}".format(w)
print "r = {}".format(r)


beta = 0.55
sol = opt.root(steady_state, np.array([.1,.1]), \
               args = [beta, sigma, alpha, A, delta, nvec])
               
print sol

