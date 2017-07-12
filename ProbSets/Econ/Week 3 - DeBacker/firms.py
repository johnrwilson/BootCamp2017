# import packages
import matplotlib.pyplot as plt
import numpy as np
from numba import jit, jitclass
import scipy.integrate as integrate
from scipy.stats import norm
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import optimize as opt


@jit
def find_vmat_fn(V,ez,sizek,z_grid,betafirm):
    Vmat = np.zeros((len(z_grid), sizek, sizek))  # initialize Vmat matrix
    for i in range(sizek):  # loop over k
        for j in range(sizek):  # loop over k'
            for k in range(len(z_grid)):
                Vmat[k, i, j] = ez[k, i, j] + betafirm * V[k][j]
    return Vmat


@jit
def find_valuefunc_fn(ez,VFtol,VFmaxiter,V,PF,z_grid,pi,betafirm,sizek):
    V=np.copy(V)
    PF=np.copy(PF)
    start_time = time.clock()
    Vstore = np.zeros((sizek, VFmaxiter))
    VFdist = 7.0
    VFiter = 1
    while VFdist > VFtol and VFiter < VFmaxiter:
        TV = np.copy(V)
        Vmat = find_vmat_fn(V,ez,sizek,z_grid,betafirm)
        for l in range(len(z_grid)):
            weights = pi[l]
            Vmat1 = np.average(Vmat,axis=0,weights=weights)
            Vstore[:, VFiter] = V[l].reshape(sizek,)  # store value function at each
            # iteration for graphing later
            V[l] = Vmat1.max(axis=1)  # apply max operator to Vmat (to get V(k))
            PF[l] = np.argmax(Vmat1, axis=1)  # find the index of the optimal k'
        VFdist = (np.absolute(V - TV)).max()  # check distance between value
        # function for this iteration and value function from past iteration
        VFiter += 1

    VFI_time = time.clock() - start_time
    if VFiter < VFmaxiter:
        #print('Value function converged after this many iterations:', VFiter)
        pass
    else:
        print('Value function did not converge')
    #print('VFI took ', VFI_time, ' seconds to solve')
    #print("VFdist was ", VFdist)
    return V, PF, Vstore
    
    
    


class firm(object):


    def __init__(self, rho=0.7605, mu=0.0, sigma_eps=0.213, alpha_k=0.297, alpha_l=0.65, \
                delta=0.154, psi=1.08, w=0.7, beta=0.96, z = 1.0, N=9, h=6.616):
        self.rho=rho
        self.mu=mu
        self.sigma_eps=sigma_eps
        self.alpha_k=alpha_k
        self.alpha_l=alpha_l
        self.delta=delta
        self.psi=psi
        self.w=w
        self.r= (1. / beta) - 1
        self.betafirm= beta
        self.VFtol = 1e-6
        self.VFmaxiter = 3000
        self.N=9
        self.z=z
        self.h = h
            
    def make_z(self,N=9):
        sigma_z = self.sigma_eps / ((1 - self.rho ** 2) ** (1 / 2))
        z_cutoffs = (sigma_z * norm.ppf(np.arange(N + 1) / N)) + self.mu

        # compute grid points for z
        z_grid = ((N * sigma_z * (norm.pdf((z_cutoffs[:-1] - self.mu) / sigma_z)
                                      - norm.pdf((z_cutoffs[1:] - self.mu) / sigma_z)))
                      + self.mu)
        z_grid=np.exp(z_grid)

        # define function that we will integrate
        def integrand(x, sigma_z, sigma_eps, rho, mu, z_j, z_jp1):
            val = (np.exp((-1 * ((x - mu) ** 2)) / (2 * (sigma_z ** 2)))
                    * (norm.cdf((z_jp1 - (mu * (1 - rho)) - (rho * x)) / sigma_eps)
                       - norm.cdf((z_j - (mu * (1 - rho)) - (rho * x)) / sigma_eps)))
            
            return val

        # compute transition probabilities
        pi = np.empty((N, N))
        for i in range(N):
            for j in range(N):
                results = integrate.quad(integrand, z_cutoffs[i], z_cutoffs[i + 1],
                                         args = (sigma_z, self.sigma_eps, self.rho, self.mu,
                                                 z_cutoffs[j], z_cutoffs[j + 1]))
                pi[i,j] = (N / np.sqrt(2 * np.pi * sigma_z ** 2)) * results[0]
        self.pi=pi
        self.z_grid=z_grid
        
        
    def make_kgrid(self,dens=1):
        N = self.N
        kstar = ((((1 / self.betafirm - 1 + self.delta) * ((self.w / self.alpha_l) **
                                         (self.alpha_l / (1 - self.alpha_l)))) /
                                         (self.alpha_k * (self.z ** (1 / (1 - self.alpha_l))))) **
                                         ((1 - self.alpha_l) / (self.alpha_k + self.alpha_l - 1)))
        kbar = 2*kstar
        lb_k = 0.001
        ub_k = kbar
        krat = np.log(lb_k / ub_k)
        numb = np.ceil(krat / np.log(1 - self.delta))
        K = np.zeros(int(numb * dens))
        # we'll create in a way where we pin down the upper bound - since
        # the distance will be small near the lower bound, we'll miss that by little
        for j in range(int(numb * dens)):
            K[j] = ub_k * (1 - self.delta) ** (j / dens)
        kvec = K[::-1]
        sizek = kvec.shape[0]
        self.kvec=kvec
        self.sizek = sizek
        self.V = np.zeros((N,sizek))  # initial guess at value function
        self.PF = np.zeros((N,sizek))
        
    def e(self):
        e = np.zeros((len(self.z_grid),self.sizek, self.sizek))
        for h in range(len(self.z_grid)):
            op = ((1 - self.alpha_l) * ((self.alpha_l / self.w) ** (self.alpha_l / (1 - self.alpha_l))) *\
              ((self.kvec ** self.alpha_k) ** (1 / (1 - self.alpha_l)))) *\
              (self.z_grid[h] ** (1 / (1 - self.alpha_l)))
            for i in range(self.sizek):
                for j in range(self.sizek):
                        e[h, i, j] = (op[i] - self.kvec[j] + ((1 - self.delta) * self.kvec[i]) -
                                   ((self.psi / 2) * ((self.kvec[j] - ((1 - self.delta) * self.kvec[i])) ** 2)
                                    / self.kvec[i]))
        self.ez = e

    def find_vmat(self):
        self.Vmat = find_vmat_fn(self.V,self.ez,self.sizek,self.z_grid,self.betafirm)
        
    def find_valuefunc(self, plot = False):
        self.e()
        self.V, self.PF, self.Vstore = find_valuefunc_fn(self.ez,self.VFtol,self.VFmaxiter,self.V,self.PF,\
                            self.z_grid,self.pi,self.betafirm,self.sizek)
        if plot:
            for i in range(len(self.z_grid)):
                plt.plot(self.kvec, self.V[i], label=str(self.z_grid[i]))
            plt.xlabel("Initial capital 'k'")
            plt.ylabel("Value function")
            plt.legend()
            plt.show()
        
        
    def step2(self):
        LD = np.zeros((len(self.z_grid),self.sizek))
        for i in range(len(self.z_grid)):
            for j in range(len(self.kvec)):
                LD[i,j] = (self.alpha_l/self.w) ** (1. / (1 - self.alpha_l)) * \
                        self.z_grid[i] ** (1. / (1 - self.alpha_l)) * self.kvec[j] ** (self.alpha_k / (1 - self.alpha_l))
        self.LD = LD
        
        ID = np.zeros((len(self.z_grid),self.sizek))
        for i in range(len(self.z_grid)):
            for j in range(len(self.kvec)):
                policy = self.PF[i,j]
                kprime = self.kvec[int(policy)]
                ID[i,j] = kprime -  (1 - self.delta) * self.kvec[i]
        self.ID = ID
        
        #If the earnings grid must be flattened into a true zxk grid, do so here
        
    
    def stationary_distribution(self, plot = False):

        # for finding the stationary distribution
        Pi = self.pi
        PF = self.PF  # the firm's policy function: k'(z,k)
        sizez, sizek = PF.shape

        # To plot the stationary distribution
        z = self.z_grid
        kvec = self.kvec # the grid of capital

        '''
        ------------------------------------------------------------------------
        Compute the stationary distribution of firms over (k, z)
        ------------------------------------------------------------------------
        SDtol     = tolerance required for convergence of SD
        SDdist    = distance between last two distributions
        SDiter    = current iteration
        SDmaxiter = maximium iterations allowed to find stationary distribution
        Gamma     = stationary distribution
        HGamma    = operated on stationary distribution
        ------------------------------------------------------------------------
        '''
        Gamma = np.ones((sizez, sizek)) * (1 / (sizek * sizez))
        SDtol = 1e-12
        SDdist = 7
        SDiter = 0
        SDmaxiter = 1000
        while SDdist > SDtol and SDmaxiter > SDiter:
            HGamma = np.zeros((sizez, sizek))
            for i in range(sizez):  # z
                for j in range(sizek):  # k
                    for m in range(sizez):  # z'
                        HGamma[m, int(PF[i, j])] = \
                            HGamma[m, int(PF[i, j])] + Pi[i, m] * Gamma[i, j]
            SDdist = (np.absolute(HGamma - Gamma)).max()
            Gamma = HGamma
            SDiter += 1

        if SDiter < SDmaxiter:
            #print('Stationary distribution converged after this many iterations: ',
            #     SDiter)
            self.Gamma = Gamma
        else:
            print('Stationary distribution did not converge')
        # plot the SD
        if plot:
            # Plot the stationary distribution over k
            fig, ax = plt.subplots()
            ax.plot(kvec, Gamma.sum(axis=0))
            plt.xlabel('Size of Capital Stock')
            plt.ylabel('Density')
            plt.title('Stationary Distribution over Capital')
            plt.show()

            fig, ax = plt.subplots()
            ax.plot(np.log(z), Gamma.sum(axis=1))
            plt.xlabel('Log Productivity')
            plt.ylabel('Density')
            plt.title('Stationary Distribution over Productivity')
            
            
            # Stationary distribution in 3D
            zmat, kmat = np.meshgrid(kvec, np.log(z))
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(kmat, zmat, Gamma, rstride=1, cstride=1, cmap=cm.Blues,
                            linewidth=0, antialiased=False)
            ax.view_init(elev=20., azim=20)  # to rotate plot for better view
            ax.set_xlabel(r'Log Productivity')
            ax.set_ylabel(r'Capital Stock')
            ax.set_zlabel(r'Density')
            plt.show()
        
        
    def find_aggregates(self):
        Gamma = self.Gamma
        LD = self.LD
        ID = self.ID
        LDagg = 0
        Iagg = 0
        PSIagg = 0
        Yagg = 0
        for i in range(len(self.z_grid)):
            for j in range(self.sizek):
                LDagg += (LD[i,j] * Gamma[i,j])
                Iagg += (ID[i,j] * Gamma[i,j])
                PSIagg += ((self.psi / 2.) * (ID[i,j] ** 2) / self.kvec[j])
                Yagg += (self.z_grid[i] * (self.kvec[j] ** self.alpha_k) * (LD[i,j] ** self.alpha_l))
        self.LDagg = LDagg
        self.Iagg = Iagg
        self.PSIagg = PSIagg
        self.Yagg = Yagg
        self.Cagg = Yagg - Iagg - PSIagg
        self.LSagg = -(self.w) / (self.h * self.Cagg)
        
    def check_equality(self):
        return self.LSagg - self.LDagg
        
        
    def run_test(self, w, plot = False):
        self.w = w
        #print("Making discrete z grid...")
        self.make_z()
        #print("Making k grid...")
        self.make_kgrid()
        #print("Finding value function...")
        self.find_valuefunc(plot = plot)
        #print("Finding labor and investment demand...")
        self.step2()
        #print("Calculating the stationary distribution...")
        self.stationary_distribution(plot=plot)
        #print("Finding aggregates...")
        self.find_aggregates()
        #print("Finished running.")
        return self.check_equality()
        
if __name__ == "__main__":
    t = firm()
    print(opt.brentq(t.run_test,0.9,.91))