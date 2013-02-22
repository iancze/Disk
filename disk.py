#!/usr/bin/env python
#Module to handle the radiative transfer

import numpy as np
import matplotlib.pyplot as plt
import constants as const

class Disk:
    def __init__(self,M_star,r0,T_r0,q,Sigma0,p,m0,X0):
        self.M_star = M_star
        self.r0 = r0 #characteristic radius
        self.T_r0 = T_r0 #Normalization of Temperature
        self.q = q #Exponent of Temperature power law
        self.Sigma0 = Sigma0 #Surface density normalization
        self.p = p #Exponent of surface denisty power law
        self.m0 = m0 #mean molecular weight of gas
        self.X0 = X0 #gas fraction of CO

    def T(self,r):
        return self.T_r0 * (r/self.r0)**(-1.*self.q)

    def h(self,r):
        return np.sqrt(2. * r**3 * const.k * self.T(r) / (const.G * self.M_star * self.m0))

    def Sigma(self,r):
        return  self.Sigma0 * (r/self.r0)**(-1.*self.p)

    def rho0(self,r):
        return self.Sigma(r)/(np.sqrt(2 * np.pi) * self.h(r))

    def rho(self,r,z):
        return self.rho0(r) * np.exp(-1. * z**2. / (2. * self.h(r)**2.))

    def plot_T(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        rs = np.logspace(np.log10(const.AU),np.log10(3 *self.r0))
        ax.loglog(rs/const.AU,self.T(rs))
        ax.set_xlabel("r [AU]")
        ax.set_ylabel("T [K]")
        fig.savefig("Plots/Test_Disk/temperature_log.png")

    def plot_h(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        rs = np.logspace(np.log10(const.AU),np.log10(3 *self.r0))
        ax.loglog(rs/const.AU,self.h(rs)/const.AU)
        ax.set_xlabel("r [AU]")
        ax.set_ylabel("h [AU]")
        fig.savefig("Plots/Test_Disk/h_log.png")

    def plot_rho(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        r = np.linspace(const.AU, 3 * self.r0)
        #r = np.logspace(np.log10(const.AU),np.log10(3 * self.r0))
        #z = np.arange(-0.4, 1.1, 0.01)
        z = np.linspace(const.AU, 100.*const.AU)
        r_grid, z_grid = np.meshgrid(r, z)
        rho_grid = np.log10(self.rho(r_grid, z_grid))
        CS = ax.contour(r_grid/const.AU, z_grid/const.AU, rho_grid,20)
        cbar = plt.colorbar(CS)
        #ax.set_xlabel("r [AU]")
        #ax.set_ylabel("T [K]")
        fig.savefig("Plots/Test_Disk/density.png")


def main():
    pass

if __name__=="__main__":
    main()
