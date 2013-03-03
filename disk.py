#!/usr/bin/env python
#Module to handle the radiative transfer

import numpy as np
import matplotlib.pyplot as plt
import constants as const

class Disk:
    def __init__(self,model,M_star,r0,T_r0,q,Sigma0,p):
        self.model = model #Reference to model
        self.M_star = M_star
        self.r0 = r0 #characteristic radius
        self.T_r0 = T_r0 #Normalization of Temperature
        self.q = q #Exponent of Temperature power law
        self.Sigma0 = Sigma0 #Surface density normalization
        self.p = p #Exponent of surface denisty power law
        self.soft = 0.1*const.AU #softening term for density and temp profiles

    def T(self,r):
        return self.T_r0 * ((r + self.soft)/self.r0)**(-1.*self.q)

    def h(self,r):
        '''Returns in [cm]'''
        return np.sqrt(2. * (r+self.soft)**3 * const.k * self.T(r) / (const.G * self.M_star * const.m0))

    def Sigma(self,r):
        return  self.Sigma0 * ((r+self.soft)/self.r0)**(-1.*self.p)

    def rho0(self,r):
        return self.Sigma(r)/(np.sqrt(2 * np.pi) * self.h(r))

    def rho(self,r,z):
        '''Returns in [g/cm^3]'''
        return self.rho0(r) * np.exp(-1. * z**2. / (2. * self.h(r)**2.))

    def v_phi(self,r):
        return np.sqrt(const.G * self.M_star/ (r+self.soft))

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
        r = np.linspace(0.1, 10*const.AU)#self.r0)
        #r = np.logspace(np.log10(const.AU),np.log10(3 * self.r0))
        #z = np.arange(-0.4, 1.1, 0.01)
        z = np.linspace(0.0, const.AU)
        r_grid, z_grid = np.meshgrid(r, z)
        rho_grid = np.log10(self.rho(r_grid, z_grid))
        levels = np.linspace(-14,-8,num=12)
        CS = ax.contourf(r_grid/const.AU, z_grid/const.AU, rho_grid,levels=levels)
        ax.set_xlabel(r"$r$ (AU)")
        ax.set_ylabel(r"$z$ (AU)")
        ax.set_title(r"$\log_{10}(\rho)$, [$\log_{10}({\rm g/cm}^3)$]" )
        cbar = plt.colorbar(CS)
        #ax.set_xlabel("r [AU]")
        #ax.set_ylabel("T [K]")
        fig.savefig("Plots/Test_Disk/density.png")


def main():
    pass

if __name__=="__main__":
    main()
