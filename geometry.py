#!/usr/bin/env python

import numpy as np
import constants as const
import radiation as rad
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve

#Testing whether the geometry is correct in determining starting positions on the box.

#Box is total width (in z direction) of 2w and total length (in y direction) of 2l
global w,l
w = 200. * const.AU
l = 500. * const.AU

class Orientation:
    def __init__(self,model,theta,distance):
        self.model = model #reference to model class
        self.theta = theta #inclination. 0 is face-on, pi/2 is edge-on
        self.distance = distance
        self.theta_c = np.arctan(l/w)
        self.delta_c = w * np.sin(self.theta) - l * np.cos(self.theta)
        self.delta_limit = w * np.sin(self.theta) + l * np.cos(self.theta)
        self.slope = 1./np.tan(self.theta)
        pass

    def center_line(self):
        '''return a function z = f(y) that is the line which passes through (y,z)=(0,0) and is parallel to line of sight. Angle with z axis will be theta.'''
        return lambda y: self.slope * y

    def line_through_point(self,yz):
        '''given a point yz, generate a line function (parallel to line of sight) that passes through this point'''
        y0,z0 = yz
        intercept = z0 - self.slope * y0
        return lambda y: self.slope * y + intercept

    def graph_on_ax(self,ax):
        ax.grid(True)
        ax.axhline(0,color="k")
        ax.axvline(0,color="k")
        ax.add_patch(plt.Rectangle((-l,-w),2*l,2*w,fill=False))
        ax.set_xlim(-1.5*l,1.5*l)
        ax.set_ylim(-1.5*w,1.5*w)
        ax.set_xlabel("y")
        ax.set_ylabel("z")
        ys = np.linspace(0,6.0)
        ax.plot(ys,self.center_line()(ys),"k")

    def __str__(self):
        return """
        theta = {theta:.2f}; {theta_deg:.2f}
        theta_c = {theta_c:.2f}, {theta_c_deg:.2f}
        delta_c = {delta_c:.2f}
        delta_high = {delta_high:.2f}
        delta_low = {delta_low:.2f}
        slope = {slope:.2f}""".format(theta=self.theta,theta_deg=(self.theta * 180./np.pi),theta_c_deg=(self.theta_c * 180./np.pi),theta_c=self.theta_c,delta_c=self.delta_c,delta_high=self.delta_limit,delta_low=(-1.*self.delta_limit),slope=self.slope)

class LineOfSight:
    def __init__(self,model,orientation,delta,alpha,nu):
        '''alpha_phys and delta_phys are the actual projected distances at the location of the disk.'''
        self.model = model #Link to the model
        self.orn = orientation #Link to the orientation object
        self.alpha = alpha
        self.alpha_phys = self.alpha * self.orn.distance * (const.AU/const.pc)
        self.delta = delta
        self.delta_phys = self.delta * self.orn.distance * (const.AU/const.pc)
        self.x0 = self.alpha_phys
        self.xf = self.x0
        self.nu = nu
        self.k_nu = rad.k_nu(self.nu,self.model.beta)
        self.calc_y0_z0()
        self.disk = self.orn.model.disk
        self.tau = np.vectorize(self.tau_general) #vectorized tau function
        self.ss = np.arange(0.,self.s_finish,0.1*const.AU)
        self.grid = self.model.grid
        #self.walk_along_grid()

    def calc_y0_z0(self):
        '''Given the inital geometry via the Orientation object, calculate y0, z0, yf, zf, and s_finish.'''
        if self.delta_phys > self.orn.delta_c:
            self.z0 = w
            self.y0 = w * np.tan(self.orn.theta) - self.delta_phys/np.cos(self.orn.theta)
            #k2 is the critical length, see pg 15, 2/21/13
            k2 = (self.y0 + l)/np.tan(self.orn.theta)
            if k2 > 2. * w:
                self.zf = -w
                self.yf = self.y0 - 2 * w * np.tan(self.orn.theta)
            else:
                self.zf = self.z0 - k2
                self.yf = -l
        else:
            self.y0 = l
            self.z0 = self.delta_phys/np.sin(self.orn.theta) + l/np.tan(self.orn.theta)
            k2 = (self.z0 + w) * np.tan(self.orn.theta)
            if k2 > 2. * l:
                self.yf = -l
                self.zf = self.z0 - 2. * l/np.tan(self.orn.theta)
            else:
                self.yf = self.y0 - k2
                self.zf = -w
        self.s_finish = (self.y0 - self.yf)/np.sin(self.orn.theta)

    def walk_along_grid(self):
        x_grid = self.grid.x_grid
        y_grid = self.grid.y_grid
        z_grid = self.grid.z_grid

        #Define initial starting cell
        self.i = np.max(np.where(self.x0 > x_grid)[0])
        self.j = np.max(np.where(self.y0 > y_grid)[0])
        self.k = np.max(np.where(self.z0 > z_grid)[0])

        print("i0:%i; j0:%i; k0:%i" % (self.i,self.j,self.k))
        print()

        s_total = 0.

        #move to first wall boundary
        guess_y = (y_grid[1] - y_grid[0])/2.
        guess_z = (z_grid[1] - z_grid[0])/2.
        y_func = lambda s,j: self.y(s) - y_grid[j]
        z_func = lambda s,k: self.z(s) - z_grid[k]
        j_track = []
        k_track = []

        while(s_total < self.s_finish and self.j >= 0 and self.k >= 0):
            j_track.append(self.j)
            k_track.append(self.k)
            s_max_y = fsolve(y_func,s_total + guess_y,(self.j))[0]
            s_max_z = fsolve(z_func,s_total + guess_z,(self.k))[0]
            if s_max_y < s_max_z:
                s_total = s_max_y
                self.j = self.j - 1
            else:
                s_total = s_max_z
                self.k = self.k - 1

            print(s_max_y,s_max_z)
            print("j:%i; k:%i" % (self.j,self.k))
            print("Complete: %.2f" % (s_total/self.s_finish))
            print()
        js = np.array(j_track)
        ks = np.array(k_track)
        plt.plot(y_grid[js]/const.AU,z_grid[ks]/const.AU,"o")
        plt.show()

    def plot_los(self):
        ss = np.linspace(0,self.s_finish)
        ys = self.y(ss)
        zs = self.z(ss)
        plt.plot(ys/const.AU,zs/const.AU)
        plt.show()

    def x(self,s):
        return self.x0 * np.ones_like(s)

    def y(self,s):
        return self.y0 - s * np.sin(self.orn.theta)

    def z(self,s):
        return self.z0 - s * np.cos(self.orn.theta)

    def cartesian(self,s):
        x = self.x0 * np.ones_like(s)
        y = self.y0 - s * np.sin(self.orn.theta)
        z = self.z0 - s * np.cos(self.orn.theta)
        return np.array([x,y,z])

    def polar(self,s):
        x = self.x0 * np.ones_like(s)
        y = self.y0 - s * np.sin(self.orn.theta)
        z = self.z0 - s * np.cos(self.orn.theta)
        r = np.sqrt(x**2. + y**2.)
        phi = np.arctan2(y,x)
        return np.array([r,z,phi])

    def coords(self,s):
        '''Function to return Cartesian and Cyclindrical Corods'''
        x = self.x0 * np.ones_like(s)
        y = self.y0 - s * np.sin(self.orn.theta)
        z = self.z0 - s * np.cos(self.orn.theta)
        r = np.sqrt(x**2. + y**2.)
        phi = np.arctan2(y,x)
        return (np.array([x,y,z]),np.array([r,z,phi]))

    def tau_dis(self,s):
        ind = np.where(self.ss <= s)[0]
        return np.trapz(self.K_nus[ind],self.ss[ind])

    def plot_path(self,ax,fmt="bo"):
        ss = np.linspace(0,self.s_finish,num=10)
        pos_array = self.cartesian(ss)
        ys = pos_array[1,:]
        zs = pos_array[2,:]
        #print(ys)
        #print(zs)
        ax.plot(ys,zs,"bo")
        return

    def graph_on_ax(self,ax,fmt="b"):
        '''ax is a matplotlib axes object.'''
        ys = np.linspace(-6.0,6.0)
        ax.plot(self.y0,self.z0,"bo")
        ax.plot(ys,self.orn.line_through_point((self.y0,self.z0))(ys),fmt)
        return

    def check_total_path_length(self):
        distance = np.sqrt((self.yf - self.y0)**2 + (self.zf - self.z0)**2)
        print("Distance={distance:.2f} ; s_finish={s_finish:.2f}".format(distance=distance,s_finish=self.s_finish))

    def v_k(self,pol):
        return self.disk.v_phi(pol[0]) * np.cos(pol[2]) * np.sin(self.orn.theta)

    def K_nu(self,pol,T,rho):
        Z = np.sqrt(1.0 + ((2. * T)/self.model.T1)**2) #Partition function
        n_l = self.x0 * rho/const.m0 * self.model.gl * np.exp(-self.model.El/(const.k * T)) / Z
        K_dust = rho * self.k_nu
        K_CO = n_l * self.sigma_nu(pol,T)
        return K_dust + K_CO

    def K_nu_tau(self,s):
        pol = self.polar(s)
        T = self.disk.T(pol[0])
        rho = self.disk.rho(pol[0],pol[1])
        Z = np.sqrt(1.0 + ((2. * T)/self.model.T1)**2) #Partition function
        n_l = self.x0 * rho/const.m0 * self.model.gl * np.exp(-self.model.El/(const.k * T)) / Z
        K_dust = rho * self.k_nu
        K_CO = n_l * self.sigma_nu(pol,T)
        return K_dust + K_CO

    def sigma_nu(self,pol,T):
        return self.model.sigma0 * self.phi_nu(pol,T) * (1. - np.exp(-const.h * self.nu/(const.k * T)))

    def phi_nu(self,pol,T):
        Delta_V = np.sqrt(2. * const.k * T/const.m_CO +  self.model.vturb**2)
        v_obs = const.c/self.model.nu0 * (self.nu - self.model.nu0)
        Delta_v = self.v_k(pol) - v_obs 
        return const.c/(self.model.nu0 * np.sqrt(np.pi) * Delta_V) * np.exp(-1. * Delta_v**2/Delta_V**2)

    def tau_general(self,s):
        return quad(self.K_nu_tau,0.0,s)[0]

    def S(self,s):
        '''Used only for plotting S along line of sight'''
        car,(r,z,phi) = self.coords(s)
        T = self.disk.T(r)
        S = rad.S(T,self.nu)
        return S

    def dI(self,s):
        pol = self.polar(s)
        T = self.disk.T(pol[0])
        S = rad.S(T,self.nu)
        rho = self.disk.rho(pol[0],pol[1])
        K_nu = self.K_nu(pol,T,rho)
        tau = self.tau(s)
        return S * np.exp(-tau) * K_nu

    def plot_K_nu(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ss = np.linspace(0,self.s_finish,num=1000)
        ax.semilogy(ss/const.AU,self.K_nu(ss))
        ax.set_xlabel(r"$s$ (AU)")
        ax.set_ylabel(r"$K_\nu(s)\quad[{\rm cm}^{-1}]$")
        fig.subplots_adjust(left=0.20)
        fig.savefig("Plots/Test_Rad/k_nu_vs_s.png")

    def plot_tau(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ss = np.linspace(0,self.s_finish,num=100)
        #func = lambda s: self.tau(s)
        #func = np.vectorize(func)
        taus = self.tau(ss)
        ax.set_xlabel(r"$s$ (AU)")
        ax.set_ylabel(r"$\tau(s)$")
        ax.semilogy(ss/const.AU, taus)
        fig.savefig("Plots/Test_Rad/tau_vs_s.png")

    def plot_S(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ss = np.linspace(0,self.s_finish,num=100)
        ax.set_xlabel(r"$s$ (AU)")
        ax.set_ylabel(r"$S_\nu(s)\quad \left[ \frac{{\rm erg}}{{\rm cm}^2\cdot {\rm s\cdot ster\cdot Hz}} \right]$")
        ax.semilogy(ss/const.AU, self.S(ss)) #requires rad.S
        fig.subplots_adjust(left=0.20)
        fig.savefig("Plots/Test_Rad/S_vs_s.png")

    def plot_dI(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ss = np.linspace(0,self.s_finish,num=100)
        #dIs = np.array([self.dI(i) for i in ss])
        dIs = self.dI(ss)
        ax.set_xlabel(r"$s$ (AU)")
        ax.set_ylabel(r"$dI(s)\quad \left[ \frac{{\rm erg}}{{\rm cm}^3\cdot {\rm s\cdot ster\cdot Hz}} \right] $")
        ax.semilogy(ss/const.AU,dIs)
        fig.subplots_adjust(left=0.20)
        fig.savefig("Plots/Test_Rad/dI_vs_s.png")

    def __str__(self):
        return """
        alpha = {alpha:.2f}
        delta = {delta:.2f}
        x0 = {x0:.2f}; y0 = {y0:.2f}; z0 = {z0:.2f}
        xf = {xf:.2f}; yf = {yf:.2f}; zf = {zf:.2f}""".format(alpha=self.alpha,delta=self.delta,x0=self.x0,y0=self.y0,z0=self.z0,xf=self.xf,yf=self.yf,zf=self.zf)

class Grid:
    def __init__(self):
        self.N_xy = 10
        self.N_z = 5

        self.x_grid = np.linspace(-l,l,num=self.N_xy)
        self.y_grid = self.x_grid.copy()
        self.z_grid = np.linspace(-w,w,num=self.N_z)

        self.x_cells = (self.x_grid[0:-1] + self.x_grid[1:])/2.
        self.y_cells = self.x_cells.copy()
        self.z_cells = (self.z_grid[0:-1] + self.z_grid[1:])/2.

        self.cells = np.empty((self.N_xy - 1,self.N_xy - 1,self.N_z - 1,3))
        for i in range(self.N_xy - 1):
            for j in range(self.N_xy-1):
                for k in range(self.N_z-1):
                    self.cells[i,j,k] = np.array([self.x_cells[i],self.y_cells[j],self.z_cells[k]])

