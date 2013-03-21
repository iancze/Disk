#!/usr/bin/env python

import numpy as np
import constants as const
import radiation as rad
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fsolve

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
    def __init__(self,model,orientation,alpha,delta):
        '''alpha_phys and delta_phys are the actual projected distances at the location of the disk.'''
        self.model = model #Link to the model
        self.orn = orientation #Link to the orientation object
        self.alpha = alpha
        self.alpha_phys = self.alpha * self.orn.distance * (const.AU/const.pc)
        self.delta = delta
        self.delta_phys = self.delta * self.orn.distance * (const.AU/const.pc)
        self.x0 = self.alpha_phys
        self.xf = self.x0
        self.calc_y0_z0()
        self.disk = self.orn.model.disk
        #self.ss = np.arange(0.,self.s_finish,0.1*const.AU)
        self.grid = self.model.grid

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

    def integrate(self):
        x_grid = self.grid.x_grid
        y_grid = self.grid.y_grid
        z_grid = self.grid.z_grid

        #Define initial starting cell
        self.i = np.max(np.where(self.x0 > x_grid)[0])
        self.j = np.max(np.where(self.y0 > y_grid)[0])
        self.k = np.max(np.where(self.z0 > z_grid)[0])

        s_total = 0.
        tau_total = 0.
        I_total = 0.

        js = []
        ks = []

        guess_y = (y_grid[1] - y_grid[0])/2.
        guess_z = (z_grid[1] - z_grid[0])/2.
        y_func = lambda s,j: self.y(s) - y_grid[j]
        z_func = lambda s,k: self.z(s) - z_grid[k]

        while(s_total < self.s_finish and self.j >= 0 and self.k >= 0):
            K = self.grid.K[self.i,self.j,self.k]
            S = self.grid.S[self.i,self.j,self.k]
            js.append(self.j)
            ks.append(self.k)
            s_max_y,info_dict,flag_y,msg = fsolve(y_func,s_total + guess_y,(self.j),full_output=True)
            s_max_z,info_dict,flag_z,msg = fsolve(z_func,s_total + guess_z,(self.k),full_output=True)
            s_max_y = s_max_y[0]
            s_max_z = s_max_z[0]
#            if flag_y != 1:
#                print("Not converged y",self.j,self.k)
#            if flag_z != 1:
#                print("Not converged z",self.j,self.k)
            if s_max_y < s_max_z:
                ds = s_max_y - s_total
                s_total = s_max_y
                self.j = self.j - 1
            else:
                ds = s_max_z - s_total
                s_total = s_max_z
                self.k = self.k - 1

            delta_tau = K * ds
            delta_I = K * S * np.exp(-tau_total) * ds
            I_total = I_total + delta_I
            tau_total = tau_total + delta_tau

        def plot_path():
            js = np.array(js)
            ks = np.array(ks)
                
            y_cells = self.grid.y_cells/const.AU
            z_cells = self.grid.z_cells/const.AU
            y_grid = y_grid/const.AU
            z_grid = z_grid/const.AU
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ss = np.linspace(0,self.s_finish)
            ys = self.y(ss)/const.AU
            zs = self.z(ss)/const.AU
            ax.plot(ys,zs)
            for i in y_grid:
                ax.axvline(i,color="k",ls=":")
            for i in z_grid:
                ax.axhline(i,color="k",ls=":")
            for i in range(len(z_cells)):
                plt.plot(y_cells,z_cells[i]*np.ones_like(y_cells),"ko")
            plt.plot(y_grid[js], z_grid[ks],"o")
            width = y_cells[1] - y_cells[0]
            height = z_cells[1] - z_cells[0]
            for i in range(len(y_grid[js])):
                ax.add_patch(plt.Rectangle((y_grid[js][i],z_grid[ks][i]),width=width,height=height,alpha=0.2))
            ax.set_xlim(-500,500)
            ax.set_xlabel("y")
            ax.set_ylabel("z")
            plt.show()

        print("Complete: %.2f" % (s_total/self.s_finish))
        return I_total

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

    def __str__(self):
        return """
        alpha = {alpha:.2f}
        delta = {delta:.2f}
        x0 = {x0:.2f}; y0 = {y0:.2f}; z0 = {z0:.2f}
        xf = {xf:.2f}; yf = {yf:.2f}; zf = {zf:.2f}""".format(alpha=self.alpha,delta=self.delta,x0=self.x0,y0=self.y0,z0=self.z0,xf=self.xf,yf=self.yf,zf=self.zf)


class Grid:
    def __init__(self,model):
        self.model = model
        self.N_xy = 500
        self.N_z = 100
        self.flat_shape = ((self.N_xy-1)**2 * (self.N_z -1),3)
        self.full_shape = ((self.N_xy-1),(self.N_xy -1), (self.N_z - 1),3)
        self.full_shape_one = ((self.N_xy-1),(self.N_xy -1), (self.N_z - 1),1)


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

        self.flat_grid = self.cells.reshape(self.flat_shape)
        self.polar()

    def polar(self):
        self.cylin_grid = self.flat_grid.copy()
        x = self.flat_grid[:,0]
        y = self.flat_grid[:,1]
        z = self.flat_grid[:,2]
        self.cylin_grid[:,0] = np.sqrt(x**2 + y**2) #r
        self.cylin_grid[:,1] = z                    #z
        self.cylin_grid[:,2] = np.arctan2(y,x)      #phi

    def plot_flat_grid(self):
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        xs = self.flat_grid[:,0]/const.AU
        ys = self.flat_grid[:,1]/const.AU
        zs = self.flat_grid[:,2]/const.AU
        ax.scatter(xs,ys,zs)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.show()

    def calc_S(self):
        self.model.radiation.S()
        self.S = self.model.radiation.S_cart

    def calc_K(self):
        self.model.radiation.K()
        self.K = self.model.radiation.K_cart
