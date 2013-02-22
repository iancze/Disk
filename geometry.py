#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#Testing whether the geometry is correct in determining starting positions on the box.

#Box is total width (in z direction) of 2w and total length (in y direction) of 2l
global w,l
w = 2.
l = 4.

class Orientation:
    def __init__(self,theta):
        self.theta = theta
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
    def __init__(self,orientation,delta,alpha):
        self.orn = orientation #Link to the orientation object
        self.alpha = alpha
        self.delta = delta
        self.x0 = self.alpha
        self.xf = self.x0
        self.calc_y0_z0()
        pass

    def calc_y0_z0(self):
        '''Given the inital geometry via the Orientation object, calculate y0, z0, yf, zf, and s_finish.'''
        if self.delta > self.orn.delta_c:
            self.z0 = w
            self.y0 = w * np.tan(self.orn.theta) - self.delta/np.cos(self.orn.theta)
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
            self.z0 = self.delta/np.sin(self.orn.theta) + l/np.tan(self.orn.theta)
            k2 = (self.z0 + w) * np.tan(self.orn.theta)
            if k2 > 2. * l:
                self.yf = -l
                self.zf = self.z0 - 2. * l/np.tan(self.orn.theta)
            else:
                self.yf = self.y0 - k2
                self.zf = -w
        self.s_finish = (self.y0 - self.yf)/np.sin(self.orn.theta)

    def cartesian(self,s):
        x = self.x0 * np.ones_like(s)
        y = self.y0 - s * np.sin(self.orn.theta)
        z = self.z0 - s * np.cos(self.orn.theta)
        return np.array([x,y,z])

    def polar(self,s):
        x,y,z = self.cartesian(s)
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

    def dI(self,s):
        car,cyl = self.coords(s)
        #rho = 
        #T = 
        #v_phi = 

    def __str__(self):
        return """
        alpha = {alpha:.2f}
        delta = {delta:.2f}
        x0 = {x0:.2f}; y0 = {y0:.2f}; z0 = {z0:.2f}
        xf = {xf:.2f}; yf = {yf:.2f}; zf = {zf:.2f}""".format(alpha=self.alpha,delta=self.delta,x0=self.x0,y0=self.y0,z0=self.z0,xf=self.xf,yf=self.yf,zf=self.zf)

