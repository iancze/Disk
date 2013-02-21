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
        self.calc_y0_z0()
        pass

    def calc_y0_z0(self):
        if self.delta > self.orn.delta_c:
            self.z0 = w
            if self.delta > 0:
                self.y0 = w * np.tan(self.orn.theta) - self.delta/np.cos(self.orn.theta)
            else:
                self.y0 = w * np.tan(self.orn.theta) + self.delta/np.cos(self.orn.theta)
        else:
            self.y0 = l
            self.z0 = l/np.tan(self.orn.theta) - self.delta/np.sin(self.orn.theta)

    def graph_on_ax(self,ax):
        '''ax is a matplotlib axes object.'''
        #ax.plot([1,2],[1,2])
        ax.grid(True)
        ax.axhline(0,color="k")
        ax.axvline(0,color="k")
        ax.add_patch(plt.Rectangle((-l,-w),2*l,2*w,fill=False))
        ax.set_xlim(-1.5*l,1.5*l)
        ax.set_ylim(-1.5*w,1.5*w)
        ax.set_xlabel("y")
        ax.set_ylabel("z")
        ys = np.linspace(-6.0,6.0)
        ax.plot(ys,self.orn.center_line()(ys))
        #ax.plot(ys,self.line_through_point(
        ax.plot(self.y0,self.z0,"bo")
        return

    def __str__(self):
        return """
        alpha = {alpha:.2f}
        delta = {delta:.2f}
        x0 = {x0:.2f}; y0 = {y0:.2f}; z0 = {z0:.2f}""".format(alpha=self.alpha,delta=self.delta,x0=self.x0,y0=self.y0,z0=self.z0)
