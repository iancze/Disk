#!/usr/bin/env python

import numpy as np

#Testing whether the geometry is correct in determining starting positions on the box.

#Box is total width (in z direction) of 2w and total length (in y direction) of 2l
global w,l
w = 2.
l = 4.

class StartingCoords:
    def __init__(self,theta,delta,alpha):
        self.theta = theta
        self.alpha = alpha
        self.delta = delta
        self.x0 = self.alpha
        self.delta_limit = w * np.sin(self.theta) + l * np.cos(self.theta)
        self.theta_c = np.arctan(l/w)
        self.delta_c = w * np.sin(self.theta) - l * np.cos(self.theta)

    def __str__(self):
        return """
        theta = {theta:.2f}
        theta_c = {theta_c:.2f}
        alpha = {alpha:.2f}
        delta = {delta:.2f}
        delta_c = {delta_c:.2f}
        delta_high = {delta_high:.2f}
        delta_low = {delta_low:.2f}
        x0 = {x0:.2f}""".format(theta=self.theta,theta_c=self.theta_c,alpha=self.alpha,delta=self.delta,delta_c=self.delta_c,delta_high=self.delta_limit,delta_low=(-1.*self.delta_limit),x0=self.x0)

