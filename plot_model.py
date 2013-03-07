#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure()
ax = fig.add_subplot(111)
alpha_grid = np.load("alpha.npy")
delta_grid = np.load("delta.npy")
img = np.load("img45_2kms.npy")
levels = np.linspace(5e-13,1.3e-12,num=20)
CS = ax.contourf(alpha_grid,delta_grid,img,levels=levels)
ax.set_xlabel(r'$\Delta\alpha$ ["]')
ax.set_ylabel(r'$\Delta\delta$ ["]')
ax.set_title(r"${}^{13}{\rm CO}\;J = 1 \rightarrow 0\;\theta=45^\circ\;\left[ \frac{{\rm erg}}{{\rm cm}^2\cdot {\rm s\cdot ster\cdot Hz}} \right]$") 
cbar = plt.colorbar(CS)

plt.show()
