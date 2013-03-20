#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure()
ax = fig.add_subplot(111)
alpha_grid = np.load("alpha.npy")
delta_grid = np.load("delta.npy")
img = np.transpose(np.load("img45.npy"))
#print(img)
levels = np.linspace(1e-11,2e-11,num=20)
CS = ax.contourf(img)#,levels=levels)
ax.set_xlabel(r'$\Delta\alpha$ ["]')
ax.set_ylabel(r'$\Delta\delta$ ["]')
ax.set_title(r"${}^{13}{\rm CO}\;J = 1 \rightarrow 0\;\theta=45^\circ\;\left[ \frac{{\rm erg}}{{\rm cm}^2\cdot {\rm s\cdot ster\cdot Hz}} \right]$") 
cbar = plt.colorbar(CS)

plt.show()


