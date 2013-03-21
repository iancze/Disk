#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

alpha_grid = np.load("alpha.npy")
delta_grid = np.load("delta.npy")
#img = [np.transpose(np.load("img_45__2.npy")),np.transpose(np.load("img_45_0.npy")),np.transpose(np.load("img_45_2.npy")),np.transpose(np.load("img_45_4.npy"))]
img = [np.transpose(np.load("img_90__2.npy")),np.transpose(np.load("img_90_0.npy")),np.transpose(np.load("img_90_2.npy")),np.transpose(np.load("img_90_4.npy"))]

fig = plt.figure(figsize=(6,6))
grid = AxesGrid(fig,111,nrows_ncols = (2,2), axes_pad = 0.2, label_mode="1",share_all=True,cbar_location="top",cbar_mode="each",cbar_size="7%",cbar_pad="2%",)
for i in range(4):
    im = grid[i].imshow(np.log10(img[i]),interpolation='gaussian',extent=[-0.5,0.5,-0.5,0.5])
    grid.cbar_axes[i].colorbar(im)

#CS0 = grid[0].contourf(alpha_grid,delta_grid,np.log10(img[0]))
#CS1 = grid[1].contourf(alpha_grid,delta_grid,np.log10(img[1]))
#CS2 = grid[2].contourf(alpha_grid,delta_grid,np.log10(img[2]))
#CS3 = grid[3].contourf(alpha_grid,delta_grid,np.log10(img[3]))

#print(img)
#CS = ax.contourf(alpha_grid,delta_grid,np.log10(img))#,levels=levels)
#CS = ax.contourf(img)#,levels=levels)
#ax.set_xlabel(r'$\Delta\alpha$ ["]')
#ax.set_ylabel(r'$\Delta\delta$ ["]')
fig.text(0.5,0.975,r"${}^{13}{\rm CO}\;J = 1 \rightarrow 0\;\theta=90^\circ\;\log_{10} \left[ \frac{{\rm erg}}{{\rm cm}^2\cdot {\rm s\cdot ster\cdot Hz}} \right]$",horizontalalignment="center",verticalalignment="top") 
#cbar = plt.colorbar(CS)

plt.show()


