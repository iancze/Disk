#!/usr/bin/env python
#Super class to build the model

import disk
import constants as const

#parameters
M_star = 2.6 * const.M_sun
m0 = 2.37 * const.amu #mean molecular weight of gas
r0 = 550. * const.AU
Sigma0 = 90. #g/cm^2
p = 0.6
T_r0 = 40. #K
q = 0.5
beta = 1.0
#M_star,r0,T_r0,q,Sigma0,p,m0,X0):
disk1 = disk.Disk(M_star,r0, T_r0, q, Sigma0, p, m0, 1.0)
disk1.plot_rho()
