#!/usr/bin/env python

import constants as const
import numpy as np

def k_nu(nu,beta):
    '''Returns the dust opacity coefficient, given by Eqn 3, Isella et al. 2007'''
    return 0.01 * (const.c/( nu * 0.13))**(-1. * beta) #cm^2 g^-1

class Radiation:
    def __init__(self,model,delta_v):
        self.model = model
        self.nu = self.model.center_frequency(delta_v)
        self.k_nu = k_nu(self.nu,self.model.beta)

    def S(self):
        rs = self.model.grid.cylin_grid[:,0]
        ind_no = np.where(rs < self.model.disk.r_in)[0]
        ind_yes = np.where(rs >= self.model.disk.r_in)[0]
        r = rs[ind_yes]
        T = self.model.disk.T(r) 
        S = 2. * const.h * self.nu**3 / const.c**2 / (np.exp(const.h * self.nu / (const.k * T)) - 1.0)
        S_full = np.zeros_like(rs)
        S_full[ind_yes] = S
        self.S_cart = S_full.reshape(self.model.grid.full_shape_one)

    def K(self):
        rs = self.model.grid.cylin_grid[:,0]
        ind_no = np.where(rs < self.model.disk.r_in)[0]
        ind_yes = np.where(rs >= self.model.disk.r_in)[0]
        zs = self.model.grid.cylin_grid[:,1]
        phis = self.model.grid.cylin_grid[:,2]
        r = rs[ind_yes]
        z = zs[ind_yes]
        phi = phis[ind_yes]
        theta = self.model.orientation.theta
        rho = self.model.disk.rho(r,z)
        T = self.model.disk.T(r) 
        Z = np.sqrt(1.0 + ((2. * T)/self.model.T1)**2) 
        n_l = self.model.x0 * rho/const.m0 * self.model.gl * np.exp(-self.model.El/(const.k * T)) / Z
        K_dust = rho * self.k_nu
        Delta_V = np.sqrt(2. * const.k * T/const.m_CO +  self.model.vturb**2)
        v_obs = const.c/self.model.nu0 * (self.nu - self.model.nu0)
        v_k = self.model.disk.v_phi(r) * np.cos(phi) * np.sin(theta)
        Delta_v = v_k - v_obs 
        phi_nu = const.c/(self.model.nu0 * np.sqrt(np.pi) * Delta_V) * np.exp(-1. * Delta_v**2/Delta_V**2)
        sigma_nu = self.model.sigma0 * phi_nu * (1. - np.exp(-const.h * self.nu/(const.k * T)))
        K_CO = n_l * sigma_nu
        K = K_dust + K_CO
        K_full = np.zeros_like(rs)
        K_full[ind_yes] = K
        self.K_cart = K_full.reshape(self.model.grid.full_shape_one)

