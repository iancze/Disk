#!/usr/bin/env python

import constants as const
import numpy as np

class Channel:
    def __init__(self,nu,beta):
        pass

def k_nu(nu,beta):
    '''Returns the dust opacity coefficient, given by Eqn 3, Isella et al. 2007'''
    return 0.01 * (const.c/( nu * 0.13))**(-1. * beta) #cm^2 g^-1

def k_nu_dust(rho, nu, beta):
    return rho * k_nu(nu,beta)

def n_l(s):
    '''Boltzmann Equation'''
    return 

def S(T,nu):
    return 2. * const.h * nu**3 / const.c**2 / (np.exp(const.h * nu / (const.k * T)) - 1.0)

