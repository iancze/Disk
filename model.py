#!/usr/bin/env python
#Super class to build the model

import disk
import radiation
import geometry
import constants as const
import numpy as np
import cProfile
import pstats


#parameters
beta = 1.0
l = 1

x0_dict = {"12CO": 7.0e-5, "13CO":1.0e-6, "C18O":1.3e-7}
dipole_dict = {"12CO":1.22e-19,"13CO":1.1e-19} #e.s.u. cm
T1_dict = {"12CO":5.53, "13CO":5.29, "C18O":5.27} #Kelvin
nu_dict = {"12CO":115.27e9 ,"13CO":110.2e9, "C18O":109.78e9}#Hz
#M_star,r0,T_r0,q,Sigma0,p,m0,X0):
global disk_params
disk_params = {"M_star":2.6 * const.M_sun,
"r0":550. * const.AU, 
"T_r0":40.,#K 
"q": 0.5, 
"Sigma0": 90., #g/cm^2
"p":0.6}
orn_params = {"theta":np.pi/4.,
        "distance":150. * const.pc}
img_width = 2. #arcseconds
img_height = 2. #arcseconds

class Model:
    '''General Container class to generate a model for a given set of parameters. The parameters are:
        * For the Disk
            * theta: inclination
            * beta: dust opacity index
            * r0: scale radius
            * Sigma0: surface density
            * p: density exponent
            * T_r0: temperature
            * q: temperature exponent
            * v_turb: turbulent velocity (generally = 0, for now)
        * For a species/transition
            * l
            * species ('12CO','13CO') '''
    def __init__(self,l,species,disk_params,orn_params):
        '''disk_params is a dictionary of disk parameters'''
        self.distance = 150. * const.pc #distance to the disk from earth
        self.l = l
        self.species = species
        self.set_species_constants()
        self.disk = disk.Disk(self,**disk_params)
        self.orientation = geometry.Orientation(self,**orn_params)
        self.beta = beta
        self.vturb = 0.0

    def generate_images(self,nu):
        '''Create an appropriately-spaced array of line of sights'''
        #Create a grid of alpha, delta, based upon a resolution parameter and a range
        alphas = np.linspace(-img_width/2.,img_width/2,20)
        deltas = np.linspace(-img_height/2.,img_height/2,20)
        alpha_grid, delta_grid = np.meshgrid(alphas, deltas)
        np.save("alpha.npy", alpha_grid)
        np.save("delta.npy", delta_grid)
        #for each (alpha,delta) pair, create a LoS, integrate, and return a value
        #do this for each nu within a range which will be the channel map
        #Right now do this for (alpha,delta) = (0.0,0.0)
        intensity_array = np.zeros_like(alpha_grid)
        for i in range(len(alpha_grid)):
            for j in range(len(delta_grid)):
                los = geometry.LineOfSight(self,self.orientation,alpha_grid[i][j],delta_grid[i][j],nu)
                intensity_array[i][j] = los.integrate()

        intensity_array = np.array(intensity_array)
        np.save("img45_2kms.npy",intensity_array)
        #LoS.plot_K_nu()
        #LoS.plot_tau()
        #LoS.plot_S()
        #LoS.plot_dI()
        #print(LoS1.integrate())
        #print(LoS2.integrate())

    def set_species_constants(self):
        self.x0 = x0_dict[self.species]
        self.gl = 2 * self.l + 1
        self.dipole = dipole_dict[self.species]
        self.T1 = T1_dict[self.species]
        self.nu0 = nu_dict[self.species]
        self.El = 0.5 * (self.l * (self.l + 1)) * const.k * self.T1
        self.sigma0 = 8. * np.pi**3 * const.k * self.T1 / (const.h**2 * const.c) * (self.l + 1)**2/(2. * self.l + 1) * self.dipole**2

    def center_frequency(self,velocity):
        '''Takes in the Delta v for the channel (in cm/s) and returns the corresponding frequency (in Hz)'''
        return self.nu0 * velocity / const.c + self.nu0

def main():
    mod1 = Model(0,"13CO",disk_params,orn_params)
    nu_off = mod1.center_frequency(2.*const.kms)
    mod1.generate_images(nu_off)
   
if __name__=="__main__":
    main()
