#!/usr/bin/env python
#Super class to build the model

import disk
import radiation
import geometry
import constants as const
import numpy as np
import cProfile
import pstats
import matplotlib.pyplot as plt

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
orn_params = {"theta":90. * np.pi/180.,
        "distance":150. * const.pc}
img_width = 1. #arcseconds
img_height = 1. #arcseconds
res_elements = 75.

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
    def __init__(self,l,species,disk_params,orn_params,delta_v):
        '''disk_params is a dictionary of disk parameters'''
        self.distance = 150. * const.pc #distance to the disk from earth
        self.l = l
        self.beta = beta
        self.vturb = 0.0
        self.species = species
        self.set_species_constants()
        self.disk = disk.Disk(self,**disk_params)
        self.orientation = geometry.Orientation(self,**orn_params)
        self.nu = self.center_frequency(delta_v)
        self.k_nu = radiation.k_nu(self.nu,self.beta)
        print("Generating Grid")
        self.grid = geometry.Grid(self)

    def generate_images(self,fname):
        '''Create an appropriately-spaced array of line of sights'''
        #Create a grid of alpha, delta, based upon a resolution parameter and a range
        alphas = np.linspace(-img_width/2.,img_width/2,res_elements)
        deltas = np.linspace(-img_height/2.,img_height/2,res_elements)
        intensity_array = np.zeros((res_elements,res_elements))
        #intensity_array = np.zeros_like(alphas)
        AA,DD = np.meshgrid(alphas,deltas)
        np.save("alpha.npy",AA)
        np.save("delta.npy",DD)
        print("Saved coord grid")
        la = len(alphas)
        for i in range(la):
            for j in range(len(deltas)):
                los = geometry.LineOfSight(self,self.orientation,alphas[i],deltas[j])
                intensity_array[i,j] = los.integrate()
            print("%.2f" % ((1.0*i)/la))
        #np.save("img{theta:.0f}.npy".format(theta=orn_params["theta"]*180./np.pi),intensity_array)
        np.save(fname,intensity_array)

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

def test_grid_walker():
    mod1 = Model(0,"13CO",disk_params,orn_params,rad_params)
    mod1.grid.calc_S()
    mod1.grid.calc_K()
    los = geometry.LineOfSight(mod1,mod1.orientation,0.0,1.0)
    los.integrate()

def test_grid():
    mod1 = Model(0,"13CO",disk_params,orn_params,rad_params)
    mod1.grid.plot_flat_grid()

def main():
    #test_grid_walker()
    #test_grid()
    #names = ["img_45__2.npy","img_45_0.npy","img_45_2.npy","img_45_4.npy"]
    names = ["img_90__2.npy","img_90_0.npy","img_90_2.npy","img_90_4.npy"]
    rad_list = [-2*const.kms, 0.0,2*const.kms,4*const.kms]
    for i in range(len(names)):
        print(names[i])
        mod1 = Model(0,"13CO",disk_params,orn_params,delta_v=rad_list[i])
        print("Calculating grid")
        mod1.grid.calc_S_and_K()
        print("Starting Integration")
        mod1.generate_images(names[i])
    #nu = mod1.center_frequency(0.0)
    #los = geometry.LineOfSight(mod1,mod1.orientation,0.0,0.0)
    #los.plot_spher_vs_s()
    #los.plot_los()
    #los.walk_along_grid()
    #print(los.integrate())
    #nu_off = mod1.center_frequency(0.0)
   
if __name__=="__main__":
    main()
