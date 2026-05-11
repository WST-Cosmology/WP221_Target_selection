import numpy as np
import scipy as sp
#import pyccl as ccl
import math
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from scipy import integrate
from functools import partial
from scipy.integrate import quad, dblquad

import pyccl as ccl
import camb
from camb import model, initialpower

deltac=1.686
c_ls=300*10**3

def Vsurvey(z,dz,sky, cosmo):    
    h = cosmo['h']
    Omega     = sky*(math.pi/180)**2 # get rid of unit
    d2        = ccl.comoving_radial_distance(cosmo,1/(1+z))
    d3        = ccl.comoving_radial_distance(cosmo,1/(1+z+dz))
    return Omega/3 * (d3**3 - d2**3)*h**3 #h**3 to have Mpc/h

def D(z, cosmo):
    return ccl.growth_factor(cosmo,1/(1+z))*0.708

def f(z, cosmo):
    return ccl.growth_rate(cosmo,1/(1+z))

def Pm(k,z1, cosmo):
    return ccl.linear_matter_power(cosmo, k*h, 1/(1+z1))*h**3

def Pz0(k, cosmo):
     return ccl.linear_matter_power(cosmo, k*h, 1/(1+0))*h**3
 
def sigma_8(z1, cosmo):
     return (Pm(0.01,z1)/Pz0(0.01))**0.5*ccl.sigma8(cosmo)
 
def sigma_chi(z, cosmo):
    sigma_z=0.001
    h = cosmo['h']
    return c*(1+z)/(ccl.h_over_h0(cosmo,1/(1+z))*h*100)*sigma_z/(1+z)
        
def Dl(z, cosmo):
    return ccl.background.luminosity_distance(cosmo, 1/(1+z))

def T(k,mod, cosmo):
    '''
    return the transfer function T(k), either with camb, or bbks approx. 
    '''
    Omega_m=cosmo['Omega_m']
    Omega_b=cosmo['Omega_b']
    Omega_c=cosmo['Omega_c']
    Omega_l=1-Omega_m
    h = cosmo['h']
    ns = cosmo['n_s']
    if mod=='camb':
        
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=100*h, ombh2=Omega_b*h**2, omch2=Omega_c*h**2)
        pars.InitPower.set_params(ns=ns)
        pars.set_matter_power(redshifts=[0], kmax=10)
        results= camb.get_results(pars)
        trans = results.get_matter_transfer_data()
        
        #get kh - the values of k/h at which they are calculated
        kh = trans.transfer_data[0,:,0]
        #transfer functions for different variables, e.g. CDM density and the Weyl potential
        #CDM perturbations have grown, Weyl is O(1) of primordial value on large scales
        delta = trans.transfer_data[model.Transfer_cdm-1,:,0]
        Tk_disc=trans.transfer_data[model.Transfer_tot-1,:,0]
        Tk_cont=sp.interpolate.interp1d(kh,Tk_disc/Tk_disc[0])
        return Tk_cont(k)
    elif mod=='bbks':
        #bbks approx
        k=k*h# We want k in Mpc-1 instead of Mpc-1h so 
        q=k/(Omega_m*h**2*math.exp(-Omega_b*(1+(2*h)**0.5/Omega_m)))
        return math.log(1+2.34*q)/2.34/(q)*(1+3.89*q+(16.2*q)**2+(5.47*q)**3+(6.71*q)**4)**(-0.25)#ccl.get_transfer(math.log(k),a)
    else:
        print('mod of the transfer function unknown')

def Pm(k,z1, cosmo):
    h = cosmo['h']
    return ccl.linear_matter_power(cosmo, k*h, 1/(1+z1))*h**3 # we use a hMpc**-1 convention for PS

def deltab_test(k,fNL,bg,p,mod,z, cosmo):
        Omega_m=cosmo['Omega_m']
        h = cosmo['h']
        return 3*fNL*(bg-p)*deltac*Omega_m/(k**2*T(k,mod,cosmo)*D(z,cosmo))*(100*h/c_ls)**2