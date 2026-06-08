
import pyccl as ccl
import math
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from scipy.integrate import quad, dblquad
from scipy import interpolate

#Set up a cosmology
yourz=3
h=0.67
H0=100*h
Omega_m=0.315
c_ls=299.792*10**3
cosmo = ccl.Cosmology(Omega_c=0.27, Omega_b=Omega_m-0.27, h=h, A_s=2.1e-9, n_s=0.96)
a=1./(1+yourz)

Ell=range(1,4000)
nlim=10000

def chi(z):
    return(ccl.comoving_radial_distance(cosmo, 1/(1+z)))

def chi_3d(r,z):
    return ccl.correlations.correlation_3d(cosmo,1/(1+z),r,'delta_matter:delta_matter')

def PNL(l,z):
    #return ccl.power.nonlin_matter_power(cosmo, k=(l+0.5)/chi(z), a=1/(1+z))
    return ccl.power.nonlin_power(cosmo, k=(l+0.5)/chi(z), a=1/(1+z), p_of_k_a='delta_matter:delta_matter')

def Plin(l,z):
    #return ccl.power.linear_matter_power(cosmo, k=(l+0.5)/chi(z), a=1/(1+z))
    return ccl.power.linear_power(cosmo, k=(l+0.5)/chi(z), a=1/(1+z), p_of_k_a='delta_matter:delta_matter')

def xi_dm_theta(theta, z , typ):
    '''same as before for theta in deg instead of rp'''
    Hz=100*h*ccl.background.h_over_h0(cosmo,a=1/(1+z))

    if typ=='linear':
        Plin_ell = [Plin(l,z) for l in Ell]
        return Hz/c_ls/chi(z)**2*ccl.correlations.correlation(cosmo,ell=Ell,C_ell=Plin_ell,theta=theta,type='NN',method='Legendre')
    
    if typ=='NL':
        PNL_ell = [PNL(l,z) for l in Ell]
        return Hz/c_ls/chi(z)**2*ccl.correlations.correlation(cosmo,ell=Ell,C_ell=PNL_ell,theta=theta,type='NN',method='Legendre')
