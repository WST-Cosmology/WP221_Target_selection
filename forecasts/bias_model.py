import numpy as np
import pyccl as ccl

h=0.677
deltac=1.686
H0=100*h
c_ls=300*10**3
nlim=10000
n_s=0.968

cosmo = ccl.Cosmology(Omega_c=0.27, Omega_b=0.045, h=h, A_s=2.1e-9, n_s=n_s,transfer_function='boltzmann_camb')

def D(z):
    a = 1.0 / (1.0 + np.asarray(z))
    return cosmo.growth_factor(a) / cosmo.growth_factor(1.0)

def bias_bg(redshift, mag):
    return 1.34 * 1 / D(redshift)

def bias_lrg(redshift, mag):
    return 0.209 * (1 + redshift)**2 + 1.415 

def bias_elg(redshift, mag):
    return 0.84 * 1 / D(redshift)

def bias_magmax(redshift, mag):
    return 0.87+0.32*(25.12-mag)*redshift**1.37

def bias_qso(redshift, mag):
    return 0.237 * (1 + redshift)**2 + 0.771
    #return 0.237 * ((1 + redshift)**2 - 6.565) + 2.328

def bias_lbg(redshift, mag):
    #from https://arxiv.org/pdf/2106.09713
    def A(m): return -0.98 * (m-25) + 0.11
    def B(m): return 0.12 * (m-25) + 0.17
    return A(mag) * (1 + redshift) + B(mag) * (1 + redshift)**2


def linear_bias(redshift, mag, tracer = 'BG_faint'):

    if tracer == 'BG_faint' or tracer == 'BG_bright': return bias_bg(redshift, mag)
    if tracer == 'LRG': return bias_lrg(redshift, mag)
    if tracer == 'ELG': return bias_elg(redshift, mag)
    if 'MagMax' in tracer.split('_'): return bias_magmax(redshift, mag)
    if tracer == 'QSO_qlf': return bias_qso(redshift, mag)
    if tracer == 'LBGu': return bias_lbg(redshift, mag)
    if tracer == 'LBGg': return bias_lbg(redshift, mag)
    if tracer == 'LBGr': return bias_lbg(redshift, mag)
        
    return np.ones(len(redshift))