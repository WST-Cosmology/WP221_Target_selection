import numpy as np

def bias_bg(redshift, mag):
    return 1.5 * np.ones(len(redshift))

def bias_lrg(redshift, mag):
    return 0.209 * (1 + redshift)**2 + 1.415 

def bias_elg(redshift, mag):
    return 1.2 * np.ones(len(redshift))

def bias_magmax(redshift, mag):
    return 1.2 * np.ones(len(redshift))

def bias_lbg(redshift, mag):
    #from https://arxiv.org/pdf/2106.09713
    def A(m): return -0.98 * (m-25) + 0.11
    def B(m): return 0.12 * (m-25) + 0.17
    return A(mag) * (1 + redshift) + B(mag) * (1 + redshift)**2

def bias_qso(redshift, mag):
    return 0.237 * (1 + redshift)**2 + 0.771