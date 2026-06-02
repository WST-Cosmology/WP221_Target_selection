import numpy as np

def bias_lbg(redshift, mag):
    #from https://arxiv.org/pdf/2106.09713
    def A(m): return -0.98 * (m-25) + 0.11
    def B(m): return 0.12 * (m-25) + 0.17
    return A(mag) * (1 + redshift) + B(mag) * (1 + redshift)**2

def bz_model_desiqso(redshift, mag):
    return 0.278 * ((1 + redshift)**2 - 6.565) + 2.393