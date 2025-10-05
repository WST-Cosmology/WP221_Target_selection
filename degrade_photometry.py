import os
import numpy as np
from astropy.table import Table, vstack
import pandas as pd
import math
from scipy.special import erf

def get_mag_from_flux(fluxs):
    return 22.5 - 2.5 * np.log10(fluxs)

# magerr = 2.5/ln(10) * fluxerr/flux 
def get_mag_magerr_from_flux_fluxivar(fluxs, flux_ivars):
    mags = get_mag_from_flux(fluxs)
    magerrs = np.log(10) / 2.5 / (fluxs * np.sqrt(flux_ivars))
    return mags, magerrs

def get_flux_from_mag(mags):
    return 10 ** (-0.4 * (mags - 22.5))

# magerr = 2.5/ln(10) * fluxerr/flux
#        = 2.5/ln(10) * (flux_nsigma_depth/nsigma)/flux
#        = 2.5/ln(10)/nsigma * 10**(0.4*(mag-depth))
def get_model_magerrs(nsigma_depth, mags, depth):
    return 2.5 / np.log(10) / nsigma_depth * 10 ** (0.4 * (mags - depth))

"""
err_orig ** 2 + err_add ** 2 = err_shallow ** 2
err_add ** 2 = err_shallow ** 2 - err_orig ** 2
err_add = np.sqrt(err_shallow ** 2 - err_orig ** 2)
"""
def get_add_magerrs(orig_mags, nsigma_depth, orig_depth, shallow_depth):
    orig_model_magerrs = get_model_magerrs(nsigma_depth, orig_mags, orig_depth)
    shallow_model_magerrs = get_model_magerrs(nsigma_depth, orig_mags, shallow_depth)
    return np.sqrt(shallow_model_magerrs ** 2 - orig_model_magerrs ** 2)

def get_shallow_mags(orig_mags, orig_magerrs, nsigma_depth, orig_depth, shallow_depth):
    add_magerrs = get_add_magerrs(orig_mags, nsigma_depth, orig_depth, shallow_depth)
    shallow_mags = orig_mags.copy()
    shallow_magerrs = orig_magerrs.copy()
    sel = (orig_magerrs > 0) & (orig_magerrs < 1000)
#     shallow_mags[sel] = orig_mags[sel] + np.random.normal(size=sel.sum()) * add_magerrs[sel]
    """
    em_sh ** 2 = em_orig ** 2 + em_add ** 2
    em_add = 1.0857 * ef_add / f
    ef_add = f * em_add / 1.0857
    """
    shallow_fluxs = get_flux_from_mag(orig_mags) 
    add_fluxerrs = shallow_fluxs * add_magerrs / (2.5 / np.log(10))
    shallow_fluxs[sel] += np.random.normal(size=sel.sum()) * add_fluxerrs[sel]
    shallow_mags = get_mag_from_flux(shallow_fluxs)
    shallow_magerrs[sel] = np.sqrt(orig_magerrs[sel] ** 2 + add_magerrs[sel] ** 2)
    return shallow_mags, shallow_magerrs

def degrade(table, bands, orig_depths, shallow_depths, suff):
    np_rd_seed = 1232
    nsigma=5
    np.random.seed(np_rd_seed)
    d = table
    d.meta["RANDSEED"] = np_rd_seed
    for band, orig_depth, shallow_depth in zip(bands, orig_depths, shallow_depths):
        d.meta["{}_ORIG".format(band)] = orig_depth
        shallow_mags, shallow_magerrs = get_shallow_mags(d[band], d["{}_err".format(band)], nsigma,
        orig_depth,shallow_depth,)
        mag_err_deep = d[f"{band}_err"]
        d["{}_{}".format(band,suff)] = shallow_mags
        d["{}_{}_err".format(band,suff)] = shallow_magerrs
        d["{}_{}_adderr".format(band,suff)] = get_add_magerrs(d[band], nsigma, orig_depth, shallow_depth)
    return d

def gaussian_integral_from_a(mean, std, a):
    """
    Computes ∫_a^∞ Normal(μ, σ) dx
    """
    z = (a - mean) / (std * math.sqrt(2))
    return 0.5 * (1-erf(z))

def Probability_detect_minput(m_input, m_depth, Nsigma=1):

    F_input = 10**(-0.4 * (m_input-22.5))
    F_depth = 10**(-0.4 * (m_depth-22.5))
    sigma_F = F_depth / 5
    lowerlimit = Nsigma * sigma_F
    # Detection probability
    P_detect = gaussian_integral_from_a(F_input, sigma_F, lowerlimit)

    return P_detect
