'''
Utility functions for CRS target selections and masking from Verdier et al 2025.
These functions are designed to handle various selection criteria and masking operations
for astronomical catalogues, particularly for the Legacy Survey (LS) and other related datasets.

'''

import numpy as np
from astropy.coordinates import SkyCoord
# from astropy.wcs import WCS
# from astropy.io import fits
# from astropy.cosmology import default_cosmology
from astropy.table import Table
# import healpy as hp
# import astropy.cosmology.units as cu
import pandas as pd
from astropy import units as u

def flux_to_mag_LS(flux, flux_mw_transmission=None):

    if flux_mw_transmission is None:
        mag = 22.5 - 2.5 * np.log10(flux)
    else:
        mag = 22.5 - 2.5 * np.log10(flux / flux_mw_transmission)

    return mag


def selection_mask(ra, dec, field_limits, pixel_mask):
    """
    Apply selection based on field limits and a pixel mask to catalogue coordinates.

    The selection criteria include:
      - RA and Dec within the specified field limits.
      - The pixel mask value must be 0 (indicating an unmasked region).

    Parameters:
        ra (array_like): Array of right ascension values in degrees.
        dec (array_like): Array of declination values in degrees.
        field_limits (list): Field limits as [ra_min, ra_max, dec_min, dec_max].
        pixel_mask (PixMask): PixMask object with an apply_mask method.

    Returns:
        numpy.ndarray: Boolean array indicating which entries pass the selection.
    """
    ra_min, ra_max, dec_min, dec_max = field_limits
    sel = ((ra - ra_min) % 360 <= (ra_max - ra_min) % 360) & \
          (dec_min <= dec) & (dec <= dec_max) & \
          (pixel_mask.apply_mask(ra, dec) == 0)
    return sel


def LS_masking(data, selections=None):
    """
    Apply Legacy Survey (LS) masking on the dataset based on specified selection bits.

    This function creates a boolean mask based on the MASKBITS and WISEMASK_W1 bits provided 
    in the selections dictionary. It iterates through given bits and masks the data accordingly.

    Parameters:
        data (numpy.recarray or structured array): Input data containing mask columns.
        selections (dict, optional): Dictionary with keys 'MASKBITS' and 'WISEMASK_W1' that list 
                                     the bit positions to check. Defaults to None.

    Returns:
        numpy.ndarray: Boolean mask array where True indicates unmasked (selected) data.
    """
    MASKBITS = selections['MASKBITS']
    WISEMASK_W1 = selections['WISEMASK_W1']
    print('MASKBITS, WISEMASK_W1', MASKBITS, WISEMASK_W1, flush=True)
    mask = np.ones(len(data), dtype=bool)

    if len(MASKBITS) > 0:
        print('MASKBITS', flush=True)
        for bit in MASKBITS:
            print('bit', bit, flush=True)
            mask &= ~(data['MASKBITS'] & 2 ** bit > 0)

    if len(WISEMASK_W1) > 0:
        print('WISEMASK_W1', flush=True)
        for bit in WISEMASK_W1:
            print('bit', bit, flush=True)
            mask &= ~(data['WISEMASK_W1'] & 2 ** bit > 0)
    return mask


def LS_masking_BG(data, selections=None):
    """
    Apply Legacy Survey masking for BG data.

    For BG data, only the MASKBITS are applied, ignoring the WISEMASK_W1 mask.

    Parameters:
        data (numpy.recarray or structured array): Input data containing the MASKBITS column.
        selections (dict, optional): Dictionary with key 'MASKBITS'. Defaults to None.

    Returns:
        numpy.ndarray: Boolean mask array for BG data.
    """
    MASKBITS = selections['MASKBITS']
    print('MASKBITS', MASKBITS, flush=True)
    mask = np.ones(len(data), dtype=bool)

    if len(MASKBITS) > 0:
        print('MASKBITS', flush=True)
        for bit in MASKBITS:
            print('bit', bit, flush=True)
            mask &= ~(data['MASKBITS'] & 2 ** bit > 0)
    return mask


def LS_final_masking(data, selections=None):
    """
    Apply final Legacy Survey masking based on the data type selection.

    Different masking criteria are applied depending on the selection type:
      - 'BG' or 's0801': Uses a set of MASKBITS.
      - 'LRG' or 's0802': Uses MASKBITS and also filters based on WISEMASK_W1.
      - Otherwise, applies a default MASKBITS mask.

    Parameters:
        data (numpy.recarray or structured array): Input data containing MASKBITS and possibly WISEMASK_W1.
        selections (str, optional): A string indicating the type of masking to apply.
                                    Defaults to None.

    Returns:
        numpy.ndarray: Boolean mask array after final masking.
    """
    if selections == 'BG' or selections == 's0801':
        mask = ~(data['MASKBITS'] & 2 ** 1 > 0)
        mask &= ~(data['MASKBITS'] & 2 ** 12 > 0)
        mask &= ~(data['MASKBITS'] & 2 ** 13 > 0)
        mask &= ~(data['MASKBITS'] & 2 ** 11 > 0)

    elif selections == 'LRG' or selections == 's0802':
        mask = ~(data['MASKBITS'] & 2 ** 1 > 0)
        mask &= ~(data['MASKBITS'] & 2 ** 12 > 0)
        mask &= ~(data['MASKBITS'] & 2 ** 13 > 0)
        mask &= ~(data['MASKBITS'] & 2 ** 11 > 0)
        mask &= ~(data['WISEMASK_W1'] > 0)
    else:
        mask = ~(data['MASKBITS'] & 2 ** 1 > 0)
    return mask


def crs_bg_selection(cat, mag_r_lim=19.25, fib_mag_r_col='fib_mag_r', mag_r_col='mag_r',
                     mag_g_col='mag_g', mag_z_col='mag_z', fibtotal_mag_col='fibtot_mag_r',
                     raw_mag_r_col='raw_mag_r', if_mask: bool = False, r_z_cut_ratio=0.9,
                     r_z_cut=0.2, colour_cut_v2: bool = True):

    """
    Apply selection criteria to a catalogue of astronomical objects.
    The selection is based on various magnitude and colour cuts, as well as
    additional criteria such as the presence of Gaia data and the number of observations.
    The function also applies a mask to the data if specified.
    Parameters:
        cat (astropy.table.Table): Input catalogue containing astronomical data.
        mag_r_lim (float): Magnitude limit for the r-band. Default is 19.25.
        fib_mag_r_col (str): Column name for fiber magnitude in r-band. Default is 'fib_mag_r'.
        mag_r_col (str): Column name for magnitude in r-band. Default is 'mag_r'.
        mag_g_col (str): Column name for magnitude in g-band. Default is 'mag_g'.
        mag_z_col (str): Column name for magnitude in z-band. Default is 'mag_z'.
        fibtotal_mag_col (str): Column name for total fiber magnitude in r-band. Default is 'fibtot_mag_r'.
        raw_mag_r_col (str): Column name for raw magnitude in r-band. Default is 'raw_mag_r'.
        if_mask (bool): If True, apply the LS final masking. Default is False.
        r_z_cut_ratio (float): Ratio for the r-z colour cut. Default is 0.9.
        r_z_cut (float): Cut value for the r-z colour cut. Default is 0.2.
        colour_cut_v2 (bool): If True, apply the new colour cut. Default is True.
    Returns:
        numpy.ndarray: Boolean mask array indicating which entries pass the selection.
    """
    
    if not isinstance(cat, pd.DataFrame):
        try:
            cat = cat.to_pandas()
        except AttributeError:
            raise TypeError("Input data must be a pandas DataFrame or an astropy table.")
        
    if mag_r_col in cat.colnames:
        mag_r = cat[mag_r_col]
    else:
        mag_r = flux_to_mag_LS(cat['FLUX_R'],
                               flux_mw_transmission=cat['MW_TRANSMISSION_R'])
        cat[mag_r_col] = mag_r

    if fib_mag_r_col in cat.colnames:
        fib_mag_r = cat[fib_mag_r_col]

    else:
        fib_mag_r = flux_to_mag_LS(cat['FIBERFLUX_R'],
                                   flux_mw_transmission=cat['MW_TRANSMISSION_R'])
        cat[fib_mag_r_col] = fib_mag_r

    if mag_g_col in cat.colnames:
        mag_g = cat[mag_g_col]
    else:
        mag_g = flux_to_mag_LS(cat['FLUX_G'],
                               flux_mw_transmission=cat['MW_TRANSMISSION_G'])
        cat[mag_g_col] = mag_g

    if mag_z_col in cat.colnames:
        mag_z = cat[mag_z_col]
    else:
        mag_z = flux_to_mag_LS(cat['FLUX_Z'],
                               flux_mw_transmission=cat['MW_TRANSMISSION_Z'])
        cat[mag_z_col] = mag_z

    if fibtotal_mag_col in cat.colnames:
        fibtot_mag_r = cat[fibtotal_mag_col]
    else:
        fibtot_mag_r = flux_to_mag_LS(cat['FIBERTOTFLUX_R'],
                                      flux_mw_transmission=cat['MW_TRANSMISSION_R'])
        cat[fibtotal_mag_col] = fibtot_mag_r

    if raw_mag_r_col in cat.colnames:
        raw_mag_r = cat[raw_mag_r_col]
    else:
        raw_mag_r = flux_to_mag_LS(cat['FLUX_R'],
                                   flux_mw_transmission=cat['MW_TRANSMISSION_R'])
        cat[raw_mag_r_col] = raw_mag_r

    sel = mag_r < mag_r_lim

    if 'REF_CAT' in cat.colnames:
        not_in_gaia = cat['REF_CAT'] == '  '

        sel &= ((cat['GAIA_PHOT_G_MEAN_MAG'] - raw_mag_r) > 0.6) & ~not_in_gaia
        sel |= not_in_gaia

    if if_mask:
        sel &= LS_final_masking(cat, selections='BG')

    # Fibre magnitude cut
    mask = mag_r < 17.8
    mask &= fib_mag_r < (22.9 + (mag_r - 17.8))
    mask1 = (mag_r > 17.8) & (mag_r < 20)
    mask1 &= fib_mag_r < 22.9
    sel &= (mask | mask1)

    # cut spurious objects
    mask = (mag_g-mag_r < -1) | (mag_g-mag_r > 4)
    mask |= (mag_r-mag_z < -1) | (mag_r-mag_z > 4)
    sel &= ~mask

    if all(key in cat.colnames for key in ['NOBS_G', 'NOBS_R', 'NOBS_Z']):
        mask = (cat['NOBS_G'] == 0) | (
            cat['NOBS_R'] == 0) | (cat['NOBS_Z'] == 0)
        sel &= ~mask

    mask = (mag_r > 12) & (fibtot_mag_r < 15)
    sel &= ~mask

    # Additional photometric cut
    if colour_cut_v2:
        sel &= mag_r- mag_z < (0.93*(mag_g-mag_r)-0.27)
        sel &= mag_r- mag_z > (0.4*(mag_g-mag_r)+0.07)
        sel &= mag_g - mag_r > 0.95


    else:
        sel &= mag_r-mag_z < (r_z_cut_ratio*(mag_g-mag_r)-r_z_cut)

    if all(key in cat.colnames for key in ['PARALLAX', 'PARALLAX_IVAR', 'PMRA', 'PMRA_IVAR', 'PMDEC', 'PMDEC_IVAR']):

        snr_par = np.abs(cat['PARALLAX']*np.sqrt(cat['PARALLAX_IVAR']))
        snr_pmra = np.abs(cat['PMRA']*np.sqrt(cat['PMRA_IVAR']))
        snr_pmdec = np.abs(cat['PMDEC']*np.sqrt(cat['PMDEC_IVAR']))
        mask = cat['PARALLAX'] != 0
        mask &= (snr_pmra > 3) | (snr_pmra > 3) | (snr_pmdec > 3)
        sel &= ~mask

    return sel


def crs_bg_sel_v2_nofib(cat, mag_r_limit=19.25, mag_r_col='mag_r',
                        mag_g_col='mag_g', mag_z_col='mag_z'):
        """
        Apply selection criteria of the CRS BGS catalogue without fiber magnitude cut.
        The selection is based on various magnitude and colour cuts.
        Parameters:
            cat (astropy.table.Table): Input catalogue containing astronomical data.
            mag_r_limit (float): Magnitude limit for the r-band. Default is 19.25.
            mag_r_col (str): Column name for magnitude in r-band. Default is 'mag_r'.
            mag_g_col (str): Column name for magnitude in g-band. Default is 'mag_g'.
            mag_z_col (str): Column name for magnitude in z-band. Default is 'mag_z'.
            
        Returns:
            numpy.ndarray: Boolean mask array indicating which entries pass the selection.

        """
            
        mag_r = cat[mag_r_col]
        mag_g = cat[mag_g_col]
        mag_z = cat[mag_z_col]

        sel = cat[mag_r_col] < mag_r_limit
        sel &= mag_r- mag_z < (0.93*(mag_g-mag_r)-0.27)
        sel &= mag_r- mag_z > (0.4*(mag_g-mag_r)+0.07)
        sel &= mag_g - mag_r > 0.95

        return sel


def desi_bgs_sel(cat, mag_r_lim=19.5, fib_mag_r_col='fib_mag_r', mag_r_col='mag_r',
                     mag_g_col='mag_g', mag_z_col='mag_z', fibtotal_mag_col='fibtot_mag_r',
                     raw_mag_r_col='raw_mag_r'):
    """

    Apply selection criteria to a catalogue of astronomical objects for DESI BGS.
    The selection is based on various magnitude and colour cuts, as well as
    additional criteria such as the presence of Gaia data and the number of observations.
    
    Parameters:
        cat (astropy.table.Table): Input catalogue containing astronomical data.
        mag_r_lim (float): Magnitude limit for the r-band. Default is 19.5.
        fib_mag_r_col (str): Column name for fiber magnitude in r-band. Default is 'fib_mag_r'.
        mag_r_col (str): Column name for magnitude in r-band. Default is 'mag_r'.
        mag_g_col (str): Column name for magnitude in g-band. Default is 'mag_g'.
        mag_z_col (str): Column name for magnitude in z-band. Default is 'mag_z'.
        fibtotal_mag_col (str): Column name for total fiber magnitude in r-band. Default is 'fibtot_mag_r'.
        raw_mag_r_col (str): Column name for raw magnitude in r-band. Default is 'raw_mag_r'.
    Returns:
        numpy.ndarray: Boolean mask array indicating which entries pass the selection.
    """
    if not isinstance(cat, pd.DataFrame):
        try:
            cat = cat.to_pandas()
        except AttributeError:
            raise TypeError("Input data must be a pandas DataFrame or an astropy table.")
    
    if mag_r_col in cat.colnames:
        mag_r = cat[mag_r_col]
    else:
        mag_r = flux_to_mag_LS(cat['FLUX_R'],
                               flux_mw_transmission=cat['MW_TRANSMISSION_R'])
        cat[mag_r_col] = mag_r

    if fib_mag_r_col in cat.colnames:
        fib_mag_r = cat[fib_mag_r_col]

    else:
        fib_mag_r = flux_to_mag_LS(cat['FIBERFLUX_R'],
                                   flux_mw_transmission=cat['MW_TRANSMISSION_R'])
        cat[fib_mag_r_col] = fib_mag_r

    if mag_g_col in cat.colnames:
        mag_g = cat[mag_g_col]
    else:
        mag_g = flux_to_mag_LS(cat['FLUX_G'],
                               flux_mw_transmission=cat['MW_TRANSMISSION_G'])
        cat[mag_g_col] = mag_g

    if mag_z_col in cat.colnames:
        mag_z = cat[mag_z_col]
    else:
        mag_z = flux_to_mag_LS(cat['FLUX_Z'],
                               flux_mw_transmission=cat['MW_TRANSMISSION_Z'])
        cat[mag_z_col] = mag_z

    if fibtotal_mag_col in cat.colnames:
        fibtot_mag_r = cat[fibtotal_mag_col]
    else:
        fibtot_mag_r = flux_to_mag_LS(cat['FIBERTOTFLUX_R'],
                                      flux_mw_transmission=cat['MW_TRANSMISSION_R'])
        cat[fibtotal_mag_col] = fibtot_mag_r

    if raw_mag_r_col in cat.colnames:
        raw_mag_r = cat[raw_mag_r_col]
    else:
        raw_mag_r = flux_to_mag_LS(cat['FLUX_R'],
                                   flux_mw_transmission=cat['MW_TRANSMISSION_R'])
        cat[raw_mag_r_col] = raw_mag_r

    sel = mag_r < mag_r_lim

    if 'REF_CAT' in cat.colnames:
        not_in_gaia = cat['REF_CAT'] == '  '

        sel &= ((cat['GAIA_PHOT_G_MEAN_MAG'] - raw_mag_r) > 0.6) & ~not_in_gaia
        sel |= not_in_gaia


    # Fibre magnitude cut
    mask = mag_r < 17.8
    mask &= fib_mag_r < (22.9 + (mag_r - 17.8))
    mask1 = (mag_r > 17.8) & (mag_r < 20)
    mask1 &= fib_mag_r < 22.9
    sel &= (mask | mask1)

    # cut spurious objects
    mask = (mag_g-mag_r < -1) | (mag_g-mag_r > 4)
    mask |= (mag_r-mag_z < -1) | (mag_r-mag_z > 4)
    sel &= ~mask

    if all(key in cat.colnames for key in ['NOBS_G', 'NOBS_R', 'NOBS_Z']):
        mask = (cat['NOBS_G'] == 0) | (
            cat['NOBS_R'] == 0) | (cat['NOBS_Z'] == 0)
        sel &= ~mask

    mask = (mag_r > 12) & (fibtot_mag_r < 15)
    sel &= ~mask


    if all(key in cat.colnames for key in ['PARALLAX', 'PARALLAX_IVAR', 'PMRA', 'PMRA_IVAR', 'PMDEC', 'PMDEC_IVAR']):

        snr_par = np.abs(cat['PARALLAX']*np.sqrt(cat['PARALLAX_IVAR']))
        snr_pmra = np.abs(cat['PMRA']*np.sqrt(cat['PMRA_IVAR']))
        snr_pmdec = np.abs(cat['PMDEC']*np.sqrt(cat['PMDEC_IVAR']))
        mask = cat['PARALLAX'] != 0
        mask &= (snr_pmra > 3) | (snr_pmra > 3) | (snr_pmdec > 3)
        sel &= ~mask

    return sel


def crs_lrg_sel(data, shift = 0.06, fib_mag_z_col='fib_mag_z', mag_r_col='mag_r',
                mag_g_col='mag_g', mag_z_col='mag_z', mag_w1_col='mag_w1', if_mask: bool = False):


    if not isinstance(data, pd.DataFrame):
        try:
            data = data.to_pandas()
        except AttributeError:
            raise TypeError("Input data must be a pandas DataFrame or an astropy table.")
        
    if mag_g_col in data.columns():
        gmag = data[mag_g_col]
    else:
        gmag = flux_to_mag_LS(data['FLUX_G'],
                                 flux_mw_transmission=data['MW_TRANSMISSION_G'])
    
    if mag_r_col in data.columns():
        rmag = data[mag_r_col]
    else:
        rmag = flux_to_mag_LS(data['FLUX_R'],
                               flux_mw_transmission=data['MW_TRANSMISSION_R'])

    if mag_z_col in data.columns():
        zmag = data[mag_z_col]
    else:
        zmag = flux_to_mag_LS(data['FLUX_Z'],
                               flux_mw_transmission=data['MW_TRANSMISSION_Z'])
    
    if mag_w1_col in data.columns():
        w1mag = data[mag_w1_col]
    else:
        w1mag = flux_to_mag_LS(data['FLUX_W1'],
                               flux_mw_transmission=data['MW_TRANSMISSION_W1'])
    
    if fib_mag_z_col in data.columns():
        fiberz = data[fib_mag_z_col]
    else:
        fiberz = flux_to_mag_LS(data['FIBERFLUX_Z'],
                                 flux_mw_transmission=data['MW_TRANSMISSION_Z'])
    
        
    if fib_mag_z_col in data.columns():
        fiberz = data[fib_mag_z_col]
    else:
        fiberz = flux_to_mag_LS(data['FIBERFLUX_Z'],
                                 flux_mw_transmission=data['MW_TRANSMISSION_Z'])

    alpha_new = 1.8
    beta_new = 17.14 - shift
    gamma_new = 16.33 - 1.5*shift

    mask_tot_ls = fiberz < 21.6
    mask_tot_ls &= (zmag - w1mag) > 0.8*(rmag - zmag) - 0.6
    mask_tot_ls &= ((gmag - w1mag) > 2.9) | ((rmag - w1mag) > 1.8)  # pq 2.97 et pas 2.9
    mask_tot_ls &= (rmag - w1mag > alpha_new*(w1mag - beta_new)) & ((rmag - w1mag > w1mag - gamma_new))

    # I would remove this cut
    #mask_tot_ls &= (photz > 0) & (photz < 1.4)

    # cut according Nobs,
    if 'NOBS_G' in data.columns():
        mask = (data['NOBS_G'] == 0) | (data['NOBS_R'] == 0) | (data['NOBS_Z'] == 0)
        mask_tot_ls &= ~mask
    
    # cut according bad photometry
    mask = (data['FLUX_IVAR_R'] < 0) | (data['FLUX_IVAR_Z'] < 0) | (data['FLUX_IVAR_W1'] < 0) | (data['FLUX_G'] < 0)
    mask_tot_ls &= ~mask

    # Additionnal cuts from the DESI LRG TS
    mask_tot_ls &= fiberz > 17.5

    if 'GAIA_PHOT_G_MEAN_MAG' in data.columns():
        mask_tot_ls &= (data['GAIA_PHOT_G_MEAN_MAG'] > 18) | (data['GAIA_PHOT_G_MEAN_MAG'] == 0) 

    if if_mask:
        mask_tot_ls &= LS_final_masking(data, selections='LRG')


    if all(key in data.colnames for key in ['PARALLAX', 'PARALLAX_IVAR', 'PMRA', 'PMRA_IVAR', 'PMDEC', 'PMDEC_IVAR']):
        snr_par = np.abs(data['PARALLAX']*np.sqrt(data['PARALLAX_IVAR']))
        snr_pmra = np.abs(data['PMRA']*np.sqrt(data['PMRA_IVAR']))
        snr_pmdec = np.abs(data['PMDEC']*np.sqrt(data['PMDEC_IVAR']))
        mask = data['PARALLAX'] != 0
        mask &= (snr_par > 3) | (snr_pmra > 3) | (snr_pmdec > 3) 
        mask_tot_ls &= ~mask
    
    return mask_tot_ls
	

def get_galactic_coords(ra, dec):
    """
    Convert equatorial coordinates (RA, Dec) to galactic coordinates (l, b).

    Parameters:
        ra (array_like): Right ascension in degrees.
        dec (array_like): Declination in degrees.

    Returns:
        tuple: Galactic coordinates (l, b) in degrees.
    """
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    galactic_coords = coords.galactic

    return galactic_coords.l.deg, galactic_coords.b.deg