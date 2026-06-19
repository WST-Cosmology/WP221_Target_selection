from astropy.table import Table, QTable, hstack, vstack, join
from astropy.coordinates import SkyCoord, match_coordinates_3d, match_coordinates_sky
import astropy.units as un
import numpy as np
import copy

def match_nearest_neghbor(base_catalog = None, ra_base = 'ra_base', dec_base = 'dec_base', label_base = '_base',
                          target_catalog = None, ra_target = 'ra_target', dec_target = 'dec_target',label_target = '_target', 
                          max_sep_arcsec = 1, force_base_format = False):

    base_catalog['id_base'] = np.arange(len(base_catalog))
    target_catalog['id_target'] = np.arange(len(target_catalog))
    if force_base_format == False:
        dat_tot = Table()
        
        base_SkyCoord = SkyCoord(ra=np.array(base_catalog[ra_base])*un.deg, dec=np.array(base_catalog[dec_base])*un.deg)
        target_SkyCoord = SkyCoord(ra=np.array(target_catalog[ra_target])*un.deg, dec=np.array(target_catalog[dec_target])*un.deg)

        idx, sep2d, sep3d = match_coordinates_sky(target_SkyCoord, base_SkyCoord, nthneighbor=1,storekdtree='kdtree_3d')
        r"""
        idx: indexes of the base catalog objects matching each target object, closest object from target
        sep2d: correspondind separation distance
        """

        mask_idx_unique_match = np.array([False for i in range(len(idx))])

        for k, idx_match in enumerate(np.unique(idx)):
            subsample_match_sep2d = sep2d[idx == idx_match]
            min_sep2d = np.min(subsample_match_sep2d*un.arcsec)
            mask = (idx == idx_match) & (sep2d*un.arcsec == min_sep2d) & (sep2d < max_sep_arcsec*un.arcsec)
            mask_idx_unique_match[(idx == idx_match) & (sep2d*un.arcsec == min_sep2d) & (sep2d < max_sep_arcsec*un.arcsec)] = True

        target_cut = target_catalog[mask_idx_unique_match]
        base_cut = base_catalog[idx[mask_idx_unique_match]]
        
        for name in base_catalog.colnames:
            dat_tot[name + label_base] = base_cut[name]
        
        for name in target_catalog.colnames:
            dat_tot[name+ label_target] = target_cut[name]
                            
        return dat_tot
