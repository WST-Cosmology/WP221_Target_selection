import numpy as np
import fitsio
import sys
from astropy.table import Table, hstack, QTable

dat_redshift_eff_2h = QTable.read("../data/efficiency_xmm_2hours_photo_z.ecsv")
dat_redshift_eff_4h = QTable.read("../data/efficiency_xmm_4hours_photo_z.ecsv")

def success_rate(z_sel, which='2h'):
    if which=='2h': 
        dat_redshift_eff = dat_redshift_eff_2h
    elif which=='4h': 
        dat_redshift_eff = dat_redshift_eff_4h

    z_tab = dat_redshift_eff['photo_z']
    f_tab = dat_redshift_eff['eff']
    mask_nan = np.invert(np.isnan(f_tab))
    w = np.interp(z_sel, z_tab[mask_nan*(z_tab>2)], f_tab[mask_nan*(z_tab>2)], left=0, right=f_tab[mask_nan*(z_tab>2)][-1])
    return w