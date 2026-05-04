import sys
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt


RA_min = 148
RA_max = 152.5
DEC_min = 0
DEC_max = 4.5

data_path = '/global/cfs/cdirs/desi/users/cpayerne/data_WP221_Target_selection/photometric_catalogs/'

cosmos_clauds_hsc = Table.read(data_path + 'COSMOS_11bands-SExtractor-Lephare.fits')

bins=[np.linspace(148, 152.5, 80), 
    np.linspace(0, 4.5, 100)]
h, x, y, _=plt.hist2d(cosmos_clauds_hsc['RA'], cosmos_clauds_hsc['DEC'],
           cmin=1, bins=bins)

print('Surface of COSMOS (deg2): ', len(h[h>1]) * (bins[0][1] - bins[0][0]) * (bins[1][1] - bins[1][0]))


xmm_clauds_hsc = Table.read(data_path + 'XMMLSS_11bands-SExtractor-Lephare.fits')

bins=[np.linspace(33.5, 38, 80), 
    np.linspace(-6.5, -3, 100)]
h, x, y, _=plt.hist2d(xmm_clauds_hsc['RA'], xmm_clauds_hsc['DEC'],
           cmin=1, bins=bins)

print('Surface of XMM (deg2): ', len(h[h>1]) * (bins[0][1] - bins[0][0]) * (bins[1][1] - bins[1][0]))

