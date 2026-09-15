import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)


def K_r_from_M1450(z, alpha_nu=-0.5, lambda_r_eff=4500.0, lambda_0=1450.0):
    z = np.asarray(z, dtype=float)
    K_z_term = -2.5 * (1.0 + alpha_nu) * np.log10(1.0 + z)
    K_color_term = 2.5 * alpha_nu * np.log10(lambda_r_eff / lambda_0)
    return K_z_term + K_color_term


def dn_dzdmdOmega(z, m, lf, Kcorr=K_r_from_M1450):
    z = np.asarray(z, dtype=float)
    m = np.asarray(m, dtype=float)
    DL = cosmo.luminosity_distance(z).to(u.pc).value
    DM = 5.0 * np.log10(DL / 10.0)
    K = 0.0 if Kcorr is None else Kcorr(z)
    M = m - DM - K
    Phi = lf(M, z)
    dV_dz_dOmega = cosmo.differential_comoving_volume(z).to(u.Mpc**3 / u.sr).value
    return dV_dz_dOmega * Phi


def qso_lf_M1450(M, z, zp=2.2):
    """
    Double power-law QSO LF in M_1450, with a PLE-style break-magnitude
    evolution about a pivot redshift zp (Croom+09 / PD16 functional form):

        M_star(z) = M_star(zp) - 2.5*(k1*(z-zp) + k2*(z-zp)**2)

    NOTE: phi_star, M_star_zp, alpha, beta, k1, k2 below are PLACEHOLDER
    coefficients illustrating the correct functional form -- replace with
    verified values from a specific paper's fitted table (e.g.
    Palanque-Delabrouille et al. 2016, Table 4; or Croom et al. 2009)
    before using this for science. Do not trust the digits as-is.
    """
    z = np.asarray(z, dtype=float)

    phi_star = 2.66e-7      # Mpc^-3 mag^-1, PLACEHOLDER (z=0 normalization)
    M_star_zp = -25.36      # M_1450 at z = zp, PLACEHOLDER
    alpha = -1.30           # faint-end slope, PLACEHOLDER
    beta = -3.11            # bright-end slope, PLACEHOLDER
    k1 = 1.0                # linear evolution coeff, PLACEHOLDER
    k2 = -=0.2               # quadratic evolution coeff, PLACEHOLDER

    M_star = M_star_zp# - 2.5 * (k1 * (z - zp) + k2 * (z - zp) ** 2)

    return phi_star / (
        10 ** (0.4 * (alpha + 1) * (M - M_star))
        + 10 ** (0.4 * (beta + 1) * (M - M_star))
    )


def dN_dzdm(z_arr, m_arr, area_sr, lf=qso_lf_M1450, Kcorr=K_r_from_M1450):
    Zg, Mg = np.meshgrid(z_arr, m_arr, indexing='ij')
    return dn_dzdmdOmega(Zg, Mg, lf, Kcorr=Kcorr) * area_sr


# --- Bin definitions ---
Z_EDGES = np.linspace(0, 6, 60)
MAG_EDGES = np.linspace(10, 28, 50)
Z_MID = 0.5 * (Z_EDGES[:-1] + Z_EDGES[1:])
MAG_MID = 0.5 * (MAG_EDGES[:-1] + MAG_EDGES[1:])
BINS = [Z_EDGES, MAG_EDGES]

# --- Evaluation grids: MUST cover the bin edges, with a little padding ---
zarr = np.linspace(1e-7, Z_EDGES[-1] + 0.05, 10001)
magarr = np.linspace(MAG_EDGES[0] - 0.5, MAG_EDGES[-1] + 0.5, 1000)

S_deg2 = 4 * np.pi * (180 / np.pi) ** 2
dNdzdm = dN_dzdm(zarr, magarr, area_sr=S_deg2 * (np.pi / 180) ** 2)

# Matrix of integrals: (z bins, magnitude bins)
Nzm = np.zeros((len(Z_EDGES) - 1, len(MAG_EDGES) - 1))
for i in range(len(Z_EDGES) - 1):
    zmask = (zarr >= Z_EDGES[i]) & (zarr < Z_EDGES[i + 1])
    for j in range(len(MAG_EDGES) - 1):
        mmask = (magarr >= MAG_EDGES[j]) & (magarr < MAG_EDGES[j + 1])
        dNdzdm_bin = dNdzdm[np.ix_(zmask, mmask)]
        int_mag = np.trapz(dNdzdm_bin, x=magarr[mmask], axis=1)
        Nzm[i, j] = np.trapz(int_mag, x=zarr[zmask])

nzm = Nzm / S_deg2

np.savez(
    '../photom_redshift_distribution/COSMOS_QSO_LF.npz',
    z_center=Z_MID,
    mag_center=MAG_MID,
    z_edges=Z_EDGES,
    mag_edges=MAG_EDGES,
    object_count=Nzm,
    surface_deg2=S_deg2,
    info_about_the_sample='QSO N(z,m) distribution based on QSO LF',
)