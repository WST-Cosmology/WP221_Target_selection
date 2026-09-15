import numpy as np
import math
from scipy.integrate import quad
from scipy.interpolate import interp1d
from functools import partial

import camb
import cosmology

nlim  = 10000

# ---------------------------------------------------------------------------
# Full cosmological parameter set
# ---------------------------------------------------------------------------
# Cosmo params beyond bias: [Mnu, H0, omega_c, omega_b, ns, As]
# Indices in the full Fisher matrix:
#   single tracer : [b,  Mnu, H0, omega_c, omega_b, ns, As]  -> 7 params
#   two tracers   : [ba, bb, Mnu, H0, omega_c, omega_b, ns, As] -> 8 params
# After marginalisation we return only the b / Mnu block.
# ---------------------------------------------------------------------------

COSMO_PARAMS = ['Mnu', 'H0', 'omega_c', 'omega_b', 'ns', 'As']

# Finite-difference step sizes for each cosmological parameter
DEFAULT_STEPS = {
    'Mnu'    : 0.06,
    'H0'     : 1.0,
    'omega_c': 0.002,
    'omega_b': 0.001,
    'ns'     : 0.01,
    'As'     : 2e-11,
}

# ---------------------------------------------------------------------------
# Planck Gaussian priors  (half-width = 1 sigma, from Table 1)
# Note: tau is not in COSMO_PARAMS (no effect on Pm) but kept here for
#       completeness if the user wants to add it as an external prior on Mnu.
# H0 prior derived from Planck h prior: sigma(H0) = 100 * sigma(h) = 0.54
# As prior: sigma(As) converted from sigma(ln As) = 0.015
#           -> sigma(As) = As_fid * 0.015 = 2.1e-9 * 0.015 ~ 3.1e-11
# ---------------------------------------------------------------------------
PLANCK_PRIORS = {
    'Mnu'    : 0.5,                        # eV  (very weak, essentially free)
    'H0'     : 100.0 * 0.0054,            # km/s/Mpc  (from sigma_h = 0.0054)
    'omega_c': 0.0012,
    'omega_b': 0.00015,
    'ns'     : 0.0042,
    'As'     : 2.1e-9 * 0.015,            # converted from sigma(ln As) = 0.015
}


# ---------------------------------------------------------------------------
# CAMB helpers
# ---------------------------------------------------------------------------

def _pars_from_dict(fid):
    """Build a CAMBparams object from a fiducial dict."""
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0    = fid['H0'],
        ombh2 = fid['omega_b'],
        omch2 = fid['omega_c'],
        mnu   = fid['Mnu'],
    )
    pars.InitPower.set_params(As=fid['As'], ns=fid['ns'])
    return pars


def _camb_pk(fid, kmax=10.0, npoints=300):
    """Return (kh, pk_z0) arrays from CAMB for a given fiducial dict."""
    pars = _pars_from_dict(fid)
    pars.set_matter_power(redshifts=[0.0], kmax=kmax * 1.1)
    results = camb.get_results(pars)
    kh, _, pk = results.get_matter_power_spectrum(
        minkh=1e-4, maxkh=kmax, npoints=npoints
    )
    return kh, pk[0]


def _build_all_derivs(cosmo, fid_dict, kmax=10.0, npoints=300, steps=None):
    """
    Build log-space interpolators for d ln Pm(k, z=0) / d theta
    for every theta in COSMO_PARAMS, via central finite differences with CAMB.
    """
    if steps is None:
        steps = DEFAULT_STEPS

    kh, pk_fid = _camb_pk(fid_dict, kmax=kmax, npoints=npoints)
    Pm_interp = interp1d(np.log(kh), np.log(pk_fid),
                         kind='cubic', fill_value='extrapolate')

    deriv_interps = {}
    for par in COSMO_PARAMS:
        step    = steps[par]
        fid_val = fid_dict[par]

        fid_p = dict(fid_dict); fid_p[par] = fid_val + step / 2
        fid_m = dict(fid_dict)
        fid_m[par] = max(fid_val - step / 2, 0.0) if par == 'Mnu' \
                     else fid_val - step / 2

        actual_step = fid_p[par] - fid_m[par]

        _, pk_p = _camb_pk(fid_p, kmax=kmax, npoints=npoints)
        _, pk_m = _camb_pk(fid_m, kmax=kmax, npoints=npoints)

        dlnPm = (np.log(pk_p) - np.log(pk_m)) / actual_step
        deriv_interps[par] = interp1d(np.log(kh), dlnPm,
                                      kind='cubic', fill_value='extrapolate')
        print(f'  derivative d ln Pm / d {par} done.')

    return deriv_interps, Pm_interp, kh


def _make_fid_dict(cosmo, Mnu_fid=0.06):
    """Extract a plain fiducial dict from the cosmology dict."""
    h = cosmo['h']
    return {
        'Mnu'    : Mnu_fid,
        'H0'     : 100.0 * h,
        'omega_c': cosmo['Omega_c'] * h**2,
        'omega_b': cosmo['Omega_b'] * h**2,
        'ns'     : cosmo['n_s'],
        'As'     : cosmo['A_s'],
    }


# ---------------------------------------------------------------------------
# Prior helper
# ---------------------------------------------------------------------------

def _apply_priors(F_full, prior_params, param_order):
    """
    Add Gaussian priors to the diagonal of F_full in-place.

    Parameters
    ----------
    F_full      : (N, N) Fisher matrix  — modified in place
    prior_params: list of str  — parameter names to apply priors to,
                                 e.g. ['omega_b', 'ns', 'As']
                                 Must be a subset of COSMO_PARAMS.
    param_order : list of str  — full ordered parameter list corresponding
                                 to rows/cols of F_full, e.g.
                                 ['b', 'Mnu', 'H0', 'omega_c', 'omega_b', 'ns', 'As']
    """
    for par in prior_params:
        if par not in PLANCK_PRIORS:
            raise ValueError(f"No Planck prior defined for '{par}'. "
                             f"Available: {list(PLANCK_PRIORS.keys())}")
        if par not in param_order:
            raise ValueError(f"Parameter '{par}' not found in param_order.")
        idx = param_order.index(par)
        sigma_prior = PLANCK_PRIORS[par]
        F_full[idx, idx] += 1.0 / sigma_prior**2
        print(f'  Planck prior on {par}: sigma = {sigma_prior:.4g}')


# ---------------------------------------------------------------------------
# Single-tracer Fisher matrix
# ---------------------------------------------------------------------------

def Mat_Fisher_1tracer_mnu_full(n, bg, z, Vsur, kmin, kmax,
                                 Mnu_fid=0.06,
                                 deriv_interps=None,
                                 Pm_interp=None,
                                 cosmo=None,
                                 Nk=200):
    """
    Full Fisher matrix for [b, Mnu, H0, omega_c, omega_b, ns, As].
    Returns the RAW (Npar x Npar) matrix — no priors, no marginalisation.
    Priors and marginalisation are applied once at the survey level.
    """
    Npar = 1 + len(COSMO_PARAMS)

    Dz = cosmology.D(z, cosmo)
    D0 = cosmology.D(0, cosmo)
    prefactor = Vsur / (4.0 * math.pi**2)

    k_nodes, wk = np.polynomial.legendre.leggauss(Nk)
    k  = 0.5 * (kmax - kmin) * k_nodes + 0.5 * (kmax + kmin)
    wk = 0.5 * (kmax - kmin) * wk

    Pm0    = np.exp(Pm_interp(np.log(k))) * (Dz / D0)**2
    P      = bg**2 * Pm0
    nP     = n * P
    weight = 2.0 * prefactor * 0.5 * (nP / (1.0 + nP))**2 * k**2 * wk

    D = np.zeros((Npar, Nk))
    D[0, :] = 2.0 / bg
    for ip, par in enumerate(COSMO_PARAMS):
        D[1 + ip, :] = deriv_interps[par](np.log(k))

    return np.einsum('ik,jk,k->ij', D, D, weight)


# ---------------------------------------------------------------------------
# Two-tracer Fisher matrix
# ---------------------------------------------------------------------------

def Mat_Fisher_2tracer_mnu_full(na, nb, ba, bb, zeff,
                                 Vsur, kmin, kmax,
                                 Mnu_fid=0.06,
                                 deriv_interps=None,
                                 Pm_interp=None,
                                 cosmo=None,
                                 Nk=200):
    """
    Full Fisher matrix for [ba, bb, Mnu, H0, omega_c, omega_b, ns, As].
    Returns the RAW (Npar x Npar) matrix — no priors, no marginalisation.
    Priors and marginalisation are applied once at the survey level.
    """
    Npar        = 2 + len(COSMO_PARAMS)
    param_order = ['ba', 'bb'] + COSMO_PARAMS

    z  = zeff
    Dz = cosmology.D(z, cosmo)
    D0 = cosmology.D(0, cosmo)
    prefactor = Vsur / (4.0 * np.pi**2)

    # Gauss-Legendre grid
    k_nodes, wk = np.polynomial.legendre.leggauss(Nk)
    k  = 0.5 * (kmax - kmin) * k_nodes + 0.5 * (kmax + kmin)
    wk = 0.5 * (kmax - kmin) * wk

    Pm0 = np.exp(Pm_interp(np.log(k))) * (Dz / D0)**2

    PA  = ba**2 * Pm0
    PB  = bb**2 * Pm0
    PAB = ba * bb * Pm0

    nA     = na * PA
    nB     = nb * PB
    nAnBX2 = na * nb * PAB**2
    one_nA = 1.0 + nA
    one_nB = 1.0 + nB
    det    = one_nA * one_nB - nAnBX2
    det2   = det**2

    Raa = (nA * one_nB / det)**2
    Rbb = (nB * one_nA / det)**2
    Rxx = na * nb * (one_nA * one_nB + nAnBX2) * PAB**2 / det2
    Rxa = na**2 * nb * one_nB * PAB**2 * PA / det2
    Rxb = nb**2 * na * one_nA * PAB**2 * PB / det2
    Rab = na**2 * nb**2 * PA * PB * PAB**2 / det2

    w = 2.0 * k**2 * wk

    dlnP_cosmo = np.array([deriv_interps[p](np.log(k)) for p in COSMO_PARAMS])

    DlnA  = np.vstack([np.full((1, Nk), 2.0/ba), np.zeros((1, Nk)), dlnP_cosmo])
    DlnB  = np.vstack([np.zeros((1, Nk)), np.full((1, Nk), 2.0/bb), dlnP_cosmo])
    DlnAB = np.vstack([np.full((1, Nk), 1.0/ba), np.full((1, Nk), 1.0/bb), dlnP_cosmo])

    F_full = np.zeros((Npar, Npar))
    for i in range(Npar):
        for j in range(i, Npar):
            integrand = (
                  0.5 * DlnA[i]  * DlnA[j]  * Raa
                + 0.5 * DlnB[i]  * DlnB[j]  * Rbb
                +       DlnAB[i] * DlnAB[j] * Rxx
                - (DlnAB[i]*DlnA[j]  + DlnA[i] *DlnAB[j]) * Rxa
                - (DlnAB[i]*DlnB[j]  + DlnB[i] *DlnAB[j]) * Rxb
                + 0.5*(DlnA[i]*DlnB[j] + DlnB[i]*DlnA[j]) * Rab
            ) * w
            F_full[i, j] = prefactor * np.sum(integrand)
            F_full[j, i] = F_full[i, j]

    return F_full


# ---------------------------------------------------------------------------
# Survey-level wrappers
# ---------------------------------------------------------------------------

def sigma_mnu_single_tracer_full(
    zarray, nz, bz,
    Area, N_degm2,
    Deltaz=0.2, kmax=0.1,
    Mnu_fid=0.06, dMnu=0.06,
    cosmo=None, Nk=200, return_F=False,
    prior_params=None,
):
    """
    Drop-in replacement for sigma_mnu_single_tracer, marginalising over
    [H0, omega_c, omega_b, ns, As] with optional Planck priors.

    Parameters
    ----------
    prior_params : list of str or None
        Subset of COSMO_PARAMS to apply Planck priors to.
        e.g. ['omega_b', 'omega_c', 'ns', 'As']  (Planck CMB priors)
        or   []   for fully free marginalisation
        or   ['H0', 'omega_b', 'omega_c', 'ns', 'As']  for all but Mnu

    Returns identical outputs as the original function.
    """
    if prior_params is None:
        prior_params = []

    if zarray[0] == 0:
        print('error: z bin cannot start at 0!')
        return

    fid_dict = _make_fid_dict(cosmo, Mnu_fid=Mnu_fid)
    print('Building CAMB derivatives for all cosmological parameters...')
    deriv_interps, Pm_interp, _ = _build_all_derivs(
        cosmo, fid_dict, kmax=max(kmax * 2, 1.0)
    )
    param_order_1t = ['b'] + COSMO_PARAMS
    print(f'Priors applied to: {prior_params if prior_params else "none"}')

    eps      = 1e-4
    dz_array = zarray[1] - zarray[0]
    Nbin     = int((zarray[-1] + dz_array - zarray[0] + eps) // Deltaz)

    list_zbin      = []
    list_sigma_b   = []
    list_sigma_mnu = []
    Flist          = []

    nz     = nz / np.sum(nz)
    size_z = len(zarray)

    for i in range(Nbin):
        imin  = i     * int((size_z + eps) // Nbin)
        imax  = (i+1) * int((size_z + eps) // Nbin)
        nzsum = np.sum(nz[imin:imax])
        if nzsum > 0:
            zbin = zarray[imin] - dz_array / 2 + Deltaz / 2
            Vsur = cosmology.Vsurvey(zbin - Deltaz / 2, Deltaz, Area, cosmo)
            kmin = 2 * math.pi / Vsur**(1.0 / 3)
            bg   = np.sum(nz[imin:imax] * bz[imin:imax]) / nzsum
            n    = nzsum * N_degm2 * Area / Vsur

            # raw full Fisher matrix, no prior, no marginalisation
            F = Mat_Fisher_1tracer_mnu_full(
                    n, bg, zbin, Vsur, kmin, kmax,
                    Mnu_fid=Mnu_fid,
                    deriv_interps=deriv_interps,
                    Pm_interp=Pm_interp,
                    cosmo=cosmo, Nk=Nk)

            # per-bin sigma (no prior, free marginalisation) for diagnostics
            Flist.append(F)
            _apply_priors(F, prior_params, param_order_1t)
            F_inv_bin = np.linalg.inv(F)
            list_zbin.append(zbin)
            list_sigma_b.append(F_inv_bin[0, 0]**0.5)
            list_sigma_mnu.append(F_inv_bin[1, 1]**0.5)

    zeff  = np.sum(zarray * nz)
    Ftot  = np.sum(np.array(Flist), axis=0)

    # --- apply prior ONCE to the total Fisher matrix ---
    _apply_priors(Ftot, prior_params, param_order_1t)

    # marginalise once over cosmo params, keep [b, Mnu]
    Ftot_inv  = np.linalg.inv(Ftot)
    Ftot_marg = np.linalg.inv(Ftot_inv[np.ix_([0, 1], [0, 1])])
    Ftot_marg_inv = np.linalg.inv(Ftot_marg)

    sigma_b_eff   = Ftot_marg_inv[0, 0]**0.5
    sigma_mnu_eff = Ftot_marg_inv[1, 1]**0.5

    if not return_F:
        return list_zbin, list_sigma_b, list_sigma_mnu, zeff, sigma_b_eff, sigma_mnu_eff
    else:
        return list_zbin, list_sigma_b, list_sigma_mnu, zeff, sigma_b_eff, sigma_mnu_eff, Ftot


def sigma_mnu_two_tracers_full(
    zarray, nza, nzb, bza, bzb,
    Area, Na_degm2, Nb_degm2,
    Deltaz=0.2, kmax=0.1,
    Mnu_fid=0.06, dMnu=0.06,
    cosmo=None, Nk=200, return_F=False,
    prior_params=None,
):
    """
    Drop-in replacement for sigma_mnu_two_tracers, marginalising over
    [H0, omega_c, omega_b, ns, As] with optional Planck priors.

    Parameters
    ----------
    prior_params : list of str or None
        Subset of COSMO_PARAMS to apply Planck priors to.
        e.g. ['omega_b', 'omega_c', 'ns', 'As']

    Returns identical outputs as the original function.
    """
    if prior_params is None:
        prior_params = []

    if zarray[0] == 0:
        print('error: z bin cannot start at 0!')
        return

    fid_dict = _make_fid_dict(cosmo, Mnu_fid=Mnu_fid)
    print('Building CAMB derivatives for all cosmological parameters...')
    deriv_interps, Pm_interp, _ = _build_all_derivs(
        cosmo, fid_dict, kmax=max(kmax * 2, 1.0)
    )
    param_order_2t = ['ba', 'bb'] + COSMO_PARAMS
    print(f'Priors applied to: {prior_params if prior_params else "none"}')

    eps      = 1e-4
    dz_array = zarray[1] - zarray[0]
    Nbin     = int((zarray[-1] + dz_array - zarray[0] + eps) // Deltaz)

    list_zbin       = []
    list_sigma_ba   = []
    list_sigma_bb   = []
    list_sigma_mnu  = []
    Flist           = []

    nza    = nza / np.sum(nza)
    nzb    = nzb / np.sum(nzb)
    size_z = len(zarray)

    for i in range(Nbin):
        imin   = i     * int((size_z + eps) // Nbin)
        imax   = (i+1) * int((size_z + eps) // Nbin)
        nzasum = np.sum(nza[imin:imax])
        nzbsum = np.sum(nzb[imin:imax])

        zbin = zarray[imin] - dz_array / 2 + Deltaz / 2
        Vsur = cosmology.Vsurvey(zbin - Deltaz / 2, Deltaz, Area, cosmo)
        kmin = 2 * math.pi / Vsur**(1.0 / 3)

        if nzasum > 0 and nzbsum > 0:
            bga = np.sum(nza[imin:imax] * bza[imin:imax]) / nzasum
            bgb = np.sum(nzb[imin:imax] * bzb[imin:imax]) / nzbsum
            na  = nzasum * Na_degm2 * Area / Vsur
            nb  = nzbsum * Nb_degm2 * Area / Vsur

            # raw full matrix, no prior
            F    = Mat_Fisher_2tracer_mnu_full(
                       na, nb, bga, bgb, zbin, Vsur, kmin, kmax,
                       Mnu_fid=Mnu_fid,
                       deriv_interps=deriv_interps,
                       Pm_interp=Pm_interp,
                       cosmo=cosmo, Nk=Nk)
            Finv = np.linalg.inv(F)
            list_zbin.append(zbin)
            list_sigma_ba.append(Finv[0, 0]**0.5)
            list_sigma_bb.append(Finv[1, 1]**0.5)
            list_sigma_mnu.append(Finv[2, 2]**0.5)
            Flist.append(F)

        elif nzasum > 0:
            bga = np.sum(nza[imin:imax] * bza[imin:imax]) / nzasum
            na  = nzasum * Na_degm2 * Area / Vsur
            F2  = Mat_Fisher_1tracer_mnu_full(
                      na, bga, zbin, Vsur, kmin, kmax,
                      Mnu_fid=Mnu_fid,
                      deriv_interps=deriv_interps,
                      Pm_interp=Pm_interp,
                      cosmo=cosmo, Nk=Nk)
            # embed [ba, Mnu] 2x2 into [ba, bb, Mnu] 3x3, then pad to full Npar
            Npar_2t = 2 + len(COSMO_PARAMS)
            F_full  = np.zeros((Npar_2t, Npar_2t))
            # ba -> index 0, Mnu -> index 2, cosmo -> 3..
            idx_map = [0] + list(range(2, Npar_2t))   # skip bb (index 1)
            for ii, r in enumerate(idx_map):
                for jj, c in enumerate(idx_map):
                    F_full[r, c] = F2[ii, jj]
            Finv2 = np.linalg.inv(F2)
            list_zbin.append(zbin)
            list_sigma_ba.append(Finv2[0, 0]**0.5)
            list_sigma_bb.append(np.inf)
            list_sigma_mnu.append(Finv2[1, 1]**0.5)
            Flist.append(F_full)

        elif nzbsum > 0:
            bgb = np.sum(nzb[imin:imax] * bzb[imin:imax]) / nzbsum
            nb  = nzbsum * Nb_degm2 * Area / Vsur
            F2  = Mat_Fisher_1tracer_mnu_full(
                      nb, bgb, zbin, Vsur, kmin, kmax,
                      Mnu_fid=Mnu_fid,
                      deriv_interps=deriv_interps,
                      Pm_interp=Pm_interp,
                      cosmo=cosmo, Nk=Nk)
            Npar_2t = 2 + len(COSMO_PARAMS)
            F_full  = np.zeros((Npar_2t, Npar_2t))
            idx_map = [1] + list(range(2, Npar_2t))   # skip ba (index 0)
            for ii, r in enumerate(idx_map):
                for jj, c in enumerate(idx_map):
                    F_full[r, c] = F2[ii, jj]
            Finv2 = np.linalg.inv(F2)
            list_zbin.append(zbin)
            list_sigma_ba.append(np.inf)
            list_sigma_bb.append(Finv2[0, 0]**0.5)
            list_sigma_mnu.append(Finv2[1, 1]**0.5)
            Flist.append(F_full)

    zeff = (
        np.sum(zarray * nza) * Na_degm2
        + np.sum(zarray * nzb) * Nb_degm2
    ) / (Na_degm2 + Nb_degm2)

    Ftot = np.sum(np.array(Flist), axis=0)

    # --- apply prior ONCE to the total Fisher matrix ---
    _apply_priors(Ftot, prior_params, param_order_2t)

    # marginalise once, keep [ba, bb, Mnu] (indices 0, 1, 2)
    Ftot_inv  = np.linalg.inv(Ftot)
    Ftot_marg = np.linalg.inv(Ftot_inv[np.ix_([0, 1, 2], [0, 1, 2])])
    Ftot_marg_inv = np.linalg.inv(Ftot_marg)

    sigma_ba_eff  = Ftot_marg_inv[0, 0]**0.5
    sigma_bb_eff  = Ftot_marg_inv[1, 1]**0.5
    sigma_mnu_eff = Ftot_marg_inv[2, 2]**0.5

    if not return_F:
        return (list_zbin,
                list_sigma_ba, list_sigma_bb, list_sigma_mnu,
                zeff, sigma_ba_eff, sigma_bb_eff, sigma_mnu_eff)
    else:
        return (list_zbin,
                list_sigma_ba, list_sigma_bb, list_sigma_mnu,
                zeff, sigma_ba_eff, sigma_bb_eff, sigma_mnu_eff, Ftot)