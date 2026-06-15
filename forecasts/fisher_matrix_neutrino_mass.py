import numpy as np
import math
from scipy.integrate import quad
from scipy.interpolate import interp1d
from functools import partial

import camb
import cosmology

nlim = 10000


# ---------------------------------------------------------------------------
# CAMB-based neutrino derivative  (built once, passed into matrix functions)
# ---------------------------------------------------------------------------

def _build_dlnPm_dSigmamnu(cosmo, Sigma_mnu_fid=0.06, dSigma=0.06,
                             kmin=1e-4, kmax=10.0, npoints=300):
    """
    Compute d ln Pm(k, z=0) / d Sigma_mnu via central finite difference
    using CAMB. Returns a log-space interpolator over k [h/Mpc].
    """
    h     = cosmo['h']
    H0    = 100.0 * h
    ombh2 = cosmo['Omega_b'] * h**2
    omch2 = cosmo['Omega_c'] * h**2
    ns    = cosmo['n_s']
    As    = cosmo['A_s']

    mnu_lo = max(Sigma_mnu_fid - dSigma / 2.0, 0.0)
    mnu_hi = Sigma_mnu_fid + dSigma / 2.0
    step   = mnu_hi - mnu_lo

    def camb_pk(mnu):
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu)
        pars.InitPower.set_params(As=As, ns=ns)
        pars.set_matter_power(redshifts=[0.0], kmax=kmax * 1.1)
        results = camb.get_results(pars)
        kh, _, pk = results.get_matter_power_spectrum(
            minkh=kmin, maxkh=kmax, npoints=npoints
        )
        return kh, pk[0]

    kh, pk_lo = camb_pk(mnu_lo)
    _,  pk_hi = camb_pk(mnu_hi)

    dlnPm = (np.log(pk_hi) - np.log(pk_lo)) / step

    deriv_interp = interp1d(np.log(kh), dlnPm,
                            kind='cubic', fill_value='extrapolate')
    return deriv_interp, kh


def _build_Pm_interp(cosmo, Sigma_mnu_fid=0.06,
                      kmin=1e-4, kmax=10.0, npoints=300):
    """
    Return a log-space interpolator for Pm(k, z=0) [(Mpc/h)^3]
    at the fiducial neutrino mass, from CAMB.
    """
    h     = cosmo['h']
    H0    = 100.0 * h
    ombh2 = cosmo['Omega_b'] * h**2
    omch2 = cosmo['Omega_c'] * h**2
    ns    = cosmo['n_s']
    As    = cosmo['A_s']

    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=Sigma_mnu_fid)
    pars.InitPower.set_params(As=As, ns=ns)
    pars.set_matter_power(redshifts=[0.0], kmax=kmax * 1.1)
    results = camb.get_results(pars)
    kh, _, pk = results.get_matter_power_spectrum(
        minkh=kmin, maxkh=kmax, npoints=npoints
    )
    return interp1d(np.log(kh), np.log(pk[0]),
                    kind='cubic', fill_value='extrapolate')


# ---------------------------------------------------------------------------
# Single-tracer Fisher matrix
# ---------------------------------------------------------------------------
# Parameters: [b,  Sigma_mnu]   (indices 0, 1)
# ---------------------------------------------------------------------------

def Mat_Fisher_1tracer_mnu(n, bg, z, Vsur, kmin, kmax,
                            Sigma_mnu_fid=0.06,
                            dlnPm_interp=None,
                            Pm_interp=None,
                            cosmo=None):
    """
    Return the 2x2 Fisher matrix for a single tracer: [b, Sigma_mnu].

    All inputs are scalars.

    Parameters
    ----------
    n             : number density  [(Mpc/h)^-3]
    bg            : galaxy bias
    z             : mean redshift of the bin
    Vsur          : survey volume  [Mpc^3]
    kmin, kmax    : k range  [h/Mpc]
    Sigma_mnu_fid : fiducial sum of neutrino masses  [eV]
    dlnPm_interp  : interpolator  log k -> d ln Pm / d Sigma_mnu  (from CAMB)
    Pm_interp     : interpolator  log k -> ln Pm(k, z=0)  (from CAMB)
    cosmo         : cosmology dict

    Power spectrum (real space, no mu):
        P(k, z) = bg^2 * Pm_fid(k, z=0) * D(z)^2 / D(0)^2

    Derivatives of ln P:
        d ln P / d b           = 2 / bg
        d ln P / d Sigma_mnu   = d ln Pm(k) / d Sigma_mnu   [k-dependent, from CAMB]
    """
    prefactor = Vsur / (4.0 * math.pi**2)

    Dz = cosmology.D(z, cosmo)
    D0 = cosmology.D(0, cosmo)

    def Pm_fid(k):
        return np.exp(Pm_interp(np.log(k)))

    def P(k):
        return bg**2 * Pm_fid(k) * (Dz / D0)**2

    def dlnP_db():
        return 2.0 / bg                              # k-independent

    def dlnP_dmnu(k):
        return float(dlnPm_interp(np.log(k)))        # k-dependent

    def Fij_integrand(k, i, j):
        pk  = P(k)
        nP  = n * pk
        w   = 0.5 * (nP / (1.0 + nP))**2

        # integrate over mu analytically: int_{-1}^{1} dmu = 2
        # (no mu dependence in real-space P, so the mu integral just gives 2)
        derivs = [dlnP_db(), dlnP_dmnu(k)]
        return prefactor * 2.0 * w * derivs[i] * derivs[j] * k**2

    def build_F():
        F = np.zeros((2, 2))
        for i in range(2):
            for j in range(i, 2):
                val = quad(partial(Fij_integrand, i=i, j=j),
                           kmin, kmax,
                           epsrel=1e-4, epsabs=1e-4, limit=nlim)[0]
                F[i, j] = val
                F[j, i] = val
        return F

    return build_F()


# ---------------------------------------------------------------------------
# Two-tracer Fisher matrix  (vectorised Gauss-Legendre over k only)
# ---------------------------------------------------------------------------
# Parameters: [ba, bb, Sigma_mnu]   (indices 0, 1, 2)
# ---------------------------------------------------------------------------

def Mat_Fisher_2tracer_mnu(
    na, nb, ba, bb, zeff,
    Vsur, kmin, kmax,
    Sigma_mnu_fid=0.06,
    dlnPm_interp=None,
    Pm_interp=None,
    cosmo=None,
    Nk=200
):
    """
    Return the 3x3 Fisher matrix for two tracers: [ba, bb, Sigma_mnu].

    No RSD, no mu integral — real-space power spectra only.
    The mu integral is performed analytically (gives a factor 2).

    Parameters
    ----------
    na, nb        : number densities  [(Mpc/h)^-3]
    ba, bb        : galaxy biases
    zeff          : effective redshift of the bin
    Vsur          : survey volume  [Mpc^3]
    kmin, kmax    : k range  [h/Mpc]
    Sigma_mnu_fid : fiducial Sigma_mnu  [eV]
    dlnPm_interp  : interpolator  log k -> d ln Pm / d Sigma_mnu
    Pm_interp     : interpolator  log k -> ln Pm(k, z=0)
    cosmo         : cosmology dict
    Nk            : number of Gauss-Legendre nodes in k

    Power spectra:
        PA  = ba^2 * Pm_fid(k, z=0) * (D(z)/D(0))^2
        PB  = bb^2 * ...
        PAB = ba*bb * ...

    Derivatives of ln PA  w.r.t. [ba, bb, Sigma_mnu]:
        [2/ba,  0,      d ln Pm / d Sigma_mnu]

    Derivatives of ln PB:
        [0,     2/bb,   d ln Pm / d Sigma_mnu]

    Derivatives of ln PAB:
        [1/ba,  1/bb,   d ln Pm / d Sigma_mnu]
    """
    z  = zeff
    Dz = cosmology.D(z, cosmo)
    D0 = cosmology.D(0, cosmo)

    prefactor = Vsur / (4.0 * np.pi**2)

    # ---------- Gauss-Legendre grid over k ----------
    k_nodes, wk = np.polynomial.legendre.leggauss(Nk)
    k  = 0.5 * (kmax - kmin) * k_nodes + 0.5 * (kmax + kmin)
    wk = 0.5 * (kmax - kmin) * wk

    # ---------- Matter power spectrum from CAMB ----------
    Pm0 = np.exp(Pm_interp(np.log(k))) * (Dz / D0)**2   # shape (Nk,)

    # ---------- k-dependent neutrino derivative ----------
    dlnPm_k = dlnPm_interp(np.log(k))                   # shape (Nk,)

    # ---------- Power spectra ----------
    PA  = ba**2 * Pm0
    PB  = bb**2 * Pm0
    PAB = ba * bb * Pm0

    # ---------- Noise and R coefficients ----------
    nA       = na * PA
    nB       = nb * PB
    PAB2     = PAB**2
    nAnB_X2  = na * nb * PAB2

    one_nA = 1.0 + nA
    one_nB = 1.0 + nB
    det    = one_nA * one_nB - nAnB_X2
    det2   = det**2

    Raa = (nA * one_nB / det)**2
    Rbb = (nB * one_nA / det)**2
    Rxx = na * nb * (one_nA * one_nB + nAnB_X2) * PAB2 / det2
    Rxa = na**2 * nb * one_nB * PAB2 * PA / det2
    Rxb = nb**2 * na * one_nA * PAB2 * PB / det2
    Rab = na**2 * nb**2 * PA * PB * PAB2 / det2

    # ---------- Derivatives of ln P ----------
    # [ba,   bb,   Sigma_mnu]
    zero = np.zeros(Nk)

    DlnA = [
        np.full(Nk, 2.0 / ba),   # d/d ba
        zero,                      # d/d bb
        dlnPm_k,                   # d/d Sigma_mnu
    ]
    DlnB = [
        zero,
        np.full(Nk, 2.0 / bb),
        dlnPm_k,
    ]
    DlnAB = [
        np.full(Nk, 1.0 / ba),
        np.full(Nk, 1.0 / bb),
        dlnPm_k,
    ]

    # ---------- Integration measure (mu integrated analytically => factor 2) ----------
    w = 2.0 * k**2 * wk

    # ---------- Fisher matrix ----------
    F = np.zeros((3, 3))

    for i in range(3):
        for j in range(i, 3):
            integrand = (
                  0.5 * DlnA[i]  * DlnA[j]  * Raa
                + 0.5 * DlnB[i]  * DlnB[j]  * Rbb
                +       DlnAB[i] * DlnAB[j] * Rxx
                - (DlnAB[i] * DlnA[j]  + DlnA[i]  * DlnAB[j]) * Rxa
                - (DlnAB[i] * DlnB[j]  + DlnB[i]  * DlnAB[j]) * Rxb
                + 0.5 * (DlnA[i] * DlnB[j] + DlnB[i] * DlnA[j]) * Rab
            ) * w

            F[i, j] = prefactor * np.sum(integrand)
            F[j, i] = F[i, j]

    return F


# ---------------------------------------------------------------------------
# Survey-level wrappers
# ---------------------------------------------------------------------------

def sigma_mnu_single_tracer(
    zarray, nz, bz,
    Area, N_degm2,
    Deltaz=0.2, kmax=0.1,
    Sigma_mnu_fid=0.06, dSigma=0.06,
    cosmo=None, return_F=False
):
    """
    Return sigma(Sigma_mnu) and sigma(b) per redshift bin and combined,
    for a single tracer (no RSD).

    Fisher parameters: [b,  Sigma_mnu]   (0, 1)

    Returns
    -------
    list_zbin       : bin centres
    list_sigma_b    : sigma(b) per bin
    list_sigma_mnu  : sigma(Sigma_mnu) per bin  [eV]
    zeff            : effective redshift
    sigma_b_eff     : combined sigma(b)
    sigma_mnu_eff   : combined sigma(Sigma_mnu)  [eV]
    [Ftot]          : total 2x2 Fisher matrix  (only if return_F=True)
    """
    if zarray[0] == 0:
        print('error: z bin cannot start at 0!')
        return

    print('Building CAMB interpolators...')
    dlnPm_interp, _ = _build_dlnPm_dSigmamnu(cosmo,
                                               Sigma_mnu_fid=Sigma_mnu_fid,
                                               dSigma=dSigma,
                                               kmax=max(kmax * 2, 1.0))
    Pm_interp = _build_Pm_interp(cosmo,
                                  Sigma_mnu_fid=Sigma_mnu_fid,
                                  kmax=max(kmax * 2, 1.0))
    print('Done.')

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
        imin = i     * int((size_z + eps) // Nbin)
        imax = (i+1) * int((size_z + eps) // Nbin)

        nzsum = np.sum(nz[imin:imax])
        if nzsum > 0:
            zbin = zarray[imin] - dz_array / 2 + Deltaz / 2
            Vsur = cosmology.Vsurvey(zbin - Deltaz / 2, Deltaz, Area, cosmo)
            kmin = 2 * math.pi / Vsur**(1.0 / 3)
            bg   = np.sum(nz[imin:imax] * bz[imin:imax]) / nzsum
            n    = nzsum * N_degm2 * Area / Vsur

            F    = Mat_Fisher_1tracer_mnu(n, bg, zbin, Vsur, kmin, kmax,
                                          Sigma_mnu_fid=Sigma_mnu_fid,
                                          dlnPm_interp=dlnPm_interp,
                                          Pm_interp=Pm_interp,
                                          cosmo=cosmo)
            Finv = np.linalg.inv(F)
            list_zbin.append(zbin)
            list_sigma_b.append(Finv[0, 0]**0.5)
            list_sigma_mnu.append(Finv[1, 1]**0.5)
            Flist.append(F)

    zeff     = np.sum(zarray * nz)
    Ftot     = np.sum(np.array(Flist), axis=0)
    Ftot_inv = np.linalg.inv(Ftot)

    sigma_b_eff   = Ftot_inv[0, 0]**0.5
    sigma_mnu_eff = Ftot_inv[1, 1]**0.5

    if not return_F:
        return list_zbin, list_sigma_b, list_sigma_mnu, zeff, sigma_b_eff, sigma_mnu_eff
    else:
        return list_zbin, list_sigma_b, list_sigma_mnu, zeff, sigma_b_eff, sigma_mnu_eff, Ftot


def sigma_mnu_two_tracers(
    zarray, nza, nzb, bza, bzb,
    Area, Na_degm2, Nb_degm2,
    Deltaz=0.2, kmax=0.1,
    Sigma_mnu_fid=0.06, dSigma=0.06,
    cosmo=None, Nk=200,
    return_F=False
):
    """
    Return sigma(Sigma_mnu), sigma(ba), sigma(bb) per redshift bin and combined,
    for two tracers (no RSD).

    Fisher parameters: [ba,  bb,  Sigma_mnu]   (0, 1, 2)

    Returns
    -------
    list_zbin         : bin centres
    list_sigma_ba     : sigma(ba) per bin
    list_sigma_bb     : sigma(bb) per bin
    list_sigma_mnu    : sigma(Sigma_mnu) per bin  [eV]
    zeff              : effective redshift
    sigma_ba_eff      : combined sigma(ba)
    sigma_bb_eff      : combined sigma(bb)
    sigma_mnu_eff     : combined sigma(Sigma_mnu)  [eV]
    [Ftot]            : total 3x3 Fisher matrix  (only if return_F=True)
    """
    if zarray[0] == 0:
        print('error: z bin cannot start at 0!')
        return

    print('Building CAMB interpolators...')
    dlnPm_interp, _ = _build_dlnPm_dSigmamnu(cosmo,
                                               Sigma_mnu_fid=Sigma_mnu_fid,
                                               dSigma=dSigma,
                                               kmax=max(kmax * 2, 1.0))
    Pm_interp = _build_Pm_interp(cosmo,
                                  Sigma_mnu_fid=Sigma_mnu_fid,
                                  kmax=max(kmax * 2, 1.0))
    print('Done.')

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
        imin = i     * int((size_z + eps) // Nbin)
        imax = (i+1) * int((size_z + eps) // Nbin)

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

            F    = Mat_Fisher_2tracer_mnu(na, nb, bga, bgb, zbin, Vsur, kmin, kmax,
                                          Sigma_mnu_fid=Sigma_mnu_fid,
                                          dlnPm_interp=dlnPm_interp,
                                          Pm_interp=Pm_interp,
                                          cosmo=cosmo, Nk=Nk)
            Finv = np.linalg.inv(F)
            list_zbin.append(zbin)
            list_sigma_ba.append(Finv[0, 0]**0.5)
            list_sigma_bb.append(Finv[1, 1]**0.5)
            list_sigma_mnu.append(Finv[2, 2]**0.5)
            Flist.append(F)

        elif nzasum > 0:   # only tracer A => embed 2x2 in 3x3
            bg  = np.sum(nza[imin:imax] * bza[imin:imax]) / nzasum
            n   = nzasum * Na_degm2 * Area / Vsur
            F2  = Mat_Fisher_1tracer_mnu(n, bg, zbin, Vsur, kmin, kmax,
                                         Sigma_mnu_fid=Sigma_mnu_fid,
                                         dlnPm_interp=dlnPm_interp,
                                         Pm_interp=Pm_interp,
                                         cosmo=cosmo)
            # params in 2x2: [ba, Sigma_mnu] => rows/cols 0, 2 of 3x3
            F3  = np.zeros((3, 3))
            for ii, r in enumerate([0, 2]):
                for jj, c in enumerate([0, 2]):
                    F3[r, c] = F2[ii, jj]
            Finv2 = np.linalg.inv(F2)
            list_zbin.append(zbin)
            list_sigma_ba.append(Finv2[0, 0]**0.5)
            list_sigma_bb.append(np.inf)
            list_sigma_mnu.append(Finv2[1, 1]**0.5)
            Flist.append(F3)

        elif nzbsum > 0:   # only tracer B => embed 2x2 in 3x3
            bg  = np.sum(nzb[imin:imax] * bzb[imin:imax]) / nzbsum
            n   = nzbsum * Nb_degm2 * Area / Vsur
            F2  = Mat_Fisher_1tracer_mnu(n, bg, zbin, Vsur, kmin, kmax,
                                         Sigma_mnu_fid=Sigma_mnu_fid,
                                         dlnPm_interp=dlnPm_interp,
                                         Pm_interp=Pm_interp,
                                         cosmo=cosmo)
            # params in 2x2: [bb, Sigma_mnu] => rows/cols 1, 2 of 3x3
            F3  = np.zeros((3, 3))
            for ii, r in enumerate([1, 2]):
                for jj, c in enumerate([1, 2]):
                    F3[r, c] = F2[ii, jj]
            Finv2 = np.linalg.inv(F2)
            list_zbin.append(zbin)
            list_sigma_ba.append(np.inf)
            list_sigma_bb.append(Finv2[0, 0]**0.5)
            list_sigma_mnu.append(Finv2[1, 1]**0.5)
            Flist.append(F3)

    zeff = (
        np.sum(zarray * nza) * Na_degm2
        + np.sum(zarray * nzb) * Nb_degm2
    ) / (Na_degm2 + Nb_degm2)

    Ftot     = np.sum(np.array(Flist), axis=0)
    Ftot_inv = np.linalg.inv(Ftot)

    sigma_ba_eff  = Ftot_inv[0, 0]**0.5
    sigma_bb_eff  = Ftot_inv[1, 1]**0.5
    sigma_mnu_eff = Ftot_inv[2, 2]**0.5

    if not return_F:
        return (list_zbin,
                list_sigma_ba, list_sigma_bb, list_sigma_mnu,
                zeff,
                sigma_ba_eff, sigma_bb_eff, sigma_mnu_eff)
    else:
        return (list_zbin,
                list_sigma_ba, list_sigma_bb, list_sigma_mnu,
                zeff,
                sigma_ba_eff, sigma_bb_eff, sigma_mnu_eff,
                Ftot)