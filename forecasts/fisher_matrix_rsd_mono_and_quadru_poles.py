import numpy as np
import math
from scipy.integrate import quad
from functools import partial

import cosmology

nlim = 10000

# Legendre polynomials
def _L0(mu): return np.ones_like(mu)
def _L2(mu): return 0.5 * (3.0 * mu**2 - 1.0)


# ---------------------------------------------------------------------------
# Multipole covariance building block
# ---------------------------------------------------------------------------

def _build_Cov_multipoles(n, P_mu_func, Nmu=200):
    """
    Compute the 2x2 covariance matrix of [P0, P2] at a given k.

    The CORRECT Gaussian covariance of the multipole estimators is:

        C_{ll'} = (2l+1)(2l'+1)/2 * int_{-1}^{1} dmu/2 * L_l(mu) * L_l'(mu)
                  * [P(k,mu) + 1/n]^2

    where P(k,mu) + 1/n is the total observed power (signal + shot noise).
    This is the Percival & White (2009) form.

    Note: do NOT use [P/(1+nP)]^2 here — that is a different (incorrect)
    normalisation that does not reduce to the standard mu-integral Fisher.
    """
    mu_nodes, wmu = np.polynomial.legendre.leggauss(Nmu)

    Pmu   = P_mu_func(mu_nodes)          # signal power  P(k, mu)
    Ptot  = Pmu + 1.0 / n               # total power   P + 1/n

    L       = [_L0(mu_nodes), _L2(mu_nodes)]
    prefacs = [1.0, 5.0]                 # (2l+1) for l=0,2

    C = np.zeros((2, 2))
    for a in range(2):
        for b in range(a, 2):
            C[a, b] = (prefacs[a] * prefacs[b] / 2.0
                       * np.sum(L[a] * L[b] * Ptot**2 * wmu / 2.0))
            C[b, a] = C[a, b]
    return C


# ---------------------------------------------------------------------------
# Multipole power spectra and their derivatives
# ---------------------------------------------------------------------------

def _multipoles_and_derivs(bs8, fs8, Pm_k, sig8_0):
    """
    Compute P0, P2 and their derivatives w.r.t. [ln(b s8), ln(f s8)].

    These are ABSOLUTE derivatives dP_l/d theta (not d ln P_l / d theta),
    because the Fisher formula uses dP @ C^{-1} @ dP, not d ln P.
    """
    base = Pm_k / sig8_0**2

    P0 = (bs8**2 + 2.0/3.0 * bs8 * fs8 + 1.0/5.0 * fs8**2) * base
    P2 = (4.0/3.0 * bs8 * fs8 + 4.0/7.0 * fs8**2) * base

    # d/d ln(bs8) = bs8 * d/d bs8
    dP0_dlnbs8 = (2.0 * bs8**2 + 2.0/3.0 * bs8 * fs8) * base
    dP2_dlnbs8 = (4.0/3.0 * bs8 * fs8) * base

    # d/d ln(fs8) = fs8 * d/d fs8
    dP0_dlnfs8 = (2.0/3.0 * bs8 * fs8 + 2.0/5.0 * fs8**2) * base
    dP2_dlnfs8 = (4.0/3.0 * bs8 * fs8 + 8.0/7.0 * fs8**2) * base

    dP0 = [dP0_dlnbs8, dP0_dlnfs8]
    dP2 = [dP2_dlnbs8, dP2_dlnfs8]

    return P0, P2, dP0, dP2


# ---------------------------------------------------------------------------
# Single-tracer Fisher matrix
# ---------------------------------------------------------------------------
# Parameters: [ln(b*sigma8),  ln(f*sigma8)]   (indices 0, 1)
#
# Equivalence with the standard mu-integral Fisher (fisher_rsd):
# ---------------------------------------------------------------
# The two approaches are equivalent when ALL multipoles are used jointly,
# provided the covariance uses (P + 1/n)^2 and the Fisher has a 1/2 prefactor:
#
#   F_ij = 1/2 * sum_{ll'} dPl/dti * C^{-1}_{ll'} * dPl'/dtj
#        = int_{-1}^1 dmu * 1/2 * (nP/(1+nP))^2 * dlnP/dti * dlnP/dtj
#
# This was verified numerically; the ratio between both methods is ~1.
#
# Degeneracy note:
# ----------------
# With only the monopole or only the quadrupole, b*s8 and f*s8 are degenerate
# (rank-1 Fisher matrix). Use fix_bias=True to get conditional (unmarginalized)
# constraints assuming the other parameter is known.
# ---------------------------------------------------------------------------

def Mat_Fisher_1tracer_multipoles(n, bg, z, Vsur, kmin, kmax,
                                   cosmo=None, Nmu=200,
                                   multipoles='both'):
    """
    2x2 Fisher matrix for a single tracer.
    Parameters: [ln(b*sigma8),  ln(f*sigma8)].

        F_ij = Vsur/(4pi^2) * int dk k^2
               * 1/2 * sum_{ll'} (dP_l/dtheta_i) [C^{-1}]_{ll'} (dP_{l'}/dtheta_j)

    where C_{ll'} = (2l+1)(2l'+1)/2 * int dmu/2 * L_l L_l' * (P+1/n)^2.
    The factor 1/2 ensures equivalence with the standard mu-integral Fisher.

    multipoles : 'monopole'   — use P0 only  (Fisher is rank-1)
                 'quadrupole' — use P2 only  (Fisher is rank-1)
                 'both'       — use P0 + P2 jointly  (full rank, default)
                                Equivalent to the standard mu-integral Fisher.
    """
    if multipoles not in ('monopole', 'quadrupole', 'both'):
        raise ValueError("multipoles must be 'monopole', 'quadrupole', or 'both'")

    if multipoles == 'monopole':
        idx = [0]
    elif multipoles == 'quadrupole':
        idx = [1]
    else:
        idx = [0, 1]

    prefactor = Vsur / (4.0 * math.pi**2)

    sig8_z = cosmology.sigma8_z(z, cosmo)
    sig8_0 = cosmology.sigma8_z(0, cosmo)
    f_z    = cosmology.f(z, cosmo)

    bs8 = bg  * sig8_z
    fs8 = f_z * sig8_z

    def P_mu(k, mu):
        return (bs8 + fs8 * mu**2)**2 * cosmology.Pm(k, 0, cosmo) / sig8_0**2

    def integrand_k(k):
        Pm_k = cosmology.Pm(k, 0, cosmo)
        P0, P2, dP0, dP2 = _multipoles_and_derivs(bs8, fs8, Pm_k, sig8_0)

        # Build full 2x2 covariance, then slice to selected multipoles
        C_full = _build_Cov_multipoles(n, lambda mu: P_mu(k, mu), Nmu=Nmu)
        C      = C_full[np.ix_(idx, idx)]
        Cinv   = np.linalg.inv(C)

        # Derivative vectors (absolute dP, not d ln P), sliced
        dPall = [
            np.array([dP0[0], dP2[0]]),   # d/d ln(bs8)
            np.array([dP0[1], dP2[1]]),   # d/d ln(fs8)
        ]
        dPvec = [d[idx] for d in dPall]

        # Factor 1/2: reality condition P(k) = P(-k)
        F_k = np.zeros((2, 2))
        for i in range(2):
            for j in range(i, 2):
                F_k[i, j] = 0.5 * dPvec[i] @ Cinv @ dPvec[j]
                F_k[j, i] = F_k[i, j]

        return F_k * k**2

    def build_F():
        F = np.zeros((2, 2))
        for i in range(2):
            for j in range(i, 2):
                def scalar_integrand(k, _i=i, _j=j):
                    return integrand_k(k)[_i, _j]
                val = quad(scalar_integrand, kmin, kmax,
                           epsrel=1e-4, epsabs=1e-4, limit=nlim)[0]
                F[i, j] = prefactor * val
                F[j, i] = F[i, j]
        return F

    return build_F()


# ---------------------------------------------------------------------------
# Two-tracer Fisher matrix  (vectorised Gauss-Legendre)
# ---------------------------------------------------------------------------

def Mat_Fisher_2tracer_multipoles(
    na, nb, ba, bb, zeff,
    Vsur, kmin, kmax,
    cosmo=None,
    Nk=150, Nmu=100
):
    """
    3x3 Fisher matrix for two tracers using monopole + quadrupole jointly.
    Parameters: [ln(bA*sigma8),  ln(bB*sigma8),  ln(f*sigma8)].
    Data vector: [P0_A, P2_A, P0_B, P2_B].

    Covariance uses (P + 1/n)^2 for each auto-spectrum and
    P_AB^2 / sqrt((P_A + 1/nA)(P_B + 1/nB)) for cross terms (Gaussian approx).
    The extra 1/2 factor is included for equivalence with the mu-integral Fisher.
    """
    z = zeff

    sig8_z = cosmology.sigma8_z(z, cosmo)
    sig8_0 = cosmology.sigma8_z(0, cosmo)
    f_z    = cosmology.f(z, cosmo)
    sig2   = cosmology.sigma_chi(z, cosmo)**2

    bAs8 = ba * sig8_z
    bBs8 = bb * sig8_z
    fs8  = f_z * sig8_z

    prefactor = Vsur / (4.0 * np.pi**2)

    k_nodes, wk   = np.polynomial.legendre.leggauss(Nk)
    mu_nodes, wmu = np.polynomial.legendre.leggauss(Nmu)

    k  = 0.5 * (kmax - kmin) * k_nodes + 0.5 * (kmax + kmin)
    wk = 0.5 * (kmax - kmin) * wk

    base = cosmology.Pm(k, 0, cosmo) / sig8_0**2   # (Nk,)

    P0A, P2A, dP0A, dP2A = _multipoles_and_derivs(bAs8, fs8, base * sig8_0**2, sig8_0)
    P0B, P2B, dP0B, dP2B = _multipoles_and_derivs(bBs8, fs8, base * sig8_0**2, sig8_0)

    mu2 = mu_nodes          # (Nmu,)
    k_c = k[:, None]        # (Nk, 1)

    damp    = np.exp(-k_c**2 * mu2**2 * sig2)
    PA_mu   = (bAs8 + fs8 * mu2**2)**2 * base[:, None] * damp   # (Nk, Nmu)
    PB_mu   = (bBs8 + fs8 * mu2**2)**2 * base[:, None] * damp
    PAB_mu  = (bAs8 + fs8 * mu2**2) * (bBs8 + fs8 * mu2**2) * base[:, None] * damp

    # CORRECT: use (P + 1/n)^2 for auto, P_AB^2 for cross
    PtotA   = PA_mu  + 1.0 / na    # (Nk, Nmu)
    PtotB   = PB_mu  + 1.0 / nb

    L0mu    = _L0(mu2)
    L2mu    = _L2(mu2)
    Lmat    = [L0mu, L2mu]
    prefacs = [1.0, 5.0]

    CovAA = np.zeros((Nk, 2, 2))
    CovBB = np.zeros((Nk, 2, 2))
    CovAB = np.zeros((Nk, 2, 2))

    for a in range(2):
        for b in range(a, 2):
            pref = prefacs[a] * prefacs[b] / 2.0
            wL   = wmu / 2.0 * Lmat[a] * Lmat[b]   # (Nmu,)
            CovAA[:, a, b] = pref * (PtotA**2  * wL).sum(axis=1)
            CovBB[:, a, b] = pref * (PtotB**2  * wL).sum(axis=1)
            # Cross-covariance: Cov(PlA, PlB) ~ int Ll Ll' PAB^2 dmu/2
            # (off-diagonal tracer block: signal only, no shot noise cross term)
            CovAB[:, a, b] = pref * (PAB_mu**2 * wL).sum(axis=1)
            CovAA[:, b, a] = CovAA[:, a, b]
            CovBB[:, b, a] = CovBB[:, b, a]
            CovAB[:, b, a] = CovAB[:, a, b]

    # Full 4x4 covariance: [P0A, P2A, P0B, P2B]
    Cov4 = np.zeros((Nk, 4, 4))
    Cov4[:, 0:2, 0:2] = CovAA
    Cov4[:, 2:4, 2:4] = CovBB
    Cov4[:, 0:2, 2:4] = CovAB
    Cov4[:, 2:4, 0:2] = CovAB

    Cinv4 = np.linalg.inv(Cov4)

    zero = np.zeros(Nk)
    dPvec = [
        np.stack([dP0A[0], dP2A[0], zero,    zero   ], axis=1),  # d/d ln(bA s8)
        np.stack([zero,    zero,    dP0B[0], dP2B[0]], axis=1),  # d/d ln(bB s8)
        np.stack([dP0A[1], dP2A[1], dP0B[1], dP2B[1]], axis=1), # d/d ln(f s8)
    ]

    # Factor 1/2: reality condition
    F = np.zeros((3, 3))
    for i in range(3):
        for j in range(i, 3):
            Fk = 0.5 * np.einsum('ki,kij,kj->k', dPvec[i], Cinv4, dPvec[j])
            F[i, j] = prefactor * np.sum(Fk * k**2 * wk)
            F[j, i] = F[i, j]

    return F


# ---------------------------------------------------------------------------
# Helper: extract sigmas safely from a Fisher matrix
# ---------------------------------------------------------------------------

def _sigmas_from_Fisher(F, fix_bias=False):
    """
    Return (sigma_bs8, sigma_fs8) from a 2x2 Fisher matrix.

    fix_bias=False (default, multipoles='both'):
        Invert the full 2x2 matrix — marginalized constraints.

    fix_bias=True (for single-multipole cases):
        The Fisher matrix is rank-1 and cannot be inverted.
        Return conditional (unmarginalized) constraints:
            sigma = 1/sqrt(F[i,i])   (other parameter fixed)
        This is the best-case constraint assuming the other param is known.
        A zero diagonal means the parameter is unconstrained (returns inf).
    """
    if fix_bias:
        sigma_bs8 = 1.0 / np.sqrt(F[0, 0]) if F[0, 0] > 0 else np.inf
        sigma_fs8 = 1.0 / np.sqrt(F[1, 1]) if F[1, 1] > 0 else np.inf
        return sigma_bs8, sigma_fs8
    else:
        Finv = np.linalg.inv(F)
        return Finv[0, 0]**0.5, Finv[1, 1]**0.5


# ---------------------------------------------------------------------------
# Survey-level wrappers
# ---------------------------------------------------------------------------

def sigma_multipoles_single_tracer(
    zarray, nz, bz,
    Area, N_degm2,
    Deltaz=0.2, kmax=0.1,
    cosmo=None, Nmu=200,
    multipoles='both',
    fix_bias=None,
    return_F=False
):
    """
    Return sigma(ln b*sigma8) and sigma(ln f*sigma8) per redshift bin
    and combined, for a single tracer.

    Parameters: [ln(b s8),  ln(f s8)]   (0, 1)

    Parameters
    ----------
    multipoles : 'monopole'   — use P0 only
                 'quadrupole' — use P2 only
                 'both'       — use P0 + P2 jointly (default)
                                Equivalent to the standard mu-integral Fisher.

    fix_bias   : bool or None
        None (default): automatically True when multipoles != 'both'.
        True  : conditional (unmarginalized) sigmas — other param fixed.
        False : invert the full 2x2 Fisher (only valid for multipoles='both').

    Returns
    -------
    list_zbin        : bin centres
    list_sigma_bs8   : sigma(ln b*sigma8) per bin
    list_sigma_fs8   : sigma(ln f*sigma8) per bin
    zeff             : effective redshift
    sigma_bs8_eff    : combined sigma(ln b*sigma8)
    sigma_fs8_eff    : combined sigma(ln f*sigma8)
    [Ftot]           : total 2x2 Fisher matrix (only if return_F=True)
    """
    if zarray[0] == 0:
        print('error: z bin cannot start at 0!')
        return

    if fix_bias is None:
        fix_bias = (multipoles != 'both')
    if fix_bias is False and multipoles != 'both':
        print("Warning: fix_bias=False with a single multipole gives a singular "
              "Fisher matrix. Setting fix_bias=True.")
        fix_bias = True

    eps      = 1e-4
    dz_array = zarray[1] - zarray[0]
    Nbin     = int((zarray[-1] + dz_array - zarray[0] + eps) // Deltaz)

    list_zbin      = []
    list_sigma_bs8 = []
    list_sigma_fs8 = []
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

            F = Mat_Fisher_1tracer_multipoles(n, bg, zbin, Vsur, kmin, kmax,
                                              cosmo=cosmo, Nmu=Nmu,
                                              multipoles=multipoles)
            sigma_bs8, sigma_fs8 = _sigmas_from_Fisher(F, fix_bias=fix_bias)

            list_zbin.append(zbin)
            list_sigma_bs8.append(sigma_bs8)
            list_sigma_fs8.append(sigma_fs8)
            Flist.append(F)

    zeff  = np.sum(zarray * nz)
    Ftot  = np.sum(np.array(Flist), axis=0)
    sigma_bs8_eff, sigma_fs8_eff = _sigmas_from_Fisher(Ftot, fix_bias=fix_bias)

    if not return_F:
        return list_zbin, list_sigma_bs8, list_sigma_fs8, zeff, sigma_bs8_eff, sigma_fs8_eff
    else:
        return list_zbin, list_sigma_bs8, list_sigma_fs8, zeff, sigma_bs8_eff, sigma_fs8_eff, Ftot


def sigma_multipoles_two_tracers(
    zarray, nza, nzb, bza, bzb,
    Area, Na_degm2, Nb_degm2,
    Deltaz=0.2, kmax=0.1,
    cosmo=None, Nk=150, Nmu=100,
    return_F=False
):
    """
    Return sigma(ln bA*s8), sigma(ln bB*s8), sigma(ln f*s8) per redshift bin
    and combined, for two tracers using monopole + quadrupole jointly.

    Parameters: [ln(bA s8),  ln(bB s8),  ln(f s8)]   (0, 1, 2)
    """
    if zarray[0] == 0:
        print('error: z bin cannot start at 0!')
        return

    eps      = 1e-4
    dz_array = zarray[1] - zarray[0]
    Nbin     = int((zarray[-1] + dz_array - zarray[0] + eps) // Deltaz)

    list_zbin       = []
    list_sigma_bAs8 = []
    list_sigma_bBs8 = []
    list_sigma_fs8  = []
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

            F    = Mat_Fisher_2tracer_multipoles(na, nb, bga, bgb, zbin, Vsur, kmin, kmax,
                                                 cosmo=cosmo, Nk=Nk, Nmu=Nmu)
            Finv = np.linalg.inv(F)
            list_zbin.append(zbin)
            list_sigma_bAs8.append(Finv[0, 0]**0.5)
            list_sigma_bBs8.append(Finv[1, 1]**0.5)
            list_sigma_fs8.append(Finv[2, 2]**0.5)
            Flist.append(F)

        elif nzasum > 0:
            bg  = np.sum(nza[imin:imax] * bza[imin:imax]) / nzasum
            n   = nzasum * Na_degm2 * Area / Vsur
            F2  = Mat_Fisher_1tracer_multipoles(n, bg, zbin, Vsur, kmin, kmax,
                                                cosmo=cosmo, Nmu=Nmu)
            F3  = np.zeros((3, 3))
            F3[0, 0] = F2[0, 0]; F3[0, 2] = F2[0, 1]
            F3[2, 0] = F2[1, 0]; F3[2, 2] = F2[1, 1]
            Finv2 = np.linalg.inv(F2)
            list_zbin.append(zbin)
            list_sigma_bAs8.append(Finv2[0, 0]**0.5)
            list_sigma_bBs8.append(np.inf)
            list_sigma_fs8.append(Finv2[1, 1]**0.5)
            Flist.append(F3)

        elif nzbsum > 0:
            bg  = np.sum(nzb[imin:imax] * bzb[imin:imax]) / nzbsum
            n   = nzbsum * Nb_degm2 * Area / Vsur
            F2  = Mat_Fisher_1tracer_multipoles(n, bg, zbin, Vsur, kmin, kmax,
                                                cosmo=cosmo, Nmu=Nmu)
            F3  = np.zeros((3, 3))
            F3[1, 1] = F2[0, 0]; F3[1, 2] = F2[0, 1]
            F3[2, 1] = F2[1, 0]; F3[2, 2] = F2[1, 1]
            Finv2 = np.linalg.inv(F2)
            list_zbin.append(zbin)
            list_sigma_bAs8.append(np.inf)
            list_sigma_bBs8.append(Finv2[0, 0]**0.5)
            list_sigma_fs8.append(Finv2[1, 1]**0.5)
            Flist.append(F3)

    zeff = (
        np.sum(zarray * nza) * Na_degm2
        + np.sum(zarray * nzb) * Nb_degm2
    ) / (Na_degm2 + Nb_degm2)

    Ftot     = np.sum(np.array(Flist), axis=0)
    Ftot_inv = np.linalg.inv(Ftot)
    sigma_bAs8_eff = Ftot_inv[0, 0]**0.5
    sigma_bBs8_eff = Ftot_inv[1, 1]**0.5
    sigma_fs8_eff  = Ftot_inv[2, 2]**0.5

    if not return_F:
        return (list_zbin,
                list_sigma_bAs8, list_sigma_bBs8, list_sigma_fs8,
                zeff,
                sigma_bAs8_eff, sigma_bBs8_eff, sigma_fs8_eff)
    else:
        return (list_zbin,
                list_sigma_bAs8, list_sigma_bBs8, list_sigma_fs8,
                zeff,
                sigma_bAs8_eff, sigma_bBs8_eff, sigma_fs8_eff,
                Ftot)