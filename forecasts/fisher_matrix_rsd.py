import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt

from scipy import integrate
from functools import partial
from scipy.integrate import quad, dblquad

import cosmology

nlim = 10000


def Mat_Fisher_1tracer_rsd(n, bg, z, Vsur, kmin, kmax, cosmo=None):
    """
    Return the 2x2 Fisher matrix for a single tracer, for RSD parameters.
    Parameters are [ln(b*sigma8), ln(f*sigma8)].

    Inputs (all scalars):
        n      : number density in (Mpc/h)^-3
        bg     : galaxy bias b(z)
        z      : mean redshift of the bin
        Vsur   : survey volume in Mpc^3
        kmin   : minimum k of the integral
        kmax   : maximum k of the integral
        cosmo  : cosmology dict

    Following Hamilton (1992) / Blake et al. conventions:
        P(k,z) = (b*s8 + f*s8*mu^2)^2 * Pm(k,z=0) / s8(z=0)^2

    Parameters:
        theta_0 = ln(b*sigma8)
        theta_1 = ln(f*sigma8)

    Derivatives of ln P:
        d ln P / d ln(b*s8) = 2*b*s8 / (b*s8 + f*s8*mu^2)
        d ln P / d ln(f*s8) = 2*f*s8*mu^2 / (b*s8 + f*s8*mu^2)
    """
    prefactor = Vsur / (4 * math.pi**2)

    sig8_z  = cosmology.sigma8_z(z, cosmo)   # sigma8(z)
    sig8_0  = cosmology.sigma8_z(0, cosmo)   # sigma8(z=0)
    f_z     = cosmology.f(z, cosmo)

    bs8 = bg * sig8_z          # b(z)*sigma8(z)
    fs8 = f_z * sig8_z         # f(z)*sigma8(z)

    def P(k, mu):
        """Redshift-space power spectrum in the bs8/fs8 parameterisation."""
        return (bs8 + fs8 * mu**2)**2 * cosmology.Pm(k, 0, cosmo) / sig8_0**2

    def dlnP_dlnbs8(mu):
        """d ln P / d ln(b*sigma8) — independent of k."""
        return 2.0 * bs8 / (bs8 + fs8 * mu**2)

    def dlnP_dlnfs8(mu):
        """d ln P / d ln(f*sigma8) — independent of k."""
        return 2.0 * fs8 * mu**2 / (bs8 + fs8 * mu**2)

    def Fij(k, mu, i, j):
        """
        Integrand of the Fisher matrix element (i,j).
        F_ij = (Vsur/4pi^2) * int dk dmu k^2 * 1/2 * (nP/(1+nP))^2 * d_i * d_j
        where d_i = d ln P / d theta_i
        """
        pk = P(k, mu)
        nP = n * pk
        weight = 0.5 * (nP / (1.0 + nP))**2

        derivs = [dlnP_dlnbs8(mu), dlnP_dlnfs8(mu)]
        return prefactor * weight * derivs[i] * derivs[j] * k**2

    def integrat_k(mu, i, j):
        return quad(partial(Fij, mu=mu, i=i, j=j),
                    kmin, kmax,
                    epsrel=1e-4, epsabs=1e-4, limit=nlim)[0]

    def build_F():
        F = np.zeros((2, 2))
        for i in range(2):
            for j in range(i, 2):
                val = quad(partial(integrat_k, i=i, j=j),
                           -1, 1,
                           epsrel=1e-4, epsabs=1e-4, limit=nlim)[0]
                F[i, j] = val
                F[j, i] = val
        return F

    return build_F()


def Mat_Fisher_2tracer_rsd(
    na, nb, ba, bb, zeff,
    Vsur, kmin, kmax,
    cosmo=None,
    Nk=100,
    Nmu=50
):
    """
    Return the 3x3 Fisher matrix for two tracers, for RSD parameters.
    Parameters are [ln(b_A*sigma8), ln(b_B*sigma8), ln(f*sigma8)],
    i.e. indices 0, 1, 2 respectively.

    Inputs:
        na, nb  : number densities in (Mpc/h)^-3
        ba, bb  : galaxy biases b_A(z), b_B(z)
        zeff    : effective redshift of the bin
        Vsur    : survey volume in Mpc^3
        kmin    : minimum k
        kmax    : maximum k
        cosmo   : cosmology dict
        Nk, Nmu : number of Gauss-Legendre nodes for k and mu

    Power spectra (Kaiser + FoG damping):
        PA  = (bA*s8 + f*s8*mu^2)^2 * Pm0/s8_0^2 * exp(-k^2 mu^2 sig_chi^2)
        PB  = (bB*s8 + f*s8*mu^2)^2 * Pm0/s8_0^2 * exp(-k^2 mu^2 sig_chi^2)
        PAB = (bA*s8 + f*s8*mu^2)*(bB*s8 + f*s8*mu^2) * Pm0/s8_0^2 * exp(...)

    Derivatives of ln P_A w.r.t. [ln(bA s8), ln(bB s8), ln(f s8)]:
        [2*bAs8/(bAs8+fs8 mu^2),  0,  2*fs8 mu^2/(bAs8+fs8 mu^2)]

    Derivatives of ln P_B:
        [0,  2*bBs8/(bBs8+fs8 mu^2),  2*fs8 mu^2/(bBs8+fs8 mu^2)]

    Derivatives of ln P_AB:
        [bAs8/(bAs8+fs8 mu^2),  bBs8/(bBs8+fs8 mu^2),
         fs8 mu^2*(1/(bAs8+fs8 mu^2) + 1/(bBs8+fs8 mu^2))]
    """
    z = zeff

    sig8_z  = cosmology.sigma8_z(z, cosmo)
    sig8_0  = cosmology.sigma8_z(0, cosmo)
    f_z     = cosmology.f(z, cosmo)
    sig2    = cosmology.sigma_chi(z, cosmo)**2

    bAs8 = ba * sig8_z
    bBs8 = bb * sig8_z
    fs8  = f_z * sig8_z

    prefactor = Vsur / (4.0 * np.pi**2)

    # ---------- Gauss-Legendre grids ----------
    k_nodes,  wk  = np.polynomial.legendre.leggauss(Nk)
    mu_nodes, wmu = np.polynomial.legendre.leggauss(Nmu)

    k  = 0.5 * (kmax - kmin) * k_nodes  + 0.5 * (kmax + kmin)
    wk = 0.5 * (kmax - kmin) * wk

    # shapes: (Nmu, Nk)
    k2  = k[None, :]
    mu2 = mu_nodes[:, None]
    wk2 = wk[None, :]
    wmu2 = wmu[:, None]

    # ---------- Base matter power spectrum ----------
    Pm0 = cosmology.Pm(k, 0, cosmo)[None, :]   # shape (1, Nk)

    # ---------- FoG damping ----------
    damp = np.exp(-k2**2 * mu2**2 * sig2)

    # ---------- Effective biases ----------
    Aa = bAs8 + fs8 * mu2**2    # shape (Nmu, Nk)
    Ab = bBs8 + fs8 * mu2**2

    Pbase = Pm0 * damp / sig8_0**2   # common Pm(k,0)/sigma8_0^2 * FoG

    PA  = Aa**2 * Pbase
    PB  = Ab**2 * Pbase
    PAB = Aa * Ab * Pbase

    # ---------- Noise terms ----------
    nA       = na * PA
    nB       = nb * PB
    PAB2     = PAB**2
    nAnB_X2  = na * nb * PAB2

    one_nA = 1.0 + nA
    one_nB = 1.0 + nB

    det  = one_nA * one_nB - nAnB_X2
    det2 = det**2

    # ---------- R coefficients (same as fnl case) ----------
    Raa = (nA * one_nB / det)**2
    Rbb = (nB * one_nA / det)**2
    Rxx = na * nb * (one_nA * one_nB + nAnB_X2) * PAB2 / det2
    Rxa = na**2 * nb * one_nB * PAB2 * PA / det2
    Rxb = nb**2 * na * one_nA * PAB2 * PB / det2
    Rab = na**2 * nb**2 * PA * PB * PAB2 / det2

    # ---------- Derivatives of ln P ----------
    # Parameters: [ln(bA s8), ln(bB s8), ln(f s8)]   =>  indices 0, 1, 2

    zero = np.zeros_like(Pbase)

    # d ln PA / d theta_i
    DlnA = [
        2.0 * bAs8 / Aa,          # d/d ln(bA s8)
        zero,                      # d/d ln(bB s8)
        2.0 * fs8 * mu2**2 / Aa,  # d/d ln(f s8)
    ]

    # d ln PB / d theta_i
    DlnB = [
        zero,
        2.0 * bBs8 / Ab,
        2.0 * fs8 * mu2**2 / Ab,
    ]

    # d ln PAB / d theta_i
    DlnAB = [
        bAs8 / Aa,
        bBs8 / Ab,
        fs8 * mu2**2 * (1.0 / Aa + 1.0 / Ab),
    ]

    # ---------- Integration measure ----------
    w = k2**2 * wk2 * wmu2

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


def sigma_rsd_single_tracer(
    zarray, nz, bz,
    Area, N_degm2,
    Deltaz=0.2, mod='camb', kmax=0.1,
    cosmo=None, return_F=False
):
    """
    Return sigma(b*sigma8) and sigma(f*sigma8) for bins of width Deltaz,
    for a single tracer.

    Returns
    -------
    list_zbin        : list of bin centres
    list_sigma_bs8   : list of sigma(ln b*sigma8) per bin
    list_sigma_fs8   : list of sigma(ln f*sigma8) per bin
    zeff             : effective redshift
    sigma_bs8_eff    : combined sigma(ln b*sigma8)
    sigma_fs8_eff    : combined sigma(ln f*sigma8)
    [Ftot]           : total Fisher matrix (only if return_F=True)
    """
    if zarray[0] == 0:
        print('error: z bin cannot start at 0!')
        return

    eps = 1e-4
    dz_array = zarray[1] - zarray[0]
    Nbin = int((zarray[-1] + dz_array - zarray[0] + eps) // Deltaz)

    list_zbin      = []
    list_sigma_bs8 = []
    list_sigma_fs8 = []
    Flist          = []

    nz = nz / np.sum(nz)
    size_z = len(zarray)

    V = 0

    for i in range(Nbin):
        imin = i     * int((size_z + eps) // Nbin)
        imax = (i+1) * int((size_z + eps) // Nbin)

        nzsum = np.sum(nz[imin:imax])
        if nzsum > 0:
            zbin  = zarray[imin] - dz_array / 2 + Deltaz / 2
            Vsur  = cosmology.Vsurvey(zbin - Deltaz / 2, Deltaz, Area, cosmo)
            V += Vsur
            kmin  = 2 * math.pi / Vsur**(1.0 / 3)
            bg    = np.sum(nz[imin:imax] * bz[imin:imax]) / nzsum
            n     = nzsum * N_degm2 * Area / Vsur

            F = Mat_Fisher_1tracer_rsd(n, bg, zbin, Vsur, kmin, kmax, cosmo=cosmo)
            Finv = np.linalg.inv(F)

            list_zbin.append(zbin)
            list_sigma_bs8.append(Finv[0, 0]**0.5)
            list_sigma_fs8.append(Finv[1, 1]**0.5)
            Flist.append(F)

    print(V)
    zeff = np.sum(zarray * nz)
    Ftot = np.sum(np.array(Flist), axis=0)
    Ftot_inv = np.linalg.inv(Ftot)
    sigma_bs8_eff = Ftot_inv[0, 0]**0.5
    sigma_fs8_eff = Ftot_inv[1, 1]**0.5

    if not return_F:
        return list_zbin, list_sigma_bs8, list_sigma_fs8, zeff, sigma_bs8_eff, sigma_fs8_eff
    else:
        return list_zbin, list_sigma_bs8, list_sigma_fs8, zeff, sigma_bs8_eff, sigma_fs8_eff, Ftot


def sigma_rsd_two_tracers(
    zarray, nza, nzb, bza, bzb,
    Area, Na_degm2, Nb_degm2,
    Deltaz=0.2, kmax=0.1,
    cosmo=None, Nk=100, Nmu=50,
    return_F=False
):
    """
    Return sigma(bA*sigma8), sigma(bB*sigma8), sigma(f*sigma8) for bins of
    width Deltaz, for two tracers.

    The Fisher matrix has parameters [ln(bA s8), ln(bB s8), ln(f s8)].

    Returns
    -------
    list_zbin         : list of bin centres
    list_sigma_bAs8   : list of sigma(ln bA*sigma8) per bin
    list_sigma_bBs8   : list of sigma(ln bB*sigma8) per bin
    list_sigma_fs8    : list of sigma(ln f*sigma8) per bin
    zeff              : effective redshift
    sigma_bAs8_eff    : combined sigma(ln bA*sigma8)
    sigma_bBs8_eff    : combined sigma(ln bB*sigma8)
    sigma_fs8_eff     : combined sigma(ln f*sigma8)
    [Ftot]            : total Fisher matrix (only if return_F=True)
    """
    if zarray[0] == 0:
        print('error: z bin cannot start at 0!')
        return

    eps = 1e-4
    dz_array = zarray[1] - zarray[0]
    Nbin = int((zarray[-1] + dz_array - zarray[0] + eps) // Deltaz)

    list_zbin       = []
    list_sigma_bAs8 = []
    list_sigma_bBs8 = []
    list_sigma_fs8  = []
    Flist           = []

    nza = nza / np.sum(nza)
    nzb = nzb / np.sum(nzb)
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

            F    = Mat_Fisher_2tracer_rsd(na, nb, bga, bgb, zbin, Vsur, kmin, kmax,
                                          cosmo=cosmo, Nk=Nk, Nmu=Nmu)
            Finv = np.linalg.inv(F)
            list_zbin.append(zbin)
            list_sigma_bAs8.append(Finv[0, 0]**0.5)
            list_sigma_bBs8.append(Finv[1, 1]**0.5)
            list_sigma_fs8.append(Finv[2, 2]**0.5)
            Flist.append(F)

        elif nzasum > 0:   # only tracer A
            bg = np.sum(nza[imin:imax] * bza[imin:imax]) / nzasum
            n  = nzasum * Na_degm2 * Area / Vsur
            F1d = Mat_Fisher_1tracer_rsd(n, bg, zbin, Vsur, kmin, kmax, cosmo=cosmo)
            # embed in 3x3: params [bA s8, bB s8, f s8] => 1D gives [bA s8, f s8]
            F3  = np.zeros((3, 3))
            F3[0, 0] = F1d[0, 0]
            F3[0, 2] = F1d[0, 1]
            F3[2, 0] = F1d[1, 0]
            F3[2, 2] = F1d[1, 1]
            Finv = np.linalg.inv(F1d)
            list_zbin.append(zbin)
            list_sigma_bAs8.append(Finv[0, 0]**0.5)
            list_sigma_bBs8.append(np.inf)
            list_sigma_fs8.append(Finv[1, 1]**0.5)
            Flist.append(F3)

        elif nzbsum > 0:   # only tracer B
            bg = np.sum(nzb[imin:imax] * bzb[imin:imax]) / nzbsum
            n  = nzbsum * Nb_degm2 * Area / Vsur
            F1d = Mat_Fisher_1tracer_rsd(n, bg, zbin, Vsur, kmin, kmax, cosmo=cosmo)
            # embed in 3x3: params [bA s8, bB s8, f s8] => 1D gives [bB s8, f s8]
            F3  = np.zeros((3, 3))
            F3[1, 1] = F1d[0, 0]
            F3[1, 2] = F1d[0, 1]
            F3[2, 1] = F1d[1, 0]
            F3[2, 2] = F1d[1, 1]
            Finv = np.linalg.inv(F1d)
            list_zbin.append(zbin)
            list_sigma_bAs8.append(np.inf)
            list_sigma_bBs8.append(Finv[0, 0]**0.5)
            list_sigma_fs8.append(Finv[1, 1]**0.5)
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
