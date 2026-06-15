import numpy as np
import math
import matplotlib.pyplot as plt

import cosmology


def compute_nbP(
    zarray, nz, bz,
    Area, N_degm2,
    k=0.1,
    Deltaz=0.2,
    cosmo=None
):
    """
    Compute n(z) * b(z)^2 * P(k, z) per redshift bin, for a single tracer.

    This is the key signal-to-noise quantity: nP >> 1 means sample-variance
    limited, nP << 1 means shot-noise limited.

    Parameters
    ----------
    zarray   : array of redshift bin centres (must not start at 0)
    nz       : array of redshift distribution (will be normalised internally)
    bz       : array of galaxy bias b(z)
    Area     : survey area  [deg^2]
    N_degm2  : number density of galaxies  [deg^-2]
    k        : wavenumber at which to evaluate P(k, z)  [h/Mpc]  (scalar or array)
    Deltaz   : redshift bin width
    cosmo    : cosmology dict

    Returns
    -------
    list_zbin : list of bin centre redshifts
    list_n    : effective number density n  [(Mpc/h)^-3]  per bin
    list_b    : effective bias per bin
    list_Pm   : matter power spectrum Pm(k, z)  [(Mpc/h)^3]  per bin
    list_nP   : n * b^2 * Pm(k, z)  [dimensionless]  per bin
    """
    if zarray[0] == 0:
        print('error: z bin cannot start at 0!')
        return

    eps      = 1e-4
    dz_array = zarray[1] - zarray[0]
    Nbin     = int((zarray[-1] + dz_array - zarray[0] + eps) // Deltaz)

    nz = nz / np.sum(nz)
    size_z = len(zarray)

    list_zbin = []
    list_n    = []
    list_b    = []
    list_Pm   = []
    list_nP   = []

    list_Vsur = []
    for i in range(Nbin):
        imin = i     * int((size_z + eps) // Nbin)
        imax = (i+1) * int((size_z + eps) // Nbin)

        nzsum = np.sum(nz[imin:imax])
        if nzsum > 0:
            zbin = zarray[imin] - dz_array / 2 + Deltaz / 2
            Vsur = cosmology.Vsurvey(zbin - Deltaz / 2, Deltaz, Area, cosmo)
            bg   = np.sum(nz[imin:imax] * bz[imin:imax]) / nzsum
            n    = nzsum * N_degm2 * Area / Vsur
            Pm   = cosmology.Pm(k, zbin, cosmo)   # scalar or array depending on k

            list_Vsur.append(Vsur)
            list_zbin.append(zbin)
            list_n.append(n)
            list_b.append(bg)
            list_Pm.append(Pm)
            list_nP.append(n * bg**2 * Pm)

    return list_zbin, list_n, list_b, list_Pm, list_nP, np.mean(np.array(list_nP))


def compute_nbP_two_tracers(
    zarray, nza, nzb, bza, bzb,
    Area, Na_degm2, Nb_degm2,
    k=0.1,
    Deltaz=0.2,
    cosmo=None
):
    """
    Compute n(z) * b(z)^2 * P(k, z) per redshift bin, for two tracers A and B.

    Parameters
    ----------
    zarray              : array of redshift bin centres (must not start at 0)
    nza, nzb            : redshift distributions for tracers A and B
    bza, bzb            : galaxy bias arrays for A and B
    Area                : survey area  [deg^2]
    Na_degm2, Nb_degm2  : number densities  [deg^-2]
    k                   : wavenumber  [h/Mpc]  (scalar or array)
    Deltaz              : redshift bin width
    cosmo               : cosmology dict

    Returns
    -------
    list_zbin  : bin centre redshifts
    list_nA    : effective number density of A  [(Mpc/h)^-3]
    list_nB    : effective number density of B  [(Mpc/h)^-3]
    list_bA    : effective bias of A
    list_bB    : effective bias of B
    list_Pm    : matter power spectrum Pm(k, z)
    list_nPA   : nA * bA^2 * Pm  per bin
    list_nPB   : nB * bB^2 * Pm  per bin
    """
    if zarray[0] == 0:
        print('error: z bin cannot start at 0!')
        return

    eps      = 1e-4
    dz_array = zarray[1] - zarray[0]
    Nbin     = int((zarray[-1] + dz_array - zarray[0] + eps) // Deltaz)

    nza = nza / np.sum(nza)
    nzb = nzb / np.sum(nzb)
    size_z = len(zarray)

    list_zbin = []
    list_nA   = []
    list_nB   = []
    list_bA   = []
    list_bB   = []
    list_Pm   = []
    list_nPA  = []
    list_nPB  = []

    for i in range(Nbin):
        imin = i     * int((size_z + eps) // Nbin)
        imax = (i+1) * int((size_z + eps) // Nbin)

        nzasum = np.sum(nza[imin:imax])
        nzbsum = np.sum(nzb[imin:imax])

        zbin = zarray[imin] - dz_array / 2 + Deltaz / 2
        Vsur = cosmology.Vsurvey(zbin - Deltaz / 2, Deltaz, Area, cosmo)
        Pm   = cosmology.Pm(k, zbin, cosmo)

        bga = np.sum(nza[imin:imax] * bza[imin:imax]) / nzasum if nzasum > 0 else 0.0
        bgb = np.sum(nzb[imin:imax] * bzb[imin:imax]) / nzbsum if nzbsum > 0 else 0.0
        na  = nzasum * Na_degm2 * Area / Vsur
        nb  = nzbsum * Nb_degm2 * Area / Vsur

        list_zbin.append(zbin)
        list_nA.append(na)
        list_nB.append(nb)
        list_bA.append(bga)
        list_bB.append(bgb)
        list_Pm.append(Pm)
        list_nPA.append(na * bga**2 * Pm)
        list_nPB.append(nb * bgb**2 * Pm)

    return list_zbin, list_nA, list_nB, list_bA, list_bB, list_Pm, list_nPA, list_nPB


# ---------------------------------------------------------------------------
# Convenience plotting functions
# ---------------------------------------------------------------------------

def plot_nbP_single(zarray, nz, bz, Area, N_degm2,
                    k_values=(0.01, 0.05, 0.1),
                    Deltaz=0.2, cosmo=None,
                    label='tracer', ax=None):
    """
    Plot nP = n * b^2 * P(k, z) as a function of redshift for several k values.

    Parameters
    ----------
    k_values : tuple/list of k values to plot  [h/Mpc]
    ax       : matplotlib axis (created if None)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))

    for k in k_values:
        zbins, _, _, _, nP = compute_nbP(
            zarray, nz, bz, Area, N_degm2,
            k=k, Deltaz=Deltaz, cosmo=cosmo
        )
        ax.plot(zbins, nP, marker='o', markersize=3,
                label=f'{label},  k={k:.3f} h/Mpc')

    ax.axhline(1.0, color='k', ls='--', lw=0.8, label='nP = 1')
    ax.set_xlabel('redshift z')
    ax.set_ylabel(r'$n\,b^2\,P(k,z)$')
    ax.set_title(r'Signal-to-noise ratio $nP$ per bin')
    ax.legend(fontsize=8)
    ax.set_yscale('log')
    plt.tight_layout()
    return ax


def plot_nbP_two_tracers(zarray, nza, nzb, bza, bzb,
                          Area, Na_degm2, Nb_degm2,
                          k_values=(0.01, 0.05, 0.1),
                          Deltaz=0.2, cosmo=None,
                          label_a='tracer A', label_b='tracer B',
                          ax=None):
    """
    Plot nP = n * b^2 * P(k, z) for two tracers as a function of redshift.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))

    colors_a = plt.cm.Blues(np.linspace(0.4, 0.9, len(k_values)))
    colors_b = plt.cm.Reds(np.linspace(0.4, 0.9, len(k_values)))

    for ki, k in enumerate(k_values):
        zbins, _, _, _, _, _, nPA, nPB = compute_nbP_two_tracers(
            zarray, nza, nzb, bza, bzb,
            Area, Na_degm2, Nb_degm2,
            k=k, Deltaz=Deltaz, cosmo=cosmo
        )
        ax.plot(zbins, nPA, marker='o', markersize=3, color=colors_a[ki],
                label=f'{label_a},  k={k:.3f} h/Mpc')
        ax.plot(zbins, nPB, marker='s', markersize=3, color=colors_b[ki],
                label=f'{label_b},  k={k:.3f} h/Mpc')

    ax.axhline(1.0, color='k', ls='--', lw=0.8, label='nP = 1')
    ax.set_xlabel('redshift z')
    ax.set_ylabel(r'$n\,b^2\,P(k,z)$')
    ax.set_title(r'Signal-to-noise ratio $nP$ per bin')
    ax.legend(fontsize=7, ncol=2)
    ax.set_yscale('log')
    plt.tight_layout()
    return ax