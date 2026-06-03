# import glob
import math
# import multiprocessing as mp
import numpy as np
from numpy.polynomial import Polynomial
from numpy.random import default_rng
from pathlib import Path
import pickle
import matplotlib.pyplot as plt
import scipy.optimize
# import subprocess
# from astropy import constants as const
from astropy.coordinates import SkyCoord
from astropy.table import Table, join
from astropy.wcs import WCS
from astropy.io import fits
# from astropy import units as u
import treecorr
# import pdb
from utils import get_builtin_cosmology
from data_io import *
import pdb
import scipy.integrate
import scipy.special

##########################################################

# Plotting parameters
size_ratio = 1
sub_width, sub_height = size_ratio*10/3, size_ratio*2.8
SMALL_SIZE = 7
MEDIUM_SIZE = 9
BIGGER_SIZE = 10

rc_default = {}
rc_default['font.family'] = 'serif'
rc_default['font.size'] = SMALL_SIZE
rc_default['axes.labelsize'] = MEDIUM_SIZE
rc_default['axes.labelweight'] = 'normal'
rc_default['axes.linewidth'] = 0.8

rc_default['axes.titlesize'] = MEDIUM_SIZE
rc_default['xtick.labelsize'] = MEDIUM_SIZE
rc_default['ytick.labelsize'] = MEDIUM_SIZE
rc_default['legend.fontsize'] = SMALL_SIZE
rc_default['figure.titlesize'] = MEDIUM_SIZE
rc_default['lines.linewidth'] = 0.8
rc_default['lines.markersize'] = 3
rc_default['figure.figsize'] = (sub_width, sub_height)
rc_default['savefig.dpi'] = 450

# Latex related
rc_default['text.usetex'] = True
rc_default['mathtext.fontset'] = 'custom'
rc_default['mathtext.rm'] = 'Bitstream Vera Sans'
rc_default['mathtext.it'] = 'Bitstream Vera Sans:italic'
rc_default['mathtext.bf'] = 'Bitstream Vera Sans:bold'

plt.rcParams.update(rc_default)
plt.style.use('tableau-colorblind10')

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

markers = ['s', 'o', 'D', 'x', '>', '<', 'p', '*', 'h', '+', 'x']



def be_fit(z, zc, alpha, beta, norm):
     """Generalised Baugh & Efstathiou (1993, eqn 7) model for N(z)."""
     return norm * z**alpha * np.exp(-(z/zc)**beta)


def nz_old(infile='WAVES-N_0p2_Z22_GalsAmbig_CompletePhotoZ.fits',
    magbins=np.linspace(15, 22, 8),
    zbins=np.linspace(0.0, 2.0, 41), outfile='NzN.pkl',
    plot=False):
    """Plot observed and predicted N(z) histograms in mag slices."""
    
    t = Table.read(infile)
    sel = (t['mag_r'] < magbins[-1]) & (t['crs_sel_v2'])
    t = t[sel]
    mag, z = t['mag_r'], t['Z']
    zcen = zbins[:-1] + 0.5*np.diff(zbins)
    zmin, zmax = zbins[0], zbins[-1]
    zp = np.linspace(zmin, zmax, 500)
    Nz_dict = {'zbins': zbins, 'zcen': zcen}
    plt.clf()
    # plt.figure(figsize=(15, 10))
    color_cycle = iter(colors)  # Use the `colors` variable defined earlier
    be_pars = np.zeros((len(magbins)-1, 4))

    for imag in range(len(magbins) - 1):
        mlo, mhi = magbins[imag], magbins[imag+1]
        sel = (magbins[imag] <= mag) * (mag < magbins[imag+1])
        color = next(color_cycle)  # Get the next color from the cycle
        counts, edges = np.histogram(z[sel], zbins)
        popt, pcov = scipy.optimize.curve_fit(
            be_fit, zcen, counts, p0=(0.5, 2.0, 1.5, 1e6), ftol=1e-3, xtol=1e-3,maxfev=10000)
        print(popt)
        be_pars[imag, :] = popt
        Nz_dict.update({imag: (mlo, mhi, counts, popt)})
        plt.stairs(counts, edges, color=color, label=f"m = {mlo}, {mhi}]")
        plt.plot(zp, be_fit(zp, *popt), color=color, ls='-')
    Nz_dict.update({'be_pars': be_pars})
    Nz_dict.update({'mbins': magbins})
    pickle.dump(Nz_dict, open(outfile, 'wb'))
    plt.legend()
    plt.xlabel('z')
    plt.ylabel('N(z)')
    plt.title('CRS-BG selection DESI DR1 observed N(z) in mag slices')
    plt.savefig('CRS_selection_nz_bg.png', bbox_inches='tight')
    plt.show()


class Nz:
    def __init__(self, data=None, magbins=None, zbins=None, be_pars=None,
                 infile=None, outfile=None, mag_column='mag_r', z_column='Z',
                 Nz_dict=None, Nz_dict_file=None, zdata=None, magdata=None):
        self.data = data
        self.magbins = magbins
        self.zbins = zbins
        self.be_pars = be_pars
        self.infile = infile
        self.outfile = outfile
        self.mag_column = mag_column
        self.z_column = z_column
        self.Nz_dict = Nz_dict
        self.Nz_dict_file = Nz_dict_file
        self.zdata = zdata
        self.magdata = magdata
        self.counts_mag = None
        self.counts = None
        self.popt_mag = None
        self.pcov_mag = None
        self.zcen = None
        self.Nz_dict_full = {}
        self.be_par_full = None

        if self.Nz_dict is None:
            if self.Nz_dict_file is not None:
                self.Nz_dict = pickle.load(open(self.Nz_dict_file, 'rb'))

        if self.data is not None:
            if zdata is None:
                if z_column in data.colnames:
                    self.zdata = data[z_column]
                elif self.zdata is not None:
                    pass
                else:
                    print(f"Warning: {z_column} not found in data columns.")

            if magdata is None:
                if mag_column in data.colnames:
                    self.magdata = data[mag_column]
                elif self.magdata is not None:
                    pass  # Use provided magdata
                else:
                    print(f"Warning: {mag_column} not found in data columns.")

        if self.Nz_dict is not None:
            if self.magbins is None:
                if 'mbins' in self.Nz_dict:
                    self.magbins = self.Nz_dict['mbins']
                else:
                    print("Warning: magbins not found in Nz_dict.")
            
            if self.zbins is None:
                if 'zbins' in self.Nz_dict:
                    self.zbins = self.Nz_dict['zbins']
                else:
                    print("Warning: zbins not found in Nz_dict.")
            if self.zcen is None:
                if 'zcen' in self.Nz_dict:
                    self.zcen = self.Nz_dict['zcen']
                else:
                    print("Warning: zcen not found in Nz_dict.")

            if self.be_pars is None:
                if 'be_pars' in self.Nz_dict:
                    self.be_pars = self.Nz_dict['be_pars']
                else:
                    print("Warning: be_pars not found in Nz_dict.")

            for imag in range(len(self.magbins) - 1):
                if imag in self.Nz_dict.keys():
                    try:
                        mlo, mhi, counts, popt, pcov = self.Nz_dict[imag]
                    except ValueError:
                        print(f"Warning: Nz_dict entry {imag} does not have the expected format.")
                        pass
                else:
                    print(f"Warning: Nz_dict entry {imag} not found.")
                    continue
            if 'nz_full' in self.Nz_dict.keys():
                self.Nz_dict_full = self.Nz_dict['nz_full']
                if 'be_pars' in self.Nz_dict_full.keys():
                    self.be_par_full = self.Nz_dict_full['be_pars']
                else:
                    print("Warning: be_par_full not found in Nz_dict.")
        else:
            self.Nz_dict = {}
        
        if self.zcen is None:
            if self.zbins is not None:
                self.zcen = self.zbins[:-1] + 0.5 * np.diff(self.zbins)
            else:
                print("Warning: zbins not provided, cannot compute zcen.")
    

    def fit_be_mag(self, magbins=None, zbins=None, plot=False, p0=(0.5, 2.0, 1.5, 1e6), 
                   ax=None, return_be_pars=True, return_dict=False, magdata=None, 
                   zdata=None, show=True, outfile='NzN_mag_slices.plk', save=True, 
                   keep=True, start_bin=0, ftol=1e-3, maxfev=10000, verbose=True,
                   line_style='-', plot_be=True, plot_step=True, colourmap=None):

        if magbins is None:
            try:
                magbins = self.magbins
            except AttributeError:
                raise ValueError("magbins must be provided or set in the Nz object.")
            
        if zbins is None:
            try:
                zbins = self.zbins
                
            except AttributeError:
                raise ValueError("zbins must be provided or set in the Nz object.")
            
        if magdata is None:
            try:
                magdata = self.magdata
            except AttributeError:
                raise ValueError("magdata must be provided or set in the Nz object.")
        if zdata is None:
            try:
                zdata = self.zdata
            except AttributeError:
                raise ValueError("zdata must be provided or set in the Nz object.")
        
        mag, z = magdata, zdata
        zcen = zbins[:-1] + 0.5*np.diff(zbins)
        zmin, zmax = zbins[0], zbins[-1]
        zp = np.linspace(zmin, zmax, 500)
        be_pars = np.zeros((len(magbins)-1, 4))

        Nz_dict = {'zbins': zbins, 'zcen': zcen}
        if plot:
            if ax is None:
                fig, ax = plt.subplots()
        if colourmap is None:
            colourmap = colors

        for imag in range(start_bin, len(magbins) - 1):
            mlo, mhi = magbins[imag], magbins[imag + 1]
            sel = (mag >= mlo) & (mag < mhi)
            counts, edges = np.histogram(z[sel], bins=zbins)
            popt, pcov = scipy.optimize.curve_fit(be_fit, zcen, counts, p0=p0, ftol=ftol, maxfev=maxfev)
            be_pars[imag, :] = popt
            Nz_dict.update({imag: (mlo, mhi, counts, popt, pcov)})
            color = colourmap[imag]  # Use the next color from the cycle
            if plot:
                if plot_step:
                    ax.stairs(counts, edges, label=fr'${mlo:.2f} < mag < {mhi:.2f}$', alpha=0.5, color=color)
                if plot_be:
                    ax.plot(zp, be_fit(zp, *popt), lw=2, color=color, ls=line_style)
        Nz_dict.update({'be_pars': be_pars,
                        'magbins': magbins})
        if verbose:
            print("Baugh & Efstathiou parameters for magnitude slices:")
            for i, pars in enumerate(be_pars):
                print(f"Slice {i}: zc={pars[0]:.3f}, alpha={pars[1]:.3f}, beta={pars[2]:.3f}, norm={pars[3]:.3e}")
        if save:
            pickle.dump(Nz_dict, open(outfile, 'wb'))
            if verbose:
                print(f"Saved N(z) data to {outfile}")
            
        if plot:
            if show:
                ax.legend()
                ax.set_xlabel('Redshift (z)')
                ax.set_ylabel('N(z)')
                ax.set_title('Observed N(z) in Magnitude Slices')
                plt.tight_layout()
                plt.show()
        
        if keep:
            self.Nz_dict = Nz_dict
            self.magbins = magbins
            self.zbins = zbins
            self.be_pars = be_pars
            self.zcen = zcen
            
        if return_dict:
            if return_be_pars:
                return Nz_dict, be_pars
            else:
                return Nz_dict
        elif return_be_pars:
            return be_pars
    
    
    def _be_fit(z, zc, alpha, beta, norm):
        """
        Generalised Baugh & Efstathiou (1993, eqn 7) model for N(z).
        Args:
            z: Redshift.
            zc: Characteristic redshift.
            alpha: Power-law index.
            beta: Exponential decay parameter.
            norm: Normalization factor.
        Returns:
            N(z): Number counts at redshift z.
        """
        return norm * z**alpha * np.exp(-(z/zc)**beta)


    def get_be_pars(self, imag=None):
        """
        Get Baugh & Efstathiou parameters for a specific magnitude slice.

        Args:
            imag: Index of the magnitude slice (default: None, returns all).

        Returns:
            be_pars: Parameters for the specified magnitude slice or all if imag is None.
        """
        if self.be_pars is not None:
            if imag is None:
                return self.be_pars
            else:
                return self.be_pars[imag]
        else:
            raise ValueError("Baugh & Efstathiou parameters not set. Run fit_be_mag first.")
    

    def get_be_fit(self, z, imag=None):
        """
        Get Baugh & Efstathiou fit for a specific redshift and magnitude slice.

        Args:
            z: Redshift value or array.
            imag: Index of the magnitude slice (default: None, uses all).

        Returns:
            fit: Baugh & Efstathiou fit for the specified redshift and magnitude slice.
        """
        if self.be_pars is not None:
            if imag is None:
                return self._be_fit(z, *self.be_par_full)
            else:
                return self._be_fit(z, *self.be_pars[imag])
        else:
            raise ValueError("Baugh & Efstathiou parameters not set. Run fit_be_mag first.") 


    def save_be_pars(self, outfile='be_pars.pkl'):
        """
        Save Baugh & Efstathiou parameters to a file.

        Args:
            outfile: Output file name (default: 'be_pars.pkl').
        """
        if self.be_pars is not None:
            with open(outfile, 'wb') as f:
                pickle.dump(self.be_pars, f)
            print(f"Baugh & Efstathiou parameters saved to {outfile}")
        else:
            raise ValueError("Baugh & Efstathiou parameters not set. Run fit_be_mag first.")
    

    def save_Nz_dict(self, outfile='NzN_mag_slices.plk'):
        """
        Save N(z) dictionary to a file.

        Args:
            outfile: Output file name (default: 'NzN_mag_slices.pkl').
        """
        if self.Nz_dict is not None:
            with open(outfile, 'wb') as f:
                pickle.dump(self.Nz_dict, f)
            print(f"N(z) dictionary saved to {outfile}")
        else:
            raise ValueError("N(z) dictionary not set. Run fit_be_mag first.")
    

    def save(self, outfile='Nz_object.pkl'):
        """
        Save the Nz object to a file.

        Args:
            outfile: Output file name (default: 'Nz_object.pkl').
        """
        with open(outfile, 'wb') as f:
            pickle.dump(self, f)
        print(f"Nz object saved to {outfile}")


    def fit_be(self, zbins=None, plot=False, p0=(0.5, 2.0, 1.5, 1e6), 
               ax=None, return_be_pars=True, return_dict=False, 
               zdata=None, show=True, outfile='NzN.plk', save=True, 
               keep=True, start_bin=0, ftol=1e-3, maxfev=10000, verbose=True,
               line_style='-', plot_be=True, plot_step=True, colourmap=None):

        if zbins is None:
            try:
                zbins = self.zbins
            except AttributeError:
                raise ValueError("zbins must be provided or set in the Nz object.")
        if zdata is None:
            try:
                zdata = self.zdata
            except AttributeError:
                raise ValueError("zdata must be provided or set in the Nz object.")
        z = zdata
        zcen = zbins[:-1] + 0.5 * np.diff(zbins)
        zmin, zmax = zbins[0], zbins[-1]
        zp = np.linspace(zmin, zmax, 500)

        Nz_dict_full = {'zbins': zbins, 'zcen': zcen}
        be_pars = np.zeros((1, 4))
        if plot:
            if ax is None:
                fig, ax = plt.subplots()
        counts, edges = np.histogram(z, bins=zbins)
        popt, pcov = scipy.optimize.curve_fit(be_fit, zcen, counts, p0=p0, ftol=ftol, maxfev=maxfev)
        be_pars[0, :] = popt
        Nz_dict_full.update({'be_pars': be_pars, 'zbins': zbins, 'counts': counts, 'edges': edges})
        if save:
            pickle.dump(Nz_dict_full, open(outfile, 'wb'))
            if verbose:
                print(f"Saved N(z) data to {outfile}")
        if colourmap is not None:
            color = colourmap[0] 
        else:
            color = colors[0] # Use the next color from the cycle
        if plot:
            if plot_step:
                ax.stairs(counts, edges, alpha=0.5, color=color)
            if plot_be:
                ax.plot(zp, be_fit(zp, *popt), lw=2, color=color, ls=line_style)
        if keep:
            self.Nz_dict_full = Nz_dict_full
            self.be_par_full = be_pars
            self.Nz_dict.update({'nz_full': Nz_dict_full})

        if return_dict:
            if return_be_pars:
                return Nz_dict_full, be_pars
            else:
                return Nz_dict_full
            


def w_lum_Nz_fit(cosmo, thetabins, m, r0, gamma, be_pars, plotint=0, pdf=None, plot_den=0, verbose=0):
    """
    Compute w(theta) for luminosity-dependent xi(r) correlation length and index.

    This function evaluates Maddox+1996 eqn 31 using a Baugh & Efstathiou (BE) fit
    to number counts N(z) and supplied xi(r) parameters for each magnitude slice.

    Args:
        cosmo: Cosmology object with methods for distance calculations.
        thetabins: Array of angular separation bins (in degrees).
        m: Magnitude slice.
        r0: Correlation length.
        gamma: Power-law index for xi(r).
        be_pars: Parameters for the BE fit.
        plotint: If > 0, enables plotting of intermediate results.
        pdf: PDF object for saving plots (optional).
        plot_den: If > 0, plots the denominator function.
        verbose: Verbosity level for debugging.

    Returns:
        List of w(theta) values for each angular separation bin.
    """
    def denfun(z):
        """Denominator of Maddox+ eqn 31."""
        return be_fit(z, *be_pars)

    def xifun(z1, z2):
        """Numerator of Maddox+ eqn 31."""
        z = 0.5 * (z1 + z2)
        x1 = cosmo.dc(z1)
        x2 = cosmo.dc(z2)
        x12 = np.sqrt(x1**2 + x2**2 - 2 * x1 * x2 * np.cos(np.deg2rad(theta)))
        return denfun(z1) * denfun(z2) * (r0 / x12)**gamma

    zmin, zmax = cosmo._z[0], cosmo._z[-1]
    z = np.linspace(zmin, zmax, 100)

    if pdf:
        if plot_den:
            plt.figure()
            plt.plot(z, denfun(z))
            plt.xlabel('Redshift')
            plt.ylabel('Denominator Function')
            plt.title(f'Magnitude = {m:.2f}')
            pdf.savefig()
            plt.close()

        plt.figure()
        z1, z2 = np.meshgrid(z, z)
        plt.imshow(xifun(z1, z2), extent=(zmin, zmax, zmax, zmin), norm='log')
        plt.xlabel('z2')
        plt.ylabel('z1')
        plt.title(f'Magnitude = {m:.2f}, Theta = {theta:.2f}')
        plt.colorbar()
        pdf.savefig()
        plt.close()

    denom = scipy.integrate.quad(denfun, zmin, zmax, epsrel=0.01)[0]**2
    if verbose > 0:
        print('Denominator:', m, denom)

    wbins = []
    for theta in thetabins:
        num = scipy.integrate.dblquad(xifun, zmin, zmax, lambda _: zmin, lambda _: zmax, epsrel=0.01)[0]
        if np.isnan(num):
            pdb.set_trace()
        if verbose > 0:
            print('Numerator:', m, theta, num)

        wbins.append(num / denom)
    return wbins


def w_a(cosmo, be_pars, gamma=1.7, r0=5.0, eps=0, plotint=0):
    """
    Compute w(theta) amplitude for a power-law xi(r) and given selection functions.

    Args:
        cosmo: Cosmology object with methods for distance calculations.
        be_pars: Parameters for the BE fit.
        gamma: Power-law index for xi(r).
        r0: Correlation length.
        eps: Evolution parameter.
        plotint: If > 0, enables plotting of intermediate results.

    Returns:
        Amplitude of w(theta).
    """
    def denfun(z):
        """Denominator of Maddox+ eqn 35a."""
        return be_fit(z, *be_pars)

    def xifun(z):
        """Numerator of Maddox+ eqn 35a."""
        x = cosmo.dc(z)
        return x**(1 - gamma) * denfun(z)**2 * (1 + z)**(gamma - 3 - eps) / cosmo.dxdz(z)

    zmin, zmax = cosmo._z[0], cosmo._z[-1]
    gfac = (math.pi**0.5 * scipy.special.gamma((gamma - 1) / 2) * r0**gamma /
            scipy.special.gamma(gamma / 2))

    if plotint:
        zp = np.linspace(zmin, zmax, 100)
        plt.figure()
        plt.plot(zp, denfun(zp))
        plt.xlabel('z')
        plt.ylabel('Denominator Function')
        plt.show()

        plt.figure()
        plt.plot(zp, xifun(zp))
        plt.xlabel('z')
        plt.ylabel('Numerator Function')
        plt.show()

    num = scipy.integrate.quad(xifun, zmin, zmax, epsabs=1e3, epsrel=1e-3)[0]
    den = scipy.integrate.quad(denfun, zmin, zmax, epsabs=1e3, epsrel=1e-3)[0]**2
    B = num / den
    A = gfac * B
    return A


def limber_scale(cosmo, be_pars, be_pars_ref, gamma1=1.7, gamma2=4, r0=5.0, eps=0,
                 gamma_1_ref=1.7, gamma_2_ref=4):
    """
    Scale w(theta) to a reference depth for a two power-law model.

    Args:
        cosmo: Cosmology object with methods for distance calculations.
        be_pars: Parameters for the BE fit for the current depth.
        be_pars_ref: Parameters for the BE fit for the reference depth.
        gamma1: Power-law index for the first model.
        gamma2: Power-law index for the second model.
        r0: Correlation length.
        eps: Evolution parameter.
        gamma_1_ref: Reference power-law index for the first model.
        gamma_2_ref: Reference power-law index for the second model.

    Returns:
        Tuple of (log10(theta scale), log10(w(theta) scale)).
    """
    A = w_a(cosmo, be_pars, gamma=gamma_1_ref, r0=r0, eps=eps)
    B = w_a(cosmo, be_pars, gamma=gamma_2_ref, r0=r0, eps=eps)
    Aref = w_a(cosmo, be_pars_ref, gamma=gamma_1_ref, r0=r0, eps=eps)
    Bref = w_a(cosmo, be_pars_ref, gamma=gamma_2_ref, r0=r0, eps=eps)
    dlgA, dlgB = np.log10(A/Aref), np.log10(B/Bref)
    dlgt = (dlgA - dlgB) / (gamma_2_ref - gamma_1_ref)
    dlgw = (gamma_1_ref-1)*dlgt - dlgA
    return dlgt, dlgw


def apply_scaling(dlgt, dlgw, theta, wtheta, error=None):
    """
    Apply scaling to theta and w(theta) based on computed scaling factors.

    Args:
        dlgt: Log10(theta scale).
        dlgw: Log10(w(theta) scale).
        theta: Array of angular separations.
        wtheta: Array of w(theta) values.
        error: Array of errors for w(theta) (optional).

    Returns:
        Tuple of (scaled theta, scaled w(theta), scaled error).
    """
    tscale = 10**dlgt
    wscale = 10**dlgw

    theta_scaled = theta * tscale
    wtheta_scaled = wscale * wtheta
    error_scaled = error * wscale if error is not None else None

    return theta_scaled, wtheta_scaled, error_scaled

def compute_absolute_mag_single(mag, z, coeffs, band_name, kc=None, cosmo=None, responses=None):
    """
    Compute absolute magnitude given apparent magnitude, redshift, and K-correction coefficients.

    Args:
        mag (float): Apparent magnitude.
        z (float): Redshift.
        coeffs (array-like): K-correction coefficients (shape must be (1, 5)).
        band_name (str): Name of the photometric band. Should match one of the responses in the Kcorrect object.
        kc (Kcorrect, optional): Kcorrect object for computing K-corrections. If None, it will be initialized.
        cosmo (Cosmology, optional): Cosmology object for distance modulus calculation. Defaults to Planck18.
        responses (list, optional): List of response functions for initializing Kcorrect.

    Returns:
        abs_mag: float: Absolute magnitude(s).
    """

    if kc is None:
        from kcorrect import Kcorrect
        try:
            kc = Kcorrect(responses=responses)
        except Exception as e:
            raise ValueError(
                "Kcorrect object could not be initialized. "
                "Ensure responses are correctly defined."
            ) from e
 
    if cosmo is None:
        cosmo = get_builtin_cosmology('Planck18')   

    coeffs = np.asarray(coeffs).reshape(1, 5)  # Ensure coefficients are in the correct shape
    band_idx = kc.responses.index(band_name)  # Get the index of the specified band



    # Compute distance modulus
    distmod = cosmo.distmod(z).value

    # Compute K-correction for the specified band
    if type(z) is not np.ndarray:
        z = np.array([z])
    k_correction = kc.kcorrect(redshift=z, coeffs=coeffs)[0, band_idx]

    # Calculate absolute magnitude
    abs_mag = mag - distmod - k_correction
    return abs_mag

def compute_absolute_mag(mag, zgrid, coeffs, band_name, kc=None, cosmo=None, responses=None,
                         redshift_ranges=None, verbose=0):
    """
    Compute absolute magnitudes given apparent magnitudes, redshift grid, 
    K-correction coefficients, and cosmology.

    Args:
        mag (float): Apparent magnitude.
        zgrid (array-like): Redshift grid.
        coeffs (array-like): K-correction coefficients (shape must be (5,) or (1, 5)).
        band_name (str): Name of the photometric band. Should match one of the responses in the Kcorrect object.
        kc (Kcorrect, optional): Kcorrect object for computing K-corrections. If None, it will be initialized.
        cosmo (Cosmology, optional): Cosmology object for distance modulus calculation. Defaults to Planck18.
        responses (list, optional): List of response functions for initializing Kcorrect.
        redshift_ranges (list of tuples, optional): List of (zmin, zmax) tuples for each redshift range.
            If provided, the function will compute absolute magnitudes only for the specified ranges.

    Returns:
        abs_mag (float or array-like): Absolute magnitude(s).
    """
    if kc is None:
        from kcorrect import Kcorrect
        try:
            kc = Kcorrect(responses=responses)
        except Exception as e:
            raise ValueError(
                "Kcorrect object could not be initialized. "
                "Ensure responses are correctly defined."
            ) from e

    if cosmo is None:
        cosmo = get_builtin_cosmology('Planck18')

    zgrid = np.atleast_1d(zgrid)

    coeffs = np.asarray(coeffs)
    if coeffs.shape == (5,):
        coeffs = coeffs.reshape(1, 5)
    elif coeffs.shape != (1, 5):
        raise ValueError(f"coeffs must be shape (5,) or (1, 5), got {coeffs.shape}")

    # Repeat coefficients for each redshift
    coeffs_repeated = np.tile(coeffs, (len(zgrid), 1))  # shape (N, 5)

    try:
        band_idx = kc.responses.index(band_name)
    except ValueError:
        raise ValueError(f"Band '{band_name}' not found in kc.responses: {kc.responses}")

    # K-corrections
    k_correction = kc.kcorrect(redshift=zgrid, coeffs=coeffs_repeated)[:, band_idx]

    # Distance modulus
    distmod = cosmo.distmod(zgrid).value

    # Absolute magnitude
    abs_mag = mag - distmod - k_correction
    
    if len(abs_mag) == 1:
        abs_mag = abs_mag[0]

    return abs_mag


def wtheta_estimate_nz(
    cosmo, be_pars, thetabins, mag, r0, gamma, Nz=None,
    responses=None, redshift_ranges=None, verbose=0):
    """
    Estimate w(theta) using N(z) and cosmological parameters.

    Args:
        cosmo (Cosmology): Cosmology object with methods for distance calculations.
        be_pars (array-like): Parameters for the Baugh & Efstathiou fit.
        thetabins (array-like): Array of angular separation bins (in degrees).
        mag (float): Magnitude slice.
        r0 (float): Correlation length.
        gamma (float): Power-law index for xi(r).
        Nz (Nz, optional): Nz object containing N(z) data.
        responses (list, optional): List of response functions for K-correction.
        redshift_ranges (list of tuples, optional): List of (zmin, zmax) tuples for redshift ranges.
        verbose (int, optional): Verbosity level for debugging. Default is 0.

    Returns:
        list: w(theta) values for each angular separation bin.
    """
    def denfun(z):
        """Denominator of Maddox+ eqn 31."""
        return be_fit(z, *be_pars)

    def xifun(theta, z1, z2):
        """Numerator of Maddox+ eqn 31."""
        z = 0.5 * (z1 + z2)
        x1 = cosmo.dc(z1)
        x2 = cosmo.dc(z2)
        cos_theta = np.cos(np.deg2rad(theta))
        x12_squared = x1**2 + x2**2 - 2 * x1 * x2 * cos_theta
        x12 = np.sqrt(x12_squared)
        return denfun(z1) * denfun(z2) * (r0 / x12)**gamma

    # Determine redshift range
    if redshift_ranges is not None:
        zmin, zmax = redshift_ranges[0], redshift_ranges[-1]
    else:
        try:
            zmin, zmax = cosmo._z[0], cosmo._z[-1]
        except AttributeError:
            raise ValueError(
                "Cosmology object does not have _z attribute. "
                "Ensure it is properly initialized with redshift data."
            )

    # Precompute denominator
    denom = scipy.integrate.quad(denfun, zmin, zmax, epsrel=0.01)[0]**2
    if verbose > 0:
        print(f"Denominator: {mag}, {denom}")

    # Precompute cos(theta) values for efficiency
    cos_theta_values = np.cos(np.deg2rad(thetabins))

    # Compute w(theta) for each angular separation bin
    wbins = []
    for cos_theta in cos_theta_values:
        def xifun_integrated(z1, z2):
            """Numerator of Maddox+ eqn 31 integrated over z1 and z2."""
            x1 = cosmo.dc(z1)
            x2 = cosmo.dc(z2)
            x12_squared = x1**2 + x2**2 - 2 * x1 * x2 * cos_theta
            x12 = np.sqrt(x12_squared)
            return denfun(z1) * denfun(z2) * (r0 / x12)**gamma

        num = scipy.integrate.dblquad(
            xifun_integrated, zmin, zmax, lambda _: zmin, lambda _: zmax, epsrel=0.01
        )[0]
        if np.isnan(num):
            pdb.set_trace()
        if verbose > 0:
            print(f"Numerator: {mag}, {cos_theta}, {num}")

        wbins.append(num / denom)

    return wbins
