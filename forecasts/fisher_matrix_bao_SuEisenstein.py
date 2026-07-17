import numpy as np
import scipy as sp
#import pyccl as ccl
import math
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from scipy import integrate
from functools import partial
from scipy.integrate import quad, dblquad

import pyccl as ccl
import camb
from camb import model, initialpower
import cosmology

deltac=1.686
c_ls=300*10**3
nlim=10000


#(na,nb,ba,bb,zeff,p,Vsur,kmin,kmax,mod, cosmo=None):
def Mat_Fisher_BAO(n,beff,zeff,Vsur, cosmo=None):
    '''
    return the 2x2pt Fisher matrix for BAO,
    following the Sue&Eisenstein equation.
    0,0 for Da
    1,1 for H
    '''
    kmin=2*math.pi/Vsur**(1/3)
    kmax=0.1*cosmology.D(0,cosmo)/cosmology.D(zeff,cosmo)#/h
    
    z=zeff
    bg=beff


    Omega_m=cosmo['Omega_m']
    Omega_b=cosmo['Omega_b']
    Omega_c=cosmo['Omega_c']
    Omega_l=1-Omega_m
    h = cosmo['h']
    ns = cosmo['n_s']



    
    ksilk=1.6*(Omega_b*h**2)**0.52*(Omega_m*h**2)**0.73*(1+(10.4*Omega_m*h**2)**(-0.95))/h
    A0=1#0.48
    G=cosmology.D(z,cosmo)*0.758
    Sigma_s=1/ksilk
    #print(Sigma_s)
    Sigma_perp=9.4*cosmology.sigma_8(z,cosmo)/0.9
    Sigma_parr=Sigma_perp*(1+cosmology.f(z,cosmo))
    
    def R(k,mu):
        return (bg+cosmology.f(z,cosmo)*mu**2)**2*math.exp(-k**2*mu**2*cosmology.sigma_chi(z,cosmo)**2)
        
    def Fij_BAO(k,mu,int2):
        nk=n
        if int2==11:
            return (mu**2-1)**2*k**2*math.exp(-2*(k*Sigma_s)**1.4)/(cosmology.Pm(k,z,cosmo)/cosmology.Pm(0.2,z,cosmo)+1/(nk*cosmology.Pm(0.2,z,cosmo)*R(k,mu)))**2*math.exp(-k**2*(1-mu**2)*Sigma_perp**2-k**2*mu**2*Sigma_parr**2)
        if int2==12:
            return (mu**2-1)*mu**2*k**2*math.exp(-2*(k*Sigma_s)**1.4)/(cosmology.Pm(k,z,cosmo)/cosmology.Pm(0.2,z,cosmo)+1/(nk*cosmology.Pm(0.2,z,cosmo)*R(k,mu)))**2*math.exp(-k**2*(1-mu**2)*Sigma_perp**2-k**2*mu**2*Sigma_parr**2)
        if int2==22:
            return mu**4*k**2*math.exp(-2*(k*Sigma_s)**1.4)/(cosmology.Pm(k,z,cosmo)/cosmology.Pm(0.2,z,cosmo)+1/(nk*cosmology.Pm(0.2,z,cosmo)*R(k,mu)))**2*math.exp(-k**2*(1-mu**2)*Sigma_perp**2-k**2*mu**2*Sigma_parr**2)
            #print(Fij_fnl_bg(0.1,0.5,11))
    def integrat_Fk_BAO(muint,int1):
        if int1 ==11:
            return integrate.quad(partial(Fij_BAO,mu=muint,int2=int1),kmin,kmax,epsrel=0.00001,epsabs=0.000001,limit=nlim)[0]
        if int1 ==12:
            return integrate.quad(partial(Fij_BAO,mu=muint,int2=int1),kmin,kmax,epsrel=0.00001,epsabs=0.00001,limit=nlim)[0]
        if int1 ==22:
            return integrate.quad(partial(Fij_BAO,mu=muint,int2=int1),kmin,kmax,epsrel=0.00001,epsabs=0.00001,limit=nlim)[0]
    def F_integrat_BAO():
        f11=Vsur*A0**2*quad(partial(integrat_Fk_BAO,int1=11),0, 1,epsrel=0.00001,epsabs=0.00001,limit=nlim)[0]
        f12=Vsur*A0**2*quad(partial(integrat_Fk_BAO,int1=12),0, 1,epsrel=0.00001,epsabs=0.00001,limit=nlim)[0]
        f22=Vsur*A0**2*quad(partial(integrat_Fk_BAO,int1=22),0, 1,epsrel=0.00001,epsabs=0.00001,limit=nlim)[0]
        return np.array([[f11,f12],[f12,f22]])

    #[[G,H],[I,J]]=np.linalg.inv(F_integrat_BAO())
    #mat=np.linalg.inv(F_integrat_BAO())
    #print('constrains on Da and H ', G**0.5,' and ',J**0.5)
    #return([G**0.5,J**0.5])
    return F_integrat_BAO()



def sigma_Da_H_single_tracer(zarray,nz,bz,Area,N_degm2,Deltaz=0.2, cosmo=None):
    '''
    return 3 array zbins, sigma(Da), sigma(H) corresponding to measurements for bins of width Dz for which sigma is finite. 
    also returns 3 values, zeff, sigma(Da)_eff,sigma(H)_eff corresponding to the combined measurements. 
    
    zarray,nz,bz are 3 arrays, nz is the normalised redshift distribution, bz the galaxy bias
    zarray is the center of the bin, so it can't start at 0!!!!!
    Area, the area in deg2
    N_degm2 the number of spec-z per deg2.
    
    '''
    if zarray[0]==0:
        print('error convention of z bin cannot start at 0!')
        return 0,0,0,0
        
    eps=0.0001
    dz_array=(zarray[1]-zarray[0])
    Nbin= int((zarray[-1]+dz_array-zarray[0]+eps)//Deltaz)
    #print(Nbin,' bins')
    
    list_zbin=[]
    list_sigma_Da=[]
    list_sigma_H=[]
    Flist=[]
    
    #normalised the nz
    nz=nz/np.sum(nz)
    size_z=len(zarray)
    
    for i in range(Nbin):
        imin,imax=i*int((size_z+eps)//Nbin),(i+1)*int((eps+size_z)//Nbin)
        nzsum=np.sum(nz[imin:imax])
        if nzsum>0: #if they are galaxies for this bin, let's do a forecast 
            
            zbin=(zarray[imin]-dz_array/2+Deltaz/2)
            #print(round(zbin,2),' zbin')
            Vsur=cosmology.Vsurvey(zbin-Deltaz/2,Deltaz,Area,cosmo)
            
            bg=np.sum(nz[imin:imax]*bz[imin:imax])/np.sum(nz[imin:imax])
            
            n=nzsum*N_degm2*Area/Vsur
            
            list_zbin.append(zbin)
            F=Mat_Fisher_BAO(n,bg,zbin,Vsur, cosmo=cosmo)
            Flist.append(F)
            list_sigma_Da.append((np.linalg.inv(F)[0][0])**0.5)
            list_sigma_H.append((np.linalg.inv(F)[1][1])**0.5)
        
    zeff=np.sum(zarray*nz)
    Flist=np.array(Flist)
    Ftot=np.sum(Flist,axis=0)
    sigma_Da_eff=(np.linalg.inv(Ftot)[0][0])**0.5
    sigma_H_eff=(np.linalg.inv(Ftot)[1][1])**0.5
    return list_zbin, list_sigma_Da,list_sigma_H, zeff, sigma_Da_eff, sigma_H_eff

def sigma_w0_wa_eff(zeff, sigma_Da_eff, sigma_H_eff, cosmo, w0=-1.0, wa=0.0, corr_DaH=0.0, rel_step=1e-2):

    Omega_m_fid = cosmo['Omega_m']
    Omega_b_fid = cosmo['Omega_b']
    h_fid       = cosmo['h']
    n_s_fid     = cosmo['n_s']
    sigma8_fid  = cosmo['sigma8']

    params_fid = {
        'w0'     : w0,
        'wa'     : wa,
        'Omega_m': Omega_m_fid,
    }

    def _make_cosmo(p):
        Om  = p['Omega_m']
        Oc  = Om - Omega_b_fid          # keep Omega_b fixed, vary Omega_cdm
        return ccl.Cosmology(
            Omega_c  = Oc,
            Omega_b  = Omega_b_fid,
            h        = h_fid,
            n_s      = n_s_fid,
            sigma8   = sigma8_fid,
            w0       = p['w0'],
            wa       = p['wa'],
        )

    def _observables(p):
        c = _make_cosmo(p)
        a = 1.0 / (1.0 + zeff)
        # comoving angular diameter distance [Mpc]
        chi  = ccl.comoving_radial_distance(c, a)   # Mpc
        Da   = chi / (1.0 + zeff)                   # Mpc
        # Hubble parameter [km/s/Mpc]
        H    = ccl.h_over_h0(c, a) * 100.0 * h_fid # km/s/Mpc
        return np.array([Da, H])

    param_names = ('w0', 'wa')
    J = np.zeros((2, len(param_names)))

    for j, par in enumerate(param_names):
        step = rel_step

        p_plus  = dict(params_fid); p_plus[par]  += step
        p_minus = dict(params_fid); p_minus[par] -= step

        J[:, j] = (_observables(p_plus) - _observables(p_minus)) / (2.0 * step)

    O_fid = _observables(params_fid)
    Da_fid, H_fid_val = O_fid

    abs_sigma_Da = sigma_Da_eff * Da_fid 
    abs_sigma_H  = sigma_H_eff  * H_fid_val 

    C_obs = np.array([
        [abs_sigma_Da**2, corr_DaH * abs_sigma_Da * abs_sigma_H],
        [corr_DaH * abs_sigma_Da * abs_sigma_H,  abs_sigma_H**2],])
    C_obs_inv = np.linalg.inv(C_obs)

    Fisher_mat = J.T @ C_obs_inv @ J

    return Fisher_mat

