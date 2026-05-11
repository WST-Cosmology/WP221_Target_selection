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

def Mat_Fisher_1tracer_fnl_bias(n,bg,z,p,Vsur,kmin,kmax,mod, cosmo=None):
    '''
    return the 2x2 fisher matrix for a single tracer. 0,0 for fnl, 1,1 for bias. 
    all input are scalars
    n density deg-2
    bg the bias
    z the mean redshift of the bin
    p the p scalars in the NG model
    Vsur the volume of the survey in Mpc^3
    kmin kmax ranges of the integral
    '''
    #fiducial
    fNL=0
    prefactor=Vsur/(4*math.pi**2)

    Omega_m=cosmo['Omega_m']
    Omega_b=cosmo['Omega_b']
    Omega_c=cosmo['Omega_c']
    Omega_l=1-Omega_m
    h = cosmo['h']
    ns = cosmo['n_s']
    
    def deltab(k):
        return 3*fNL*(bg-p)*deltac*Omega_m/(k**2*cosmology.T(k,mod,cosmo)*cosmology.D(z,cosmo))*(100*h/c_ls)**2

    def P(k,z1,mu):
        return (bg+deltab(k)+cosmology.f(z1,cosmo)*mu*mu)**2*cosmology.Pm(k,z1,cosmo)

    def ddb_df(k):
        return 3*(bg-p)*deltac*Omega_m/(k**2*cosmology.T(k,mod,cosmo)*cosmology.D(z,cosmo))*(100*h/c_ls)**2
    
    def ddb_db(k):#derivative of deltab by bg
        return 3*fNL*deltac*Omega_m/(k**2*cosmology.T(k,mod,cosmo)*cosmology.D(z,cosmo))*(100*h/c_ls)**2
        
    def Fij_fnl_bg(k,mu,int2):
        pk=P(k,z,mu)
        ddbdf=ddb_df(k)
        ddbdb=ddb_db(k)
        if int2==11:
            return prefactor*2*(n*pk/(n*pk+1))**2*ddbdf**2*k**2/(bg+deltab(k)+cosmology.f(z,cosmo)*mu*mu)**2
        if int2==12:
            return prefactor*2*(n*pk/(n*pk+1))**2*ddbdf*(1+ddbdb)*k**2/(bg+deltab(k)+cosmology.f(z,cosmo)*mu*mu)**2
        if int2==22:
            return prefactor*2*(n*pk/(n*pk+1))**2*(1+ddbdb)**2*k**2/(bg+deltab(k)+cosmology.f(z,cosmo)*mu*mu)**2
        
    def integrat_Fk_fnl_bg(muint,int1):
        return integrate.quad(partial(Fij_fnl_bg,mu=muint,int2=int1),kmin,kmax,epsrel=0.0001,epsabs=0.0001,limit=nlim)[0]
            
    def F_integrat_fnl_bg():
        f11=quad(partial(integrat_Fk_fnl_bg,int1=11),-1, 1,epsrel=0.0001,epsabs=0.0001,limit=nlim)[0]
        f12=quad(partial(integrat_Fk_fnl_bg,int1=12),-1, 1,epsrel=0.0001,epsabs=0.0001,limit=nlim)[0]
        f22=quad(partial(integrat_Fk_fnl_bg,int1=22),-1, 1,epsrel=0.0001,epsabs=0.0001,limit=nlim)[0]
        return np.array([[f11,f12],[f12,f22]])
            
    return F_integrat_fnl_bg()

def Mat_Fisher_2tracer_fnl_bias(na,nb,ba,bb,zeff,p,Vsur,kmin,kmax,mod, cosmo=None):
    '''na,nb,ba,bb,zeff,p,Vsur,kmin,kmax,mod
    return the Fisher matrix 3x3, i=1 and 2 being for biases, and i=0 for fnl.
    zmin<z<zmin+dz is the zrange, with effective z zeff
    na and nb are the effective number densities in (Mpc/h)**-3
    ba and bb are the effective galaxy biases
    sky, the area in deg2
    mod=='camb' or 'bbks' for evaluating transfer function
    p is the usual annoying parameter. std conventions: p=1.6 or p=1 the latter giving better sigma(fnl). 
    '''

    Omega_m=cosmo['Omega_m']
    Omega_b=cosmo['Omega_b']
    Omega_c=cosmo['Omega_c']
    Omega_l=1-Omega_m
    h = cosmo['h']
    ns = cosmo['n_s']
    z=zeff
    f_z=cosmology.f(z)
    prefactor=Vsur/(4*math.pi**2)
    
    def PA(k,mu):
        return (ba+f_z*mu**2)**2*cosmology.Pm(k,z,cosmo)*math.exp(-k**2*mu**2*cosmology.sigma_chi(z,cosmo)**2)
    def PB(k,mu):
        return (bb+f_z*mu**2)**2*cosmology.Pm(k,z,cosmo)*math.exp(-k**2*mu**2*cosmology.sigma_chi(z,cosmo)**2)
    def PAB(k,mu):
        return (ba+f_z*mu**2)*(bb+f_z*mu**2)*cosmology.Pm(k,z,cosmo)*math.exp(-k**2*mu**2*cosmology.sigma_chi(z,cosmo)**2)
    def ddb_df(k,bg):
        return 3*(bg-p)*deltac*Omega_m/(k**2*cosmology.T(k,mod,cosmo)*cosmology.D(z,cosmo))*(100*h/c_ls)**2
    
    def Dt(X,i,mu,k):
        if X=='A':
            return [2/(ba+f_z*mu**2),0,2*ddb_df(k,ba)/(ba+f_z*mu**2)][i]
        if X=='B':
            return [0,2/(bb+f_z*mu**2),2*ddb_df(k,bb)/(bb+f_z*mu**2)][i]
        if X=='AB':
            return [1/(ba+f_z*mu**2),1/(bb+f_z*mu**2),ddb_df(k,bb)/(bb+f_z*mu**2)+ddb_df(k,ba)/(ba+f_z*mu**2)][i]
        print('error DX')
    
    def Raa(k,mu):
        return (na*PA(k,mu)*(1+nb*PB(k,mu))/((1+na*PA(k,mu))*(1+nb*PB(k,mu))-na*nb*PAB(k,mu)**2))**2
    def Rbb(k,mu):
        return (nb*PB(k,mu)*(1+na*PA(k,mu))/((1+nb*PB(k,mu))*(1+na*PA(k,mu))-na*nb*PAB(k,mu)**2))**2
    def Rxx(k,mu):
        return na*nb*((1+na*PA(k,mu))*(1+nb*PB(k,mu))+na*nb*PAB(k,mu)**2)/((1+na*PA(k,mu))*(1+nb*PB(k,mu))-na*nb*PAB(k,mu)**2)**2*PAB(k,mu)**2
    def Rxa(k,mu):
        return na**2*nb*(1+nb*PB(k,mu))/((1+na*PA(k,mu))*(1+nb*PB(k,mu))-na*nb*PAB(k,mu)**2)**2*PAB(k,mu)**2*PA(k,mu)
    def Rxb(k,mu):
        return nb**2*na*(1+na*PA(k,mu))/((1+na*PA(k,mu))*(1+nb*PB(k,mu))-na*nb*PAB(k,mu)**2)**2*PAB(k,mu)**2*PB(k,mu)
    def Rab(k,mu):
        return na**2*nb**2*PA(k,mu)*PB(k,mu)*PAB(k,mu)**2/((1+na*PA(k,mu))*(1+nb*PB(k,mu))-na*nb*PAB(k,mu)**2)**2
    
    def FX(X,k,mu,i,j):
        if X=='A':
            return 1/2*Dt(X,i,mu,k)*Dt(X,j,mu,k)*Raa(k,mu)
        if X=='B':
            return 1/2*Dt(X,i,mu,k)*Dt(X,j,mu,k)*Rbb(k,mu)
        if X=='AB':
            return Dt(X,i,mu,k)*Dt(X,j,mu,k)*Rxx(k,mu)-(Dt(X,i,mu,k)*Dt('A',j,mu,k)+Dt('A',i,mu,k)*Dt(X,j,mu,k))*Rxa(k,mu)-(Dt(X,i,mu,k)*Dt('B',j,mu,k)+Dt('B',i,mu,k)*Dt(X,j,mu,k))*Rxb(k,mu)+1/2*(Dt('A',i,mu,k)*Dt('B',j,mu,k)+Dt('B',i,mu,k)*Dt('A',j,mu,k))*Rab(k,mu)
        print('error Fx')
        
    def integrand(k,mu,i,j):
        return prefactor *(FX('A',k,mu,i,j)+FX('B',k,mu,i,j)+FX('AB',k,mu,i,j))*k**2
    
    def integrat_rsd(muint,int_i,int_j):
        return integrate.quad(partial(integrand,mu=muint,i=int_i,j=int_j),kmin,kmax,epsrel=10**(-3),epsabs=10**(-4),limit=nlim)[0]
    
    def F_rsd():
        
        f = np.zeros((3, 3))
        for i in range(3):
            for j in range(i, 3):  # only compute upper triangle
                res= quad(partial(integrat_rsd, int_i=i, int_j=j),-1, 1, epsrel=10**(-3),epsabs=10**(-4), limit=nlim)[0]
                f[i, j] = res
                f[j, i] = f[i, j]
        return np.rot90(f, 2).T
        
    return F_rsd()

def sigma_fnl_single_tracer(zarray,nz,bz,Area,N_degm2,Deltaz=0.2,p=1,mod='camb',kmax=0.1, cosmo=None):
    '''
    return 2 array zbins, sigma(fnl) corresponding to measurements for bins of width Dz for which sigma is finite. 
    also returns 2 values, zeff, sigma(fnl)_eff corresponding to the combined measurements. 
    
    zarray,nz,bz are 3 arrays, nz is the normalised redshift distribution, bz the galaxy bias
    zarray is the center of the bin, so it can't start at 0!!!!!
    Area, the area in deg2
    N_degm2 the number of spec-z per deg2.
    
    mod=='camb' or 'bbks' for evaluating transfer function
    p is the usual annoying parameter. std conventions: p=1.6 or p=1 the latter giving better sigma(fnl). 
    kmax is the max mode.
    '''
    if zarray[0]==0:
        print('error convention of z bin cannot start at 0!')
        return 0,0,0,0
        
    eps=0.0001
    dz_array=(zarray[1]-zarray[0])
    Nbin= int((zarray[-1]+dz_array-zarray[0]+eps)//Deltaz)
    print(Nbin,' bins')
    
    list_zbin=[]
    list_sigma_fnl=[]
    Flist=[]
    
    #normalised the nz
    nz=nz/np.sum(nz)
    size_z=len(zarray)
    
    for i in range(Nbin):
        imin,imax=i*int((size_z+eps)//Nbin),(i+1)*int((eps+size_z)//Nbin)
        nzsum=np.sum(nz[imin:imax])
        if nzsum>0: #if they are galaxies for this bin, let's do a forecast 
            #print(imin,imax)
            zbin=(zarray[imin]-dz_array/2+Deltaz/2)
            print(round(zbin,2),' zbin')
            Vsur=cosmology.Vsurvey(zbin-Deltaz/2,Deltaz,Area,cosmo)
            #print(Vsur, ' Vsur')
            kmin=2*math.pi/Vsur**(1/3)
            #print(kmin,'kmin')
            bg=np.sum(nz[imin:imax]*bz[imin:imax])/np.sum(nz[imin:imax])
            #print(bg,' bias')
            n=nzsum*N_degm2*Area/Vsur
            #print(n,' density deg-2')

            list_zbin.append(zbin)
            F=Mat_Fisher_1tracer_fnl_bias(n,bg,zbin,p,Vsur,kmin,kmax,mod,cosmo=cosmo)
            Flist.append(F)
            list_sigma_fnl.append((np.linalg.inv(F)[0][0])**0.5)
        
    zeff=np.sum(zarray*nz)
    Flist=np.array(Flist)
    Ftot=np.sum(Flist,axis=0)
    sigma_fnl_eff=(np.linalg.inv(Ftot)[0][0])**0.5
    return list_zbin, list_sigma_fnl, zeff, sigma_fnl_eff

def sigma_fnl_two_tracers(zarray,nza,nzb,bza,bzb,Area,Na_degm2,Nb_degm2,Deltaz=0.2,p=1,mod='camb',kmax=0.1, cosmo=None):
    '''
    return 2 array zbins, sigma(fnl) corresponding to measurements for bins of width Dz for which sigma is finite. 
    also returns 2 values, zeff, sigma(fnl)_eff corresponding to the combined measurements. 
    
    zarray,nza,nzb,bza,bzb are 5 arrays, nz is the normalised redshift distribution, bz the galaxy bias
    zarray is the center of the bin, so it can't start at 0!!!!!
    Area, the area in deg2
    N_degm2 the number of spec-z per deg2.
    
    mod=='camb' or 'bbks' for evaluating transfer function
    p is the usual annoying parameter. std conventions: p=1.6 or p=1 the latter giving better sigma(fnl). 
    kmax is the max mode.
    '''
    if zarray[0]==0:
        print('error convention of z bin cannot start at 0!')
        return 0,0,0,0
        
    eps=0.0001
    dz_array=(zarray[1]-zarray[0])
    Nbin= int((zarray[-1]+dz_array-zarray[0]+eps)//Deltaz)
    print(Nbin,' bins')
    
    list_zbin=[]
    list_sigma_fnl=[]
    Flist=[]
    
    #normalised the nz
    nza=nza/np.sum(nza)
    nzb=nzb/np.sum(nzb)
    size_z=len(zarray)
    
    for i in range(Nbin):
        imin,imax=i*int((size_z+eps)//Nbin),(i+1)*int((eps*size_z)//Nbin)
        nzasum=np.sum(nza[imin:imax])
        nzbsum=np.sum(nzb[imin:imax])
        zbin=(zarray[imin]-dz_array/2+Deltaz/2)
        print(round(zbin,2),' zbin')
        Vsur=cosmology.Vsurvey(zbin-Deltaz/2,Deltaz,Area,cosmo)
        kmin=2*math.pi/Vsur**(1/3)
        if (nzasum>0 and nzbsum>0): #if they are X galaxies for this bin, let's do a double tracer forecas
            print(nzasum,nzbsum)
            bga=np.sum(nza[imin:imax]*bza[imin:imax])/np.sum(nza[imin:imax])
            bgb=np.sum(nzb[imin:imax]*bzb[imin:imax])/np.sum(nzb[imin:imax])
                
            #print(bga,bgb,' biases')
            na=nzasum*Na_degm2*Area/Vsur
            nb=nzbsum*Nb_degm2*Area/Vsur
            #print(na,nb,' density deg-2')

            list_zbin.append(zbin)
            F=Mat_Fisher_2tracer_fnl_bias(na,nb,bga,bgb,zbin,p,Vsur,kmin,kmax,mod,cosmo=cosmo)
            Flist.append(F)
            list_sigma_fnl.append((np.linalg.inv(F)[0,0])**0.5)
        elif nzasum>0:#nzbsum=0
            bg=np.sum(nza[imin:imax]*bza[imin:imax])/np.sum(nza[imin:imax])
            n=nzasum*Na_degm2*Area/Vsur
            
            list_zbin.append(zbin)
            F=Mat_Fisher_1tracer_fnl_bias(n,bg,zbin,p,Vsur,kmin,kmax,mod,cosmo=cosmo)
            #print(F)
            list_sigma_fnl.append((np.linalg.inv(F)[0][0])**0.5)
            F1=np.zeros((3,3))
            F1[0,0],F1[2,0],F1[0,2],F1[2,2]=F[0,0],F[1,0],F[0,1],F[1,1]
            Flist.append(F1)
        elif nzbsum>0:#nzasum=0
            bg=np.sum(nzb[imin:imax]*bzb[imin:imax])/np.sum(nzb[imin:imax])
            n=nzbsum*Nb_degm2*Area/Vsur

            list_zbin.append(zbin)
            F=Mat_Fisher_1tracer_fnl_bias(n,bg,zbin,p,Vsur,kmin,kmax,mod,cosmo=cosmo)
            list_sigma_fnl.append((np.linalg.inv(F)[0][0])**0.5)
            F1=np.zeros((3,3))
            F1[0,0],F1[1,0],F1[0,1],F1[1,1]=F[0,0],F[1,0],F[0,1],F[1,1]
            Flist.append(F1)
        #print(Flist[-1])
    zeff=(np.sum(zarray*nza)*Na_degm2+np.sum(zarray*nzb)*Nb_degm2)/(Na_degm2+Nb_degm2)
    Flist=np.array(Flist)
    Ftot=np.sum(Flist,axis=0)
    sigma_fnl_eff=(np.linalg.inv(Ftot)[0][0])**0.5
    return list_zbin, list_sigma_fnl, zeff,sigma_fnl_eff
