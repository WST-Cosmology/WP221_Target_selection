""
'''
Package
'''
#Computing
from astropy.io import fits
import numpy as np
import math
import treecorr
import time
import random
import pyccl as ccl
import gc
from scipy.optimize import least_squares
from scipy.integrate import quad, dblquad
from scipy import interpolate
#for saving:
import pickle

#for ploting
import matplotlib.pyplot as plt
from matplotlib import rcParams
from IPython.display import display, Math

""
"""
Parameters
"""

import sys
sample=sys.argv[1]

#sample = 'magmax'

scale_def=['theta',5,20]
z_range=[0.,1.6]
Nz=16
dz=0.1

H_array=np.array([20.35,20.85,21.35,21.7])
Zarray=np.array([z_range[0]+dz*(i+0.5) for i in range(Nz)]) #Z list for 8 bins
ialpha=1 # weithing equal accross scales



""
"""
Cosmology
"""
yourz=0.
h=0.67
H0=100*h
Omega_m=0.315
c_ls=299.792*10**3
cosmo = ccl.Cosmology(Omega_c=0.27, Omega_b=Omega_m-0.27, h=h, A_s=2.1e-9, n_s=0.96)
a=1./(1+yourz)

def ra_cosmo(z):
    return(ccl.comoving_angular_distance(cosmo, 1/(1+z)))

def D_cosmo(z):
    return ccl.growth_factor(cosmo,1/(1+z))

def chi(z):
    return(ccl.comoving_radial_distance(cosmo, 1/(1+z)))

def chi_3d(r,z):
    return ccl.correlations.correlation_3d(cosmo,1/(1+z),r,'delta_matter:delta_matter')

def PNL(l,z):
    #return ccl.power.nonlin_matter_power(cosmo, k=(l+0.5)/chi(z), a=1/(1+z))
    return ccl.power.nonlin_power(cosmo, k=(l+0.5)/chi(z), a=1/(1+z), p_of_k_a='delta_matter:delta_matter')

def Plin(l,z):
    #return ccl.power.linear_matter_power(cosmo, k=(l+0.5)/chi(z), a=1/(1+z))
    return ccl.power.linear_power(cosmo, k=(l+0.5)/chi(z), a=1/(1+z), p_of_k_a='delta_matter:delta_matter')

""
"""
Compute wm(z)
"""


def xi_dm_theta(theta,z, typ):
    '''same as before for theta in deg instead of rp'''
    Hz=100*h*ccl.background.h_over_h0(cosmo,a=1/(1+z))

    if typ=='linear':
        Plin_ell = [Plin(l,z) for l in Ell]
        return Hz/c_ls/chi(z)**2*ccl.correlations.correlation(cosmo,ell=Ell,C_ell=Plin_ell,theta=theta,type='NN',method='Legendre')
    
    if typ=='NL':
        PNL_ell = [PNL(l,z) for l in Ell]
        return Hz/c_ls/chi(z)**2*ccl.correlations.correlation(cosmo,ell=Ell,C_ell=PNL_ell,theta=theta,type='NN',method='Legendre')

def w_dm_theta(z,typ):
    '''same as before for theta  in deg instead of rp'''
    Hz=100*h*ccl.background.h_over_h0(cosmo,a=1/(1+z))
    if typ=='linear':
        P_delta = [Plin(l,z) for l in Ell]
    if typ=='NL':
        P_delta = [PNL(l,z) for l in Ell]    
    
    #weight that should be normalised to unity
    def Wtheta(theta):
        return theta**alpha
    norm=quad(Wtheta,theta_min, theta_max,epsrel=10**(-6),epsabs=10**(-6),limit=nlim)[0]
        
    def xi_dm_int(theta):
        return Hz/c_ls/chi(z)**2*ccl.correlations.correlation(cosmo,ell=Ell,C_ell=P_delta,theta=theta,type='NN',method='Legendre')

    def Integ(theta):
        return Wtheta(theta)/norm*xi_dm_int(theta)
    
    return quad(Integ,theta_min, theta_max,epsrel=10**(-6),epsabs=10**(-6),limit=nlim)[0]


Ell=range(1,4000)
nlim=10000

theta_min=scale_def[1]/60 #deg
theta_max=scale_def[2]/60 #deg

if ialpha==0:
    alpha=-1
elif ialpha==1:
    alpha=0
if ialpha==2:
    alpha=1
    

#W_dm_lin_theta=[w_dm_theta(z,typ='linear') for z in Zarray]
W_dm_NL_theta=[w_dm_theta(z,typ='NL') for z in Zarray]

w_m=W_dm_NL_theta

print('wm=',w_m)
""
"""
# Load data
"""
#path for the data files
path = '/data/astro/scratch/wdassign/WST/'

def load_catalog(filename,sample):
    """Load a FITS catalog and return relevant columns and weights."""
    with fits.open(filename) as hdul:
        data = hdul[1].data
        ra = data['ra_gal']
        dec = data['dec_gal']
        z = data['true_redshift_gal']
        if sample=='magmax':
            magn=-2.5*np.log10(data['euclid_nisp_h'])-48.6
        elif sample=='BG':
            magn=-2.5*np.log10(data['lsst_r'])-48.6
        elif sample=='LRG':
            magn=-2.5*np.log10(data['lsst_z'])-48.6
        elif sample=='ELG':
            magn=-2.5*np.log10(data['lsst_g'])-48.6
        zp=data['zp']
    return ra, dec, z, magn,zp


if sample=='magmax':
    data_file = path + 'FS2MagMax.fits'
elif sample=='BG':
    data_file = path + 'BG.fits'
elif sample=='LRG':
   data_file = path + 'LRG.fits'
elif sample=='ELG':
   data_file = path + 'ELG.fits'

ra_sample, dec_sample, z_sample, mag_sample,zp_sample = load_catalog(data_file,sample)

ra_min = np.min(ra_sample)
ra_max = np.max(ra_sample)
dec_min = np.min(dec_sample)
dec_max = 20#np.max(dec_sample)

print('ra range = %f .. %f' % (ra_min, ra_max))
print('dec range = %f .. %f' % (dec_min, dec_max))

sky_deg2=41256/(4*math.pi)*math.pi/180*(ra_max-ra_min)*(math.cos(math.pi/2-math.pi/180*dec_max)-math.cos(math.pi/2-math.pi/180*dec_min))
print(' we select', sky_deg2,' deg2')


sel_sky=((ra_sample>=ra_min)&(ra_sample<ra_max)&(dec_sample>=dec_min)&(dec_sample<dec_max))
ra_sample, dec_sample, z_sample, mag_sample,zp_sample=ra_sample[sel_sky], dec_sample[sel_sky], z_sample[sel_sky], mag_sample[sel_sky],zp_sample[sel_sky]


## Create randoms
ra_rand_sample = np.random.uniform(ra_min, ra_max, len(ra_sample)//2)
rand_sindec = np.random.uniform(np.sin(dec_min*np.pi/180), np.sin(dec_max*np.pi/180), len(ra_sample)//2)
dec_rand_sample = np.arcsin(rand_sindec)*180/np.pi


""
"""
Measure auto corr
"""


def w_xx(sample_x, scale_range, z_range, Ntheta, Nz, Npatch, z_for_rand ,Eta_rand, option):
    '''
    Code to compute the auto-correlation of sample_x,
    -------------------------------------
    Outputs:
    Wz: is a Nz x Nalpha array, where Nz is the number of z bins, and Nalpha the number of scale weighting
    Thus Wz[iz][ialpha] is the value for wxx for the z-bin iz, and the scale weighting ialpha
    
    Err is a Nz x Nalpha uncertainty array associated to WZ.

    Multi_cov is the Nz x Njkk x Nalpha array usefull to compute the full covariance!!

    For std tests, use directly WZ and Err. Usually redshift covariance is negligible.
    -------------------------------------
    Inputs:
    sample_x: name of the sample, the only option now is sample_x='photo' and 'eBOSS ELG'
    
    scale_range is an array = [convention,scale_min,scale_max], with convention ='rp' or 'theta', with scale_min/max IN ARCMIN or MPC
    z_for_rand is an array = [zmin,zmax]

    Ntheta is the number of theta-bins (default use 10)
    Nz is the number of z-bins
    Npatch the number of Jkk patches to evaluate the cov
    
    z_for_rand==True means there are redshifts for randoms
    Eta_rand is how many more randoms you want, to evaluate the DD/RR counts
    -------------------------------------
    '''
    print('evaluate the auto-corr of ', sample_x)
    print('for',z_range[0], '<z<',z_range[1],'  with ', Nz, ' bins')
    if len(option)>0:
        print('for ',option[0], option[1],'<',option[2])
    zmin=z_range[0]
    zmax=z_range[1]
    dz=(zmax-zmin)/Nz
    
    # You can add other options, with another elif
    if sample_x=='eBOSS ELG':
        ra_gal =ra_eboss_ELG_south
        dec_gal =dec_eboss_ELG_south
        z_gal =z_eboss_ELG_south
        w_gal =  weight_eboss_ELG_south
        
        ra_rand =ra_rand_eboss_ELG_south
        dec_rand =dec_rand_eboss_ELG_south
        w_rand =  weight_rand_eboss_ELG_south
        if z_for_rand==True:
            z_rand=z_rand_eboss_ELG_south
    elif sample_x=='magmax':
        Hmin=option[1]
        Hmax=option[2]
        sel_H=((mag_sample>=Hmin)&(mag_sample<Hmax))
        ra_gal =ra_sample[sel_H]
        dec_gal =dec_sample[sel_H]
        z_gal =z_sample[sel_H]
        w_gal = np.ones(len(ra_gal)) #if you do have weights, replace with np.ones(len(ra_gal))
        
        ra_rand =ra_rand_sample
        dec_rand =dec_rand_sample
        w_rand =  np.ones(len(ra_rand))
        if z_for_rand==True:
            z_rand=1
    else:
        print('sample_x not included in the code')

    print(len(ra_gal) ,' galaxies ')
    print(len(ra_rand) ,'randoms ')

    
    Wz=[]
    Cov=[]
    Err=[]
    Multi_cov=[]
    
    #create a catalog to have the Jkk patches
    cat_patch = treecorr.Catalog(ra=ra_gal,dec=dec_gal,w=w_gal,ra_units='degrees',dec_units='degrees',npatch=Npatch)

    def integ_wxx(w1):
        '''
        integrate w1 over theta, with different scales weighting, defined by List_alpha:
        W(theta)=theta**alpha/norm
        '''
        Results_alpha=[]
        for alpha in  [-1,0,1]:
            w1b=np.sum(w1*rlist**alpha*redges)
            norm=np.sum(rlist**alpha*redges)
            Results_alpha.append(w1b/norm)
        return np.array(Results_alpha)
        
    for iz in range(Nz):
        zmean_i=zmin+(iz+0.5)*dz
        zmin_i=zmin+(iz)*dz
        zmax_i=zmin+(iz+1)*dz

        if scale_range[0]=='theta':
            theta_min=scale_range[1]/60 # arcmin to degree
            theta_max=scale_range[2]/60 # arcmin to degree
        elif scale_range[0]=='rp':
            theta_min=scale_range[1]/ra_cosmo(zmean_i)*360/(2*math.pi) #from Mpc to rad, to deg
            theta_max=scale_range[2]/ra_cosmo(zmean_i)*360/(2*math.pi) #from Mpc to rad, to deg
        else:
            print('error scale conv')
            return 0,0,0

            
        print('zi=',round(zmean_i,3),' dz=',round(dz,3))
        
        sel_gal_subbin=((z_gal>=zmin_i)&(z_gal<zmax_i))
        cat_gal_subbin=treecorr.Catalog(ra=ra_gal[sel_gal_subbin],dec=dec_gal[sel_gal_subbin],w=w_gal[sel_gal_subbin],ra_units='degrees',dec_units='degrees',patch_centers=cat_patch.patch_centers)
        

        # Do we have z for randoms: 
        if z_for_rand==True:
            sel_rand_subbin=((z_rand>=zmin_i)&(z_rand<zmax_i))
        else:
            sel_rand_subbin=np.ones(len(ra_rand),dtype=bool)
        # Select Eta_rand more randoms than gal
        Ntot=len(ra_rand[sel_rand_subbin])
        Nrand=Eta_rand*len(ra_gal[sel_gal_subbin])
        Index=random.choices(range(Ntot), k=Nrand)
        
        
        ra_rand_select=ra_rand[sel_rand_subbin][Index]
        dec_rand_select=dec_rand[sel_rand_subbin][Index]
        w_rand_select=w_rand[sel_rand_subbin][Index]
    
        cat_rand_subbin=treecorr.Catalog(ra=ra_rand_select,dec=dec_rand_select,w=w_rand_select,ra_units='degrees',dec_units='degrees', patch_centers=cat_patch.patch_centers)
    
        wxx  = treecorr.NNCorrelation(min_sep=theta_min,max_sep=theta_max,nbins=Ntheta,var_method='jackknife',sep_units='degree',bin_slop=0.01,cross_patch_weight='match')
        rrxx = treecorr.NNCorrelation(min_sep=theta_min,max_sep=theta_max,nbins=Ntheta,var_method='jackknife',sep_units='degree',bin_slop=0.01,cross_patch_weight='match')
        drxx = treecorr.NNCorrelation(min_sep=theta_min,max_sep=theta_max,nbins=Ntheta,var_method='jackknife',sep_units='degree',bin_slop=0.01,cross_patch_weight='match')

        #Now I use r as a name instead of theta
        rlist=wxx.rnom
        redges=wxx.right_edges-wxx.left_edges
        #print(rlist)
        #print(redges)
            
            
        rrxx.process(cat_rand_subbin,cat_rand_subbin)
        drxx.process(cat_rand_subbin,cat_gal_subbin)
        
        wxx.process(cat_gal_subbin,cat_gal_subbin)
        wxx.calculateXi(rr=rrxx,dr=drxx)
         
        my_funct = lambda corrs: integ_wxx(corrs[0].xi)
        corrs = [wxx]
            
        ratio = my_funct(corrs)  
        cov = treecorr.estimate_multi_cov(corrs, 'jackknife', func=my_funct)
        multi_cov=treecorr.build_multi_cov_design_matrix(corrs,'jackknife', func=my_funct, comm=None)
        
        Wz.append(ratio)
        Cov.append(cov)
        Err.append([np.sqrt(cov[i][i]) for i in range(np.size(cov,0))])
        Multi_cov.append(multi_cov[0])
        
        del(cat_gal_subbin)
        del(cat_rand_subbin)
        gc.collect()
    return(np.array(Wz),np.array(Err),np.array(Multi_cov))



A0,B0,C0=w_xx(sample_x='magmax', scale_range=scale_def, z_range=z_range, Ntheta=10, Nz=Nz, Npatch=100, z_for_rand=False ,Eta_rand=5,option=['Hcut',20,20.5])
A1,B1,C1=w_xx(sample_x='magmax', scale_range=scale_def, z_range=z_range, Ntheta=10, Nz=Nz, Npatch=100, z_for_rand=False ,Eta_rand=5,option=['Hcut',20.5,21])
A2,B2,C2=w_xx(sample_x='magmax', scale_range=scale_def, z_range=z_range, Ntheta=10, Nz=Nz, Npatch=100, z_for_rand=False ,Eta_rand=5,option=['Hcut',21,21.5])
A3,B3,C3=w_xx(sample_x='magmax', scale_range=scale_def, z_range=z_range, Ntheta=10, Nz=Nz, Npatch=100, z_for_rand=False ,Eta_rand=5,option=['Hcut',21.5,21.8])


""
"""
Fit bias
"""



def model(params, z_array, m_array):
    a, b, c, d = params

    # Create Nz x Nm grids
    z, m = np.meshgrid(z_array, m_array, indexing='ij')

    return a + (b + c * m) * z**d


def residuals(params, z_array, m_array, b_array):
    return (model(params, z_array, m_array) - b_array).ravel()


# --------------------------------------------------
# Your data
# --------------------------------------------------
# z_array: shape (Nz,)
# m_array: shape (Nm,)
# b_array: shape (Nz, Nm)

# Example initial guess
p0 = [0.0, 1.0, 1.0, 1.0]


b_array=np.transpose(np.array([[np.sqrt(A0[i][ialpha]/w_m[i]*dz) for i in range(Nz)],
                  [np.sqrt(A1[i][ialpha]/w_m[i]*dz) for i in range(Nz)],
                  [np.sqrt(A2[i][ialpha]/w_m[i]*dz) for i in range(Nz)],
                  [np.sqrt(A3[i][ialpha]/w_m[i]*dz) for i in range(Nz)]]))

b_array_err=np.transpose(np.array([[0.5 * np.sqrt(A0[i][ialpha]/w_m[i]*0.1) * B0[i][ialpha] / A0[i][ialpha] for i in range(Nz)],
                                    [0.5 * np.sqrt(A1[i][ialpha]/w_m[i]*0.1) * B1[i][ialpha] / A1[i][ialpha] for i in range(Nz)],
                                    [0.5 * np.sqrt(A2[i][ialpha]/w_m[i]*0.1) * B2[i][ialpha] / A2[i][ialpha] for i in range(Nz)],
                                    [0.5 * np.sqrt(A3[i][ialpha]/w_m[i]*0.1) * B3[i][ialpha] / A3[i][ialpha] for i in range(Nz)]]))


z_array, m_array=Zarray, H_array

# Fit
result = least_squares(
    residuals,
    p0,
    args=(z_array, m_array, b_array)
)

# Best-fit parameters
a_fit, b_fit, c_fit, d_fit = result.x

print("Best-fit parameters:")
print(f"a = {a_fit}")
print(f"b = {b_fit}")
print(f"c = {c_fit}")
print(f"d = {d_fit}")

# Fitted array
b_fit_array = model(result.x, z_array, m_array)

# RMS residual
rms = np.sqrt(np.mean((b_fit_array - b_array)**2))
print(f"RMS residual = {rms}")


""
"""
Save data
"""


data = {
    '20_20p5_data': A0,
    '20_20p5_cov': B0,
    '20_20p5_multi': C0,
    '20p5_21_data': A1,
    '20p5_21_cov': B1,
    '20p5_21_multi': C1,
    '21_21p5_data': A2,
    '21_21p5_cov': B2,
    '21_21p5_multi': C2,
    '21p5_21p8_data': A3,
    '21p5_21p8_cov': B3,
    '21p5_21p8_multi': C3,
    'best_fit':[a_fit,b_fit,c_fit,d_fit]
    
}

chemin_fichier_pickle = "/nfs/pic.es/user/w/wdassign/WST/MagMax/FS2/Data_magmax_ls.pickle"
                       
try:
    # Ouvrir le fichier en mode binaire pour l'écriture
    with open(chemin_fichier_pickle, 'wb') as fichier_pickle:
        # Écrire l'array dans le fichier .pickle en utilisant pickle.dump()
        pickle.dump(data, fichier_pickle)

    print("The data was written with sucess")
except Exception as e:
    print("Error while writting the .pickle :", str(e))




""
"""
Plot results
"""


plt.style.use('/nfs/pic.es/user/w/wdassign/codes_FS2/niceplots/niceplots/euclid_stylesheet_v1.mplstyle')

plt.rcParams['text.usetex'] = False
plt.rcParams['mathtext.fontset'] = 'cm' # Font to use inside mathtext expresions
plt.rcParams["axes.formatter.use_mathtext"] = True # Use mathtext in the ticks of the axes

#rcParams['font.family'] = 'serif'
#rcParams['font.sans-serif'] = ['century']
# Créez une figure


fw, fh=plt.rcParams['figure.figsize'] # Get figure width and height

fig = plt.figure(constrained_layout=True, figsize=(fw*0.9, fh*0.9))
plt.tick_params(direction="in", top=True, right=True, zorder=10)



magn=20.25
plt.errorbar(Zarray, [np.sqrt(A0[i][ialpha]/w_m[i]*0.1) for i in range(Nz)],yerr=[0.5 * np.sqrt(A0[i][ialpha]/w_m[i]*0.1) * B0[i][ialpha] / A0[i][ialpha] for i in range(Nz)],color='darkorange',marker='o',label=r'$20.0<m_H<20.5$',linestyle='')
plt.plot(Zarray, a_fit+(b_fit+c_fit*magn)*Zarray**d_fit,color='darkorange',linestyle='--')

magn=20.75
plt.errorbar(Zarray, [np.sqrt(A1[i][ialpha]/w_m[i]*0.1) for i in range(Nz)],yerr=[0.5 * np.sqrt(A1[i][ialpha]/w_m[i]*0.1) * B1[i][ialpha] / A1[i][ialpha] for i in range(Nz)],color='red',marker='o',label=r'$20.5<m_H<21.0$',linestyle='')
plt.plot(Zarray, a_fit+(b_fit+c_fit*magn)*Zarray**d_fit,color='red',linestyle='--')

magn=21.25
plt.errorbar(Zarray, [np.sqrt(A2[i][ialpha]/w_m[i]*0.1) for i in range(Nz)],yerr=[0.5 * np.sqrt(A2[i][ialpha]/w_m[i]*0.1) * B2[i][ialpha] / A2[i][ialpha] for i in range(Nz)],color='blue',marker='o',label=r'$21.0<m_H<21.5$',linestyle='')
plt.plot(Zarray, a_fit+(b_fit+c_fit*magn)*Zarray**d_fit,color='blue',linestyle='--')

magn=21.65
plt.errorbar(Zarray, [np.sqrt(A3[i][ialpha]/w_m[i]*0.1) for i in range(Nz)],yerr=[0.5 * np.sqrt(A3[i][ialpha]/w_m[i]*0.1) * B3[i][ialpha] / A3[i][ialpha] for i in range(Nz)],color='purple',marker='o',label=r'$21.5<m_H<21.8$',linestyle='')
plt.plot(Zarray, a_fit+(b_fit+c_fit*magn)*Zarray**d_fit,color='purple',linestyle='--')

b_amp=0.84
bgrowth_model=(b_amp*D_cosmo(0)/D_cosmo(Zarray))
plt.plot(Zarray,bgrowth_model,color='dimgrey',linestyle='-.',label=r'$b(z)=0.84/D(z)$ (ELG)')



plt.plot([],[],color='black',linestyle='--',label=r'best fit $b(z,m)$')
plt.xlabel('Redshift $z$')
plt.legend()

plt.title('MagMax bias')
plt.ylabel(' Galaxy bias $b$')
plt.savefig('/nfs/pic.es/user/w/wdassign/WST/MagMax/FS2/Bias_magmax_1pannel_ls.png')
