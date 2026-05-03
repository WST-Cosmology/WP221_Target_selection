import numpy as np

def E_wst(redshift, mag, tracer = 'BG_faint'):
    if tracer == 'BG_faint': return E_wst_bg_faint(redshift, mag)
    if tracer == 'BG_bright': return E_wst_bg_bright(redshift, mag)
    if tracer == 'ELG': return E_wst_elg(redshift, mag)
    if tracer == 'LRG': return E_wst_lrg(redshift, mag)
    if tracer == 'QSO': return E_wst_qso(redshift, mag)
    if tracer == 'LBGu': return E_wst_lbg_dropout_piecewise(redshift, mag, p_min=0.7, plateau=0.8, dropout_band = 'u', return_magnitudes=False)
    if tracer == 'LBGg': return E_wst_lbg_dropout_piecewise(redshift, mag, p_min=0.5, plateau=0.8, dropout_band = 'g', return_magnitudes=False)
    if tracer == 'LBGr': return E_wst_lbg_dropout_piecewise(redshift, mag, p_min=0.3, plateau=0.8, dropout_band = 'r', return_magnitudes=False)
def n_pass_wst(redshift, mag, tracer = 'BG_faint'):
    if tracer == 'BG_faint': return n_pass_wst_bg_faint(redshift, mag)
    if tracer == 'BG_bright': return n_pass_wst_bg_bright(redshift, mag)
    if tracer == 'ELG': return n_pass_wst_elg(redshift, mag)
    if tracer == 'LRG': return n_pass_wst_lrg(redshift, mag)
    if tracer == 'QSO': return n_pass_wst_qso(redshift, mag)
    if tracer == 'LBGu': return n_pass_wst_lbg_dropout_piecewise(redshift,mag,p_min=0.7,plateau=0.8,dropout_band='u')
    if tracer == 'LBGg': return n_pass_wst_lbg_dropout_piecewise(redshift,mag,p_min=0.5,plateau=0.8,dropout_band='g')
    if tracer == 'LBGr': return n_pass_wst_lbg_dropout_piecewise(redshift,mag,p_min=0.3,plateau=0.8,dropout_band='r')

### BG (bright, faint) ####
def E_wst_bg_bright(redshift, mag):
    return 0.99 * np.ones_like(redshift)
def n_pass_wst_bg_bright(redshift,mag):
    n = np.ones(len(mag)) 
    return n
def E_wst_bg_faint(redshift, mag):
    return 0.99 * np.ones_like(redshift)
def n_pass_wst_bg_faint(redshift,mag):
    n = np.ones(len(mag)) 
    return n

########## ELG ###########
def E_wst_elg(redshift, mag):
    # [O II] doublet observable in DESI roughly for 0.6 < z < 1.6
    return 0.726 * ((redshift > 0.6) & (redshift < 1.6))
def n_pass_wst_elg(redshift,mag):
    n = np.ones(len(mag)) 
    return n

########## LRG ###########
def E_wst_lrg(redshift, mag):
    # DESI LRG spectroscopic efficiency roughly valid for 0.4 < z < 1.1
    return 0.99 * ((redshift > 0.4) & (redshift < 1.1))
def n_pass_wst_lrg(redshift,mag):
    n = np.ones(len(mag)) 
    return n

########## QSO ###########
def E_wst_qso(redshift, mag):
    # Mg II enters DESI range at z ~ 0.3
    return 0.70 * (redshift > 0.3)
def n_pass_wst_qso(redshift,mag):
    n = np.ones(len(mag)) 
    return n

########## LBG ##########
def E_mse_udrop_single_exp(m): 
    return np.maximum((-0.18*m + 4.8), 0)
def E_desi_udrop_single_exp(m): 
    return np.maximum(-0.2*(m-23.5) + 0.75, 0)
    
def E_wst_lbg_dropout(redshift, mag, dropout_band = 'u', reference = 'mse'):

    Dz_ugr = [6, 7, 8] # Mpc

    texp = 1000 #Dark time WST
    SWST=12**2 #Surface of WST mirror
    if reference=='mse':
        tMSE = 1800 #MSE time
        SMSE = 11.25**2
        alpha = np.sqrt(texp * SWST)/np.sqrt(tMSE * SMSE)
    if reference=='desi':
        tDESI= 2*60*60 #Dark time WST
        SDESI = 3.8**2 #Surface of MSE mirror
        alpha = np.sqrt(texp * SWST)/np.sqrt(tDESI * SDESI)
    
    def f(z, k=10): return 1 / (1 + np.exp(-k * (z - 2.5)))
        
    D_zu = Dz_ugr[0]
    if dropout_band == 'u': D_zx = Dz_ugr[0]
    if dropout_band == 'g': D_zx = Dz_ugr[1]
    if dropout_band == 'r': D_zx = Dz_ugr[2]

    ratio = D_zu / D_zx
    shifted_m = mag -5 * np.log10(ratio)

    efficency = (ratio ** 2) * alpha * E_mse_udrop_single_exp(shifted_m) 
    #if dropout_band == 'u':
    efficency *= f(redshift)
    return efficency

def E_wst_lbg_dropout_piecewise(redshift, mag, p_min=0.7, plateau=0.8, dropout_band = 'u', return_magnitudes=False):

    E1 = np.minimum(E_wst_lbg_dropout(redshift, mag, dropout_band=dropout_band), plateau)
    E2 = E1 + np.minimum(np.sqrt(2)*E1, plateau)*(1-E1)
    E3 = E2 + np.minimum(np.sqrt(3)*E1, plateau)*(1-E2)

    m1 = mag[np.argmin(np.abs(E1 - p_min))]
    m2 = mag[np.argmin(np.abs(E2 - p_min))]

    E = np.where(mag < m1, E1, np.where(mag < m2, E2, E3))
    if return_magnitudes:
        return E, m1, m2
    else: return E

def n_pass_wst_lbg_dropout_piecewise(redshift,mag,p_min=0.7,plateau=0.8,dropout_band='u'):

    # cumulative efficiencies
    E1_raw = E_wst_lbg_dropout(redshift, mag, dropout_band=dropout_band)
    E1 = np.minimum(np.clip(E1_raw, 0, 1), plateau)

    E2 = E1 + np.minimum(np.sqrt(2) * E1, plateau) * (1 - E1)
    E3 = E2 + np.minimum(np.sqrt(3) * E1, plateau) * (1 - E2)

    # transition magnitudes
    m1 = mag[np.argmin(np.abs(E1 - p_min))]
    m2 = mag[np.argmin(np.abs(E2 - p_min))]
    m3 = 30

    # conditional probabilities
    p2 = np.zeros_like(E1)
    mask = (1 - E1) > 0
    p2[mask] = (E2[mask] - E1[mask]) / (1 - E1[mask])

    # expected number of passes
    n = np.zeros(len(mag))

    # 1-pass regime
    mask1 = mag < m1
    n[mask1] = 1

    # 2-pass regime
    mask2 = (mag >= m1) & (mag < m2)
    n[mask2] = (1 + (1 - E1[mask2]))

    # 3-pass regime
    mask3 = (mag >= m2) & (mag < m3)
    n[mask3] = (
        1
        + (1 - E1[mask3])
        + (1 - E1[mask3]) * (1 - p2[mask3])
    )

    return n

def E_wst_lbg_dropout_piecewise_wdassignies(redshift, mag, p_min=0.7, plateau=0.8, dropout_band = 'u', return_magnitudes=False):

    E1 = np.minimum(E_wst_lbg_dropout(redshift, mag, dropout_band=dropout_band), plateau)
    E2 = E1 + np.minimum(np.sqrt(2)*E1, plateau)*(1-E1)
    E3 = E2 + np.minimum(np.sqrt(3)*E1, plateau)*(1-E2)

    m1 = mag[np.argmin(np.abs(E1 - p_min))]
    m2 = mag[np.argmin(np.abs(E2 - p_min))]
    m3 = mag[np.argmin(np.abs(E3 - p_min))]

    E = np.zeros(len(mag))
    E[mag <= m1] = E1[mag <= m1]
    E[(mag > m1)*(mag <= m2)] = p_min
    E[(mag > m2)*(mag <= m3)] = p_min
    E[(mag > m3)] = E3[(mag > m3)]
    return E
    
def n_pass_wst_lbg_dropout_piecewise_wdassignies(redshift, mag, p_min=0.7, plateau=0.8, dropout_band='u'):

    E1 = np.minimum(E_wst_lbg_dropout(redshift, mag, dropout_band=dropout_band), plateau)
    E2 = E1 + np.minimum(np.sqrt(2)*E1, plateau)*(1-E1)
    E3 = E2 + np.minimum(np.sqrt(3)*E1, plateau)*(1-E2)

    m1 = mag[np.argmin(np.abs(E1 - p_min))]
    m2 = mag[np.argmin(np.abs(E2 - p_min))]
    m3 = mag[np.argmin(np.abs(E3 - p_min))]

    eta_2 = np.clip((p_min - E1)/(np.sqrt(2)*E1), 0, 1)
    eta_3 = np.clip((p_min - E1 * (np.sqrt(2) + 1 - np.sqrt(2) * E1))/(np.sqrt(3)*E1), 0, 1)

    n = np.zeros(len(E1))

    mask_m1 = mag <= m1
    n[mask_m1] = 1  

    mask_m1_m2 = (mag > m1) * (mag <= m2)
    n[mask_m1_m2] = 1 + eta_2[mask_m1_m2]

    mask_m2_m3 = (mag > m2) * (mag <= m3)
    n[mask_m2_m3] = 2 - E1[mask_m2_m3] + eta_3[mask_m2_m3]

    mask_m3 = mag > m3
    n[mask_m3] = 3 - E3[mask_m3]

    return n


#def n_pass_wst_lbg_dropout_piecewise_old(redshift, mag, p_min=0.7, plateau=0.8, dropout_band = 'u'): 
#    E, m1, m2 = E_wst_lbg_dropout_piecewise(redshift, mag, p_min=p_min, plateau=plateau, dropout_band = dropout_band, return_magnitudes=True)
#    m3=30 
#    E1 = E_wst_lbg_dropout(redshift, mag, dropout_band=dropout_band) 
#    E1_plateau = np.minimum(E1, plateau) 
#    E1sqrt2_plateau = np.minimum(np.sqrt(2)*E1, plateau) 
#    n = np.zeros(len(mag)) 
#    n[mag < m1] = 1 
#@    n[(mag >= m1)*(mag < m2)] = 1 + (1 - E1_plateau[(mag >= m1)*(mag < m2)]) 
#    n[(mag >= m2)*(mag < m3)] = 1 + (1 - E1_plateau[(mag >= m2)*(mag < m3)]) + (1 - E1_plateau[(mag >= m2)*(mag < m3)])*(1 - E1sqrt2_plateau[(mag >= m2)*(mag < m3)]) 
#    return n
#def n_pass_wst_lbg_dropout_piecewise_wdassignies_old(redshift, mag, p_min=0.7, plateau=0.8, dropout_band = 'u'):

#    E, m1, m2 = E_wst_lbg_dropout_piecewise_wdassignies(redshift, mag, p_min=p_min, plateau=plateau, dropout_band = dropout_band, return_magnitudes=True)
#    E1 = np.minimum(E_wst_lbg_dropout(redshift, mag, dropout_band=dropout_band), plateau)

#    eta_2 = np.clip((p_min - E1)/(np.sqrt(2)*E1), 0, 1)
#    eta_3 = np.clip((p_min - E1 * (np.sqrt(2) + 1 - np.sqrt(2) * E1))/(np.sqrt(3)*E1), 0,1)
#    n = np.zeros(len(E1))
#    mask_m1 = mag <= m1
#    n[mask_m1] = 1  
#    mask_m1_m2 = (mag > m1)*(mag <= m2)
#    n[mask_m1_m2] = 1  + eta_2[mask_m1_m2]
#    mask_m2 = (mag > m2)
#    n[mask_m2] = 2 - E1[mask_m2] + eta_3[mask_m2]
#    return n