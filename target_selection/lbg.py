import numpy as np

def refined_selection_fct(tab, conv):
    
    mask  = tab[conv['u_err']] > 0
    mask *= tab[conv['g_err']] > 0
    mask *= tab[conv['r_err']] > 0

    mask *= tab[conv['u_err']] < 0.5
    mask *= tab[conv['g_err']] < 0.5
    mask *= tab[conv['r_err']] < 0.5

    mask *= abs(tab[conv['u']] - tab[conv['g']])  < 10
    mask *= abs(tab[conv['g']] - tab[conv['r']])  < 10
    mask *= abs(tab[conv['u']] - tab[conv['r']])  < 10
    return mask

def LBG_SELECTION(tab, conv, name='COSMOS_TMG_U'):

    if name=='COSMOS_TMG_U': return COSMOS_TMG_U(tab, conv)
    if name=='COSMOS_BXU_U': return COSMOS_BXU_U(tab, conv)
    if name=='COSMOS_BXU_U_ext': return COSMOS_BXU_U_ext(tab, conv)
    if name=='COSMOS_BXU_U_ext_normagcut': return COSMOS_BXU_U_ext_normagcut(tab, conv)
    if name=='COSMOS_G': return COSMOS_G(tab, conv)

def COSMOS_TMG_U_ext(tab, conv):

    mask = tab[conv['u']] - tab[conv['g']] > 0.3
    mask *= (-0.5 < tab[conv['g']] - tab[conv['r']]) * (tab[conv['g']] - tab[conv['r']] < 1.)
    mask_1 = (tab[conv['u']] - tab[conv['g']] > 2.2*(tab[conv['g']] - tab[conv['r']]) + 0.32)
    mask_2 = ((tab[conv['u']] - tab[conv['g']] > 0.9) * (tab[conv['u']] - tab[conv['g']] > 1.6*(tab[conv['g']] - tab[conv['r']]) + 0.75))
    mask *= (mask_1 + mask_2)
    mask *= (23.2 < tab[conv['r']]) * (tab[conv['r']] < 24.2) 
    return mask

def COSMOS_TMG_U(tab, conv):
    mask = tab[conv['u']] - tab[conv['g']] > 0.3
    mask *= (-0.5 < tab[conv['g']] - tab[conv['r']]) * (tab[conv['g']] - tab[conv['r']] < 1.)
    mask_1 = (tab[conv['u']] - tab[conv['g']] > 2.2*(tab[conv['g']] - tab[conv['r']]) + 0.32)
    mask_2 = ((tab[conv['u']] - tab[conv['g']] > 0.9) * (tab[conv['u']] - tab[conv['g']] > 1.6*(tab[conv['g']] - tab[conv['r']]) + 0.75))
    mask *= (mask_1 + mask_2)
    mask *= (22.5 < tab[conv['r']]) * (tab[conv['r']] < 23.75)
    return mask

def COSMOS_BXU_U_custom(tab, conv):
    mask = tab[conv['u']] - tab[conv['g']] > 0
    mask *= (-0.5 < tab[conv['g']] - tab[conv['r']]) * (tab[conv['g']] - tab[conv['r']] < 1.2)
    mask_1 = (tab[conv['u']] - tab[conv['g']]  > 2.2*(tab[conv['g']] - tab[conv['r']]) + 0.32)
    mask_2 = ((tab[conv['u']] - tab[conv['g']] > 0.9) * (tab[conv['u']] - tab[conv['g']] > 1.6*(tab[conv['g']] - tab[conv['r']]) + 0.75))
    mask *= (mask_1 + mask_2)
    return mask

def COSMOS_BXU_U_ext(tab, conv):
    mask = tab[conv['u']] - tab[conv['g']] > 0
    mask *= (-0.5 < tab[conv['g']] - tab[conv['r']]) * (tab[conv['g']] - tab[conv['r']] < 1.2)
    mask_1 = (tab[conv['u']] - tab[conv['g']]  > 2.2*(tab[conv['g']] - tab[conv['r']]) + 0.32)
    mask_2 = ((tab[conv['u']] - tab[conv['g']] > 0.9) * (tab[conv['u']] - tab[conv['g']] > 1.6*(tab[conv['g']] - tab[conv['r']]) + 0.75))
    mask *= (mask_1 + mask_2)
    mask *= (23.75 < tab[conv['r']]) * (tab[conv['r']] < 24.2)
    return mask

def COSMOS_BXU_U_ext_normagcut(tab, conv):
    mask = tab[conv['u']] - tab[conv['g']] > 0
    mask *= (-0.5 < tab[conv['g']] - tab[conv['r']]) * (tab[conv['g']] - tab[conv['r']] < 1.2)
    mask_1 = (tab[conv['u']] - tab[conv['g']]  > 2.2*(tab[conv['g']] - tab[conv['r']]) + 0.32)
    mask_2 = ((tab[conv['u']] - tab[conv['g']] > 0.9) * (tab[conv['u']] - tab[conv['g']] > 1.6*(tab[conv['g']] - tab[conv['r']]) + 0.75))
    mask *= (mask_1 + mask_2)
    return mask

def COSMOS_BXU_U(tab, conv):
    mask = tab[conv['u']] - tab[conv['g']] > 0
    mask *= (-0.5 < tab[conv['g']] - tab[conv['r']]) * (tab[conv['g']] - tab[conv['r']] < 1.2)
    mask_1 = (tab[conv['u']] - tab[conv['g']]  > 2.2*(tab[conv['g']] - tab[conv['r']]) + 0.32)
    mask_2 = ((tab[conv['u']] - tab[conv['g']] > 0.9) * (tab[conv['u']] - tab[conv['g']] > 1.6*(tab[conv['g']] - tab[conv['r']]) + 0.75))
    mask *= (mask_1 + mask_2)
    mask *= (23.75 < tab[conv['r']]) * (tab[conv['r']] < 24.5)
    return mask

def COSMOS_G(tab, conv):
    mask = tab[conv['g']] - tab[conv['r']] > 1.
    mask *= (-1.5 < tab[conv['r']] - tab[conv['i']]) * (tab[conv['r']] - tab[conv['i']] < 1.)
    mask *= (tab[conv['g']] - tab[conv['r']] > 1.5*(tab[conv['r']] - tab[conv['i']]) + 0.8) 
    mask *= (22.5 < tab[conv['i']]) * (tab[conv['i']] < 25.5)
    return mask