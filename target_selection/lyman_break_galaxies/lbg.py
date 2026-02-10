import numpy as np

def LBG_SELECTION(tab, conv, name='COSMOS_TMG_U_normagcut'):

    if name=='COSMOS_TMG_U_normagcut': return COSMOS_TMG_U_normagcut(tab, conv)
    if name=='COSMOS_BXU_U_normagcut': return  COSMOS_BXU_U_normagcut(tab, conv)
    if name=='COSMOS_BXU_U_ext_normagcut': return COSMOS_BXU_U_ext_normagcut(tab, conv)
    if name=='XMMLSS_uS_dropout_normagcut': return XMMLSS_uS_dropout_normagcut(tab, conv)
    if name=='COSMOS_G_noimagcut': return COSMOS_G_noimagcut(tab, conv)
    if name=='COSMOS_R_nozmagcut': return COSMOS_R_nozmagcut(tab, conv)

def COSMOS_TMG_U_normagcut(tab, conv):
    mask = tab[conv['u']] - tab[conv['g']] > 0.3
    mask *= (-0.5 < tab[conv['g']] - tab[conv['r']]) * (tab[conv['g']] - tab[conv['r']] < 1.)
    mask_1 = (tab[conv['u']] - tab[conv['g']] > 2.2*(tab[conv['g']] - tab[conv['r']]) + 0.32)
    mask_2 = ((tab[conv['u']] - tab[conv['g']] > 0.9) * (tab[conv['u']] - tab[conv['g']] > 1.6*(tab[conv['g']] - tab[conv['r']]) + 0.75))
    mask *= (mask_1 + mask_2)
    #mask *= (22.5 < tab[conv['r']]) * (tab[conv['r']] < 23.75)
    return mask

def COSMOS_BXU_U_ext_normagcut(tab, conv):
    mask = tab[conv['u']] - tab[conv['g']] > 0
    mask *= (-0.5 < tab[conv['g']] - tab[conv['r']]) * (tab[conv['g']] - tab[conv['r']] < 1.2)
    mask_1 = (tab[conv['u']] - tab[conv['g']]  > 2.2*(tab[conv['g']] - tab[conv['r']]) + 0.32)
    mask_2 = ((tab[conv['u']] - tab[conv['g']] > 0.9) * (tab[conv['u']] - tab[conv['g']] > 1.6*(tab[conv['g']] - tab[conv['r']]) + 0.75))
    mask *= (mask_1 + mask_2)
    return mask

def COSMOS_BXU_U_normagcut(tab, conv):
    mask = tab[conv['u']] - tab[conv['g']] > 0
    mask *= (-0.5 < tab[conv['g']] - tab[conv['r']]) * (tab[conv['g']] - tab[conv['r']] < 1.2)
    mask_1 = (tab[conv['u']] - tab[conv['g']]  > 2.2*(tab[conv['g']] - tab[conv['r']]) + 0.32)
    mask_2 = ((tab[conv['u']] - tab[conv['g']] > 0.9) * (tab[conv['u']] - tab[conv['g']] > 1.6*(tab[conv['g']] - tab[conv['r']]) + 0.75))
    mask *= (mask_1 + mask_2)
    return mask

def COSMOS_G_noimagcut(tab, conv):
    mask = tab[conv['g']] - tab[conv['r']] > 1.
    mask *= (-1.5 < tab[conv['r']] - tab[conv['i']]) * (tab[conv['r']] - tab[conv['i']] < 1.)
    mask *= (tab[conv['g']] - tab[conv['r']] > 1.5*(tab[conv['r']] - tab[conv['i']]) + 0.8) 
    #mask *= (22.5 < tab[conv['i']]) * (tab[conv['i']] < 25.5)
    return mask

def COSMOS_R_nozmagcut(tab, conv):
    mask = tab[conv['r']] - tab[conv['i']] > 1.7
    mask *= tab[conv['r']] - tab[conv['i']] < 5
    mask *= ((tab[conv['r']] - tab[conv['i']]) > 1.6 * (tab[conv['i']] - tab[conv['z']]) + 1.2)
    mask *= tab[conv['i']] - tab[conv['z']] < 1.5
    mask *= tab[conv['i']] - tab[conv['z']] > -0.5
    #mask *= ((tab[conv['r']] - tab[conv['i']]) > 1.5*(tab[conv['i']] - tab[conv['z']]) + 1)
    #mask *= (22.5 < tab[conv['i']]) * (tab[conv['i']] < 25.5)
    return mask
    
def XMMLSS_uS_dropout_normagcut(tab, conv):

    mask = tab[conv['u']] - tab[conv['g']] > 0.3
    mask *= (0 < tab[conv['g']] - tab[conv['r']]) * (tab[conv['g']] - tab[conv['r']] < 0.8)
    
    mask_1 = (tab[conv['u']] - tab[conv['g']] > 2.0*(tab[conv['g']] - tab[conv['r']]) + 0.42)
    mask_2 = ((tab[conv['u']] - tab[conv['g']] > 1.6*(tab[conv['g']] - tab[conv['r']]) + 0.55))
    mask *= (mask_1 + mask_2)
    
    #mask *= (22.7 < tab[conv['r']]) * (tab[conv['r']] < 24.2) 
    return mask