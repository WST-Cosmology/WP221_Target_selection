import numpy as np

def E_wst_udrop_single_exp(m):
    return -0.18*m + 4.8

def E_wst_lbg_dropout(m_i, D_zu, D_zx, E_wst_udrop_single_exp, dropout_band = 'u'):
    """
    Parameters:
        m_i : float
            Input magnitude
        D_zu : float
            Growth factor D(z_u)
        D_zg : float
            Growth factor D(z_g)
        E_wst_udrop : function
            Function that takes magnitude and returns E_WST^udrop

    Returns:
        float : E_WST^drop
    """
    if not dropout_band == 'u':
        ratio = D_zu / D_zx
        shifted_m = m_i -5 * np.log10(ratio)
    else:
        ratio=1
        shifted_m = m_i
    
    return (ratio ** 2) * E_wst_udrop_single_exp(shifted_m)

# define a simple udrop function