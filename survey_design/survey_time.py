import numpy as np

def passes_needed(Xtarget_init, Xfiber, X_I_want, max_passes=2000, coll=True, poiss=True):
    nleft = []
    Xremaining = Xtarget_init
    cumulative = 0.0
    nleft.append(Xtarget_init)
    for n in range(1, max_passes + 1):
        lambda_ = Xremaining/Xfiber
        p_associated = 1-np.exp(-lambda_)
        associated_this_pass = Xfiber
        if poiss: associated_this_pass *= p_associated
        observed_this_pass = associated_this_pass 
        cumulative += observed_this_pass
        Xremaining -= observed_this_pass
        nleft.append(Xremaining)
        if Xremaining < 1:
            return n, cumulative, Xremaining, nleft