import numpy as np
import math


def t_tel(m,t1exp=15,eta=1.2,Stel=Swst,SNR=10,tref=2*60, mref=24.5, SNRref=10):

    def I(m):#flux as function of magnitude. To be checked.
        I0=1# anyway we will consider ratio of flux
        return I0*10**(-0.4*m)
    
    Sdesi=math.pi*4**2 #desi is a 4 meter
    
    x=tref/t1exp *Sdesi/Stel*(I(mref)*SNR/(eta*I(m)*SNRref))**2
    
    return math.ceil(x)*t1exp #ceil is partie entiere sup
    
def Tsur(m_array,Phi_m,texp_m,m_max):
    '''
    Time allocated for the observation of a tracer up to a magnitude m_max
    
    m_array is the array of magnitude
    Phi_m the target density for a certain magnitude (per bin of magnitude) in deg-2
    texp_m the array of exposure time for a given magnitude
    m_max the maximal magnitude to be considered
    '''
    prefact=fsky*41253*eta_til/Nfib
    sel=((m_array<m_max))
    integ=np.sum(Phi_m[sel]*texp_m[sel])
    return prefact*integ