import numpy as np
import math, copy
import matplotlib.pyplot as plt
import _tracer_spectroscopic_efficiency as tracer_spectroscopic_efficiency

def passes_needed(Xtarget_init, Xfiber, X_I_want, max_passes=2000, poiss=True):
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

def plot_diagnostics(config_survey_update, max_mag = None):

    plt.subplot(141)
    for j, tracer in enumerate(config_survey_update['tracers']):
        
        plt.plot( config_survey_update[tracer + '_' + 'mag_centers'], 
                  config_survey_update[tracer + '_' + 'target_density'], f'-', 
                  color = config_survey_update['color'][j], label = tracer)
        plt.plot( config_survey_update[tracer + '_' + 'mag_centers'], 
                  config_survey_update[tracer + '_' + 'spec_density'], f'--', 
                  color = config_survey_update['color'][j])

    plt.ylabel(r'Density (deg$^{-2}$)', fontsize=12)
    plt.xlabel('maximum magnitude', fontsize=12)
    plt.yscale('log')
    plt.legend()
    
    plt.subplot(142)
    for j, tracer in enumerate(config_survey_update['tracers']):
        plt.plot( config_survey_update[tracer + '_' + 'mag_centers'], 
                 config_survey_update[tracer + '_' + 'spec_density']/config_survey_update[tracer + '_' + 'target_density'], f'',
                  color = config_survey_update['color'][j], label = tracer)
    plt.ylabel(r'Efficiency', fontsize=12)
    plt.xlabel('maximum magnitude', fontsize=12)
    plt.legend()

    plt.subplot(143)
    for j, tracer in enumerate(config_survey_update['tracers']):
        plt.plot( config_survey_update[tracer + '_' + 'mag_centers'], 
                  config_survey_update[tracer + '_' + 'calendar_time'], f'-',
                  color = config_survey_update['color'][j], label = tracer)
    plt.yscale('log')
    plt.xlabel('maximum magnitude', fontsize=12)
    plt.ylabel('Survey time per tracer', fontsize=12)
    plt.axhline(5, color='k', ls='--')

    plt.subplot(144)
    for j, tracer in enumerate(config_survey_update['tracers']):
        index = np.argmin(abs(config_survey_update[tracer + '_' + 'mag_centers'] - max_mag[j]))
        n = config_survey_update[tracer + '_' + 'spec_density'][index]
        z, x = config_survey_update[tracer + '_' + 'redshift_centers'], config_survey_update[tracer + '_' + 'spec_redshift_density'][index]
        plt.step(z ,x, f'-',  color = config_survey_update['color'][j], where='mid', label = tracer + f'_m_max={max_mag[j]}, n={n:.0f} deg-2')

        z, x = config_survey_update[tracer + '_' + 'redshift_centers'], config_survey_update[tracer + '_' + 'target_redshift_density'][index]
        plt.step(z ,x, f'--',  color = config_survey_update['color'][j], where='mid',)

    plt.ylabel(r'Density', fontsize=12)
    plt.xlabel('redshift', fontsize=12)
    plt.legend()

def Survey_design_telescope_metrics(config_survey, max_mag = None):

    S_survey = config_survey['S_survey']
    N_fibres = config_survey['N_fibres']
    S_FoV = config_survey['S_FoV']
    t_exp = config_survey['exposure_time'] 
    observational_fraction = config_survey['observation_fraction']
    config_survey_update = copy.deepcopy(config_survey)
    for i, tracer in enumerate(config_survey['tracers']):
        N_zm = np.load(config_survey['tracer_N_zm_file'][i])
        mag_centers = N_zm['mag_center']
        z_centers = N_zm['z_center']
        n_target_count = N_zm['object_count'] / N_zm['surface_deg2']
        Efficiency = np.zeros([len(z_centers), len(mag_centers)])
        n_pass = np.zeros([len(z_centers), len(mag_centers)])
        for i, z in enumerate(z_centers):
            Efficiency[i,:] = tracer_spectroscopic_efficiency.E_wst(z, mag_centers, tracer = tracer)
            n_pass[i,:] = tracer_spectroscopic_efficiency.n_pass_wst(z, mag_centers, tracer = tracer)

        n_pointings = []
        n_target = []
        n_spec = []
        n_specz_redshift = []
        n_target_redshift = []
        
        for m in mag_centers:
            n_target.append(np.sum(np.sum(n_target_count[:, mag_centers <= m], axis=1), axis=0))
            n_spec.append(np.sum(np.sum((n_target_count * Efficiency)[:, mag_centers <= m], axis=1), axis=0))
            n_pointings.append(np.sum(np.sum((n_target_count * n_pass)[:, mag_centers <= m], axis=1), axis=0))
            n_specz_redshift.append(np.sum((n_target_count * Efficiency)[:, mag_centers <= m], axis=1))
            n_target_redshift.append(np.sum((n_target_count)[:, mag_centers <= m], axis=1))
        config_survey_update[tracer + '_' + 'target_density'] = np.array(n_target)
        config_survey_update[tracer + '_' + 'spec_density'] = np.array(n_spec)
        config_survey_update[tracer + '_' + 'spec_redshift_density'] = np.array(n_specz_redshift)
        config_survey_update[tracer + '_' + 'target_redshift_density'] = np.array(n_target_redshift)
        config_survey_update[tracer + '_' + 'target_pointings'] = np.array(n_pointings)
        config_survey_update[tracer + '_' + 'fibre_time'] = np.array(n_pointings) * t_exp * (S_survey / N_fibres) / (365.25 * 24 * 3600)
        config_survey_update[tracer + '_' + 'calendar_time'] = config_survey_update[tracer + '_' + 'fibre_time'] / observational_fraction
        config_survey_update[tracer + '_' + 'mag_centers'] = np.array(mag_centers)
        config_survey_update[tracer + '_' + 'redshift_centers'] = np.array(z_centers)

    if max_mag != None:
        print(' We compute the survey completeness C(n) after n >= 1 passes')
        print(' After a maximum number of passes, the survey is said complete ')
        time = 0
        n_pointings = 0
        t_full_sky_one_exp = t_exp * (S_survey / S_FoV)
        for j, tracer in enumerate(config_survey_update['tracers']):
            index = np.argmin(abs(config_survey_update[tracer + '_' + 'mag_centers'] - max_mag[j]))
            n_pointings += config_survey_update[tracer + '_' + 'target_pointings'][index]
    
        n_fibres= N_fibres/S_FoV #density of fibres
        n, cumulative, Xremaining, nleft = passes_needed(n_pointings, n_fibres, n_pointings, max_passes=2000, poiss=True)

        number_passes, completeness = [], []
        for j in range(1,len(nleft)):
            comp = 1 - nleft[j]/n_pointings
            completeness.append(comp)
            number_passes.append(j)
            if comp > 0.999: continue

        config_survey_update['total_number_of_passes'] = np.array(number_passes)
        config_survey_update['total_survey_completeness'] = np.array(completeness)
        config_survey_update['total_survey_fibre_time'] = np.array(t_full_sky_one_exp) * np.array(number_passes) / (365.25 * 24 * 3600)
        config_survey_update['total_survey_calendar_time'] = config_survey_update['total_survey_fibre_time'] / observational_fraction
        
    return config_survey_update




































































# def t_tel(m,t1exp=15,eta=1.2,Stel=Swst,SNR=10,tref=2*60, mref=24.5, SNRref=10):

#     def I(m):#flux as function of magnitude. To be checked.
#         I0=1# anyway we will consider ratio of flux
#         return I0*10**(-0.4*m)
    
#     Sdesi=math.pi*4**2 #desi is a 4 meter
    
#     x=tref/t1exp *Sdesi/Stel*(I(mref)*SNR/(eta*I(m)*SNRref))**2
    
#     return math.ceil(x)*t1exp #ceil is partie entiere sup
    
# def Tsur(m_array,Phi_m,texp_m,m_max):
#     '''
#     Time allocated for the observation of a tracer up to a magnitude m_max
    
#     m_array is the array of magnitude
#     Phi_m the target density for a certain magnitude (per bin of magnitude) in deg-2
#     texp_m the array of exposure time for a given magnitude
#     m_max the maximal magnitude to be considered
#     '''
#     prefact=fsky*41253*eta_til/Nfib
#     sel=((m_array<m_max))
#     integ=np.sum(Phi_m[sel]*texp_m[sel])
#     return prefact*integ