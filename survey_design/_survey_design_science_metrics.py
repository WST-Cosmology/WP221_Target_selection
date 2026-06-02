import numpy as np
import math, copy
import matplotlib.pyplot as plt
import _tracer_spectroscopic_efficiency as tracer_spectroscopic_efficiency
import sys
sys.path.append('../forecasts/')
import bias_model
import fisher_matrix_bao_SuEisenstein
import fisher_matrix_local_png

def linear_bias(redshift, mag, tracer = 'BG_faint'):

    if tracer == 'QSO': return bias_model.bz_model_desiqso(redshift, mag)
    if tracer == 'LBGu': return bias_model.bias_lbg(redshift, mag)
    if tracer == 'LBGg': return bias_model.bias_lbg(redshift, mag)
    if tracer == 'LBGr': return bias_model.bias_lbg(redshift, mag)
        
    return np.ones(len(redshift))

def Survey_design_science_metrics(config_survey_update, cosmo, redshift_eval_range = None, mag_max_eval_range = None):

    S_survey = config_survey_update['S_survey']
    N_fibres = config_survey_update['N_fibres']
    S_FoV = config_survey_update['S_FoV']
    t_exp = config_survey_update['exposure_time'] 
    observational_fraction = config_survey_update['observation_fraction']
    config_survey_update = copy.deepcopy(config_survey_update)

    forecasts = {}
    
    for i, tracer in enumerate(config_survey_update['tracers']):

        mag_centers = config_survey_update[tracer + '_' + 'mag_centers']
        z_centers = config_survey_update[tracer + '_' + 'redshift_centers']
        mask_mag_max_eval_range = (mag_centers >= mag_max_eval_range[i][0])*(mag_centers <= mag_max_eval_range[i][1]) #mag_max to evaluate
        mag_centers_frame = mag_centers[mask_mag_max_eval_range]
        nz_distrib_frame = config_survey_update[tracer + '_' + 'spec_redshift_density'][mask_mag_max_eval_range]
        nspec_deg2_frame = config_survey_update[tracer + '_' + 'spec_density'][mask_mag_max_eval_range]

        print('Computing forecasts: Survey ', config_survey_update['survey_type'], ' --- tracer: ', tracer)

        forecasts[tracer+'_mag_max_eval'] = mag_centers_frame
        forecasts[tracer+'_list_zbin_Da'] = []
        forecasts[tracer+'_list_sigma_Da'] = []
        forecasts[tracer+'_sigma_Da_eff'] = []

        forecasts[tracer+'_list_zbin_H'] = []
        forecasts[tracer+'_list_sigma_H'] = []
        forecasts[tracer+'_sigma_H_eff'] = []

        forecasts[tracer+'_list_zbin_fnl'] = []
        forecasts[tracer+'_list_sigma_fnl'] = []
        forecasts[tracer+'_sigma_fnl_eff'] = []

        for j, mag_max in enumerate(mag_centers_frame):
            mask_redshift_eval_range = (z_centers >= redshift_eval_range[i][0])*(z_centers <= redshift_eval_range[i][1])
            zarray = z_centers[mask_redshift_eval_range]
            z = np.linspace(zarray[0], zarray[-1], 100)
            nz = nz_distrib_frame[j][mask_redshift_eval_range]
            #plt.plot(zarray, nz)
            nspec_deg2 = nspec_deg2_frame[j]
            bz = linear_bias(zarray, mag_max, tracer = tracer)
            # BA0 PS constraints on Da and H
            list_zbin, list_sigma_Da,list_sigma_H, zeff, sigma_Da_eff, sigma_H_eff = fisher_matrix_bao_SuEisenstein.sigma_Da_H_single_tracer(z,
                                                                                                            np.interp(z, zarray, nz),
                                                                                                            np.interp(z, zarray, bz),
                                                                                                            S_survey,nspec_deg2,
                                                                                                            Deltaz=0.2,cosmo=cosmo)

            forecasts[tracer+'_list_zbin_Da'].append(list_zbin)
            forecasts[tracer+'_list_sigma_Da'].append(list_sigma_Da)
            forecasts[tracer+'_sigma_Da_eff'].append(sigma_Da_eff)

            forecasts[tracer+'_list_zbin_H'].append(list_zbin)
            forecasts[tracer+'_list_sigma_H'].append(list_sigma_H)
            forecasts[tracer+'_sigma_H_eff'].append(sigma_H_eff)
            
            # PS constraints on fnl
            list_zbin, list_sigma_fnl, zeff, sigma_fnl_eff = fisher_matrix_local_png.sigma_fnl_single_tracer(z,
                                                                                        np.interp(z, zarray, nz),
                                                                                        np.interp(z, zarray, bz),
                                                                                        S_survey,nspec_deg2,
                                                                                        Deltaz=0.2,p=1,mod='bbks',kmax=0.1,cosmo=cosmo)
            forecasts[tracer+'_list_zbin_fnl'].append(list_zbin)
            forecasts[tracer+'_list_sigma_fnl'].append(list_sigma_fnl)
            forecasts[tracer+'_sigma_fnl_eff'].append(sigma_fnl_eff)

    return forecasts

        
        