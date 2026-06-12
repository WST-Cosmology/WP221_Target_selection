import numpy as np
import math, copy
import matplotlib.pyplot as plt
import _tracer_spectroscopic_efficiency as tracer_spectroscopic_efficiency
import sys
sys.path.append('../forecasts/')
import bias_model
import fisher_matrix_bao_SuEisenstein
import fisher_matrix_local_png
from itertools import product
import numpy as np
def linear_bias(redshift, mag, tracer = 'BG_faint'):

    if tracer == 'BG_faint' or tracer == 'BG_bright': return bias_model.bias_bg(redshift, mag)
    if tracer == 'ELG': return bias_model.bias_elg(redshift, mag)
    if tracer == 'LRG': return bias_model.bias_lrg(redshift, mag)
    if tracer == 'QSO': return bias_model.bias_qso(redshift, mag)
    if tracer == 'LBGu': return bias_model.bias_lbg(redshift, mag)
    if tracer == 'LBGg': return bias_model.bias_lbg(redshift, mag)
    if tracer == 'LBGr': return bias_model.bias_lbg(redshift, mag)
        
    return np.ones(len(redshift))

def Survey_design_science_metrics(config_survey_update, cosmo, redshift_eval_range = None, mag_max_eval_range = None, 
                                  one_mag_bin_approach = True,
                                  multi_mag_bin_approach = False):

    S_survey = config_survey_update['S_survey']
    N_fibres = config_survey_update['N_fibres']
    S_FoV = config_survey_update['S_FoV']
    t_exp = config_survey_update['exposure_time'] 
    observational_fraction = config_survey_update['observation_fraction']
    config_survey_update = copy.deepcopy(config_survey_update)

    forecasts = {}

    if one_mag_bin_approach:
    
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
                nspec_deg2 = nspec_deg2_frame[j]
                bz = linear_bias(zarray, mag_max, tracer = tracer)
                # BA0 PS constraints on Da and H
                list_zbin, list_sigma_Da,list_sigma_H, zeff, sigma_Da_eff, sigma_H_eff = fisher_matrix_bao_SuEisenstein.sigma_Da_H_single_tracer(z,
                                                                                                                np.interp(z, zarray, nz),
                                                                                                                np.interp(z, zarray, bz),
                                                                                                                S_survey,nspec_deg2,Deltaz=0.3,
                                                                                                                cosmo=cosmo)
    
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
                                                                                            S_survey,nspec_deg2,Deltaz=0.3,
                                                                                            p=1,mod='bbks',kmax=0.1,cosmo=cosmo)
                forecasts[tracer+'_list_zbin_fnl'].append(list_zbin)
                forecasts[tracer+'_list_sigma_fnl'].append(list_sigma_fnl)
                forecasts[tracer+'_sigma_fnl_eff'].append(sigma_fnl_eff)

    if multi_mag_bin_approach:
        
        for i, tracer in enumerate(config_survey_update['tracers']):
            
            if 'LBG' not in tracer: continue  #only works for LBG tracers

            mag_centers = config_survey_update[tracer + '_' + 'mag_centers']
            z_centers = config_survey_update[tracer + '_' + 'redshift_centers']
            mask_mag_max_eval_range = (mag_centers >= mag_max_eval_range[i][0])*(mag_centers <= mag_max_eval_range[i][1]) #mag_max to evaluate
            mag_centers_frame = mag_centers[mask_mag_max_eval_range]
            nz_distrib_frame = config_survey_update[tracer + '_' + 'spec_redshift_density'][mask_mag_max_eval_range]
            nspec_deg2_frame = config_survey_update[tracer + '_' + 'spec_density'][mask_mag_max_eval_range]
            
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
            
            print('Computing forecasts: Survey ', config_survey_update['survey_type'], ' --- tracer: ', tracer)
        
            mag_separation = 25
            mag_centers = config_survey_update[tracer + '_' + 'mag_centers']
            z_centers = config_survey_update[tracer + '_' + 'redshift_centers']
            mask_mag_max_eval_range = (mag_centers >= mag_max_eval_range[i][0])*(mag_centers <= mag_max_eval_range[i][1]) #mag_max to evaluate
            mag_centers_frame = mag_centers[mask_mag_max_eval_range]
            forecasts[tracer+'_mag_max_eval'] = mag_centers_frame
            nz_distrib_frame = config_survey_update[tracer + '_' + 'spec_redshift_density'][mask_mag_max_eval_range] #cumulative density
            nspec_deg2_frame = config_survey_update[tracer + '_' + 'spec_density'][mask_mag_max_eval_range]
            has_cross_mag_2bin = False
            index_max_mag_first_bin = 0
            mag_max_first_bin = 0
            mask_redshift_eval_range = (z_centers >= redshift_eval_range[i][0])*(z_centers <= redshift_eval_range[i][1])
            zarray = z_centers[mask_redshift_eval_range]
            z = np.linspace(zarray[0], zarray[-1], 100)
            for j, mag_max in enumerate(mag_centers_frame):
                if mag_max < mag_separation:
                    nz1 = nz_distrib_frame[j][mask_redshift_eval_range]
                    nspec1_deg2 = nspec_deg2_frame[j]
                    bz1 = linear_bias(zarray, mag_max, tracer = tracer)

                    list_zbin, list_sigma_fnl, zeff, sigma_fnl_eff = fisher_matrix_local_png.sigma_fnl_single_tracer(z,
                                                                                            np.interp(z, zarray, nz1),
                                                                                            np.interp(z, zarray, bz1),
                                                                                            S_survey, nspec1_deg2, Deltaz=0.3,
                                                                                            p=1,mod='bbks',kmax=0.1,cosmo=cosmo)
                    
                if mag_max >= mag_separation: 

                    if not has_cross_mag_2bin:
                        print('Now, forecasts are done considering 2 samples')
                        index_max_mag_first_bin = j - 1
                        mag_max_first_bin = mag_centers_frame[j - 1]
                        nz1 = nz_distrib_frame[j - 1][mask_redshift_eval_range]
                        nspec1_deg2 = nspec_deg2_frame[j - 1]
                        bz1 = linear_bias(zarray, mag_max_first_bin, tracer = tracer)
                        has_cross_mag_2bin = True

                    nz2 = nz_distrib_frame[j][mask_redshift_eval_range] - nz_distrib_frame[index_max_mag_first_bin][mask_redshift_eval_range] #low magn bin
                    nspec2_deg2 = nspec_deg2_frame[j] - nspec_deg2_frame[index_max_mag_first_bin]
                    bz2 = linear_bias(zarray, mag_max, tracer = tracer)
                    #print(f'1st sample; {nspec1_deg2:.0f} deg-2, 2nd sample {nspec2_deg2:.0f} deg-2)
                    

                    list_zbin, list_sigma_fnl, zeff, sigma_fnl_eff = fisher_matrix_local_png.sigma_fnl_two_tracers(z,
                                                                                            np.interp(z, zarray, nz1), np.interp(z, zarray, nz2),
                                                                                            np.interp(z, zarray, bz1), np.interp(z, zarray, bz2),
                                                                                            S_survey, nspec1_deg2, nspec2_deg2, Deltaz=0.3,
                                                                                            p=1,mod='bbks',kmax=0.1, cosmo=cosmo)

                forecasts[tracer+'_list_zbin_fnl'].append(list_zbin)
                forecasts[tracer+'_list_sigma_fnl'].append(list_sigma_fnl)
                forecasts[tracer+'_sigma_fnl_eff'].append(sigma_fnl_eff)

                    #declare sample 2
                        
            

    return forecasts

def build_total_survey_information_metrics(config_survey_update, forecasts_survey):

    tracers = config_survey_update['tracers']

    # shape of the parameter grid
    shape = tuple(
        len(forecasts_survey[f"{tracer}_mag_max_eval"])
        for tracer in tracers
    )

    total_time = np.zeros(shape)
    total_pointings_density = np.zeros(shape)
    total_efficiency = np.zeros(shape)
    Information_fnl = np.zeros(shape)
    Information_Da = np.zeros(shape)
    Information_H = np.zeros(shape)

    # iterate over all tracer index combinations
    for idx in product(*[range(n) for n in shape]):

        total_time[idx] = np.sum([config_survey_update[f"{tracer}_fibre_time"][m] for tracer, m in zip(tracers, idx)])
        total_pointings_density[idx] = np.sum([config_survey_update[f"{tracer}_target_density"][m] for tracer, m in zip(tracers, idx)])

        total_efficiency[idx] = (np.sum([config_survey_update[f"{tracer}_spec_density"][m]
                for tracer, m in zip(tracers, idx)])
                        /
                np.sum([config_survey_update[f"{tracer}_target_density"][m]
                for tracer, m in zip(tracers, idx)])
        )

        Information_fnl[idx] = np.sum([
            1.0 / forecasts_survey[f"{tracer}_sigma_fnl_eff"][m]**2
            for tracer, m in zip(tracers, idx)])

        Information_Da[idx] = np.sum([
            1.0 / forecasts_survey[f"{tracer}_sigma_Da_eff"][m]**2
            for tracer, m in zip(tracers, idx)])

        Information_H[idx] = np.sum([
            1.0 / forecasts_survey[f"{tracer}_sigma_H_eff"][m]**2
            for tracer, m in zip(tracers, idx)])

    Informations = {'mag_max_eval': [forecasts_survey[f"{tracer}_mag_max_eval"] for tracer in tracers],
                    'total_survey_time': total_time,
                    'total_pointings_density': total_pointings_density,
                    'total_survey_efficiency': total_efficiency,
                    'total_survey_fisher_information_fnl': Information_fnl,
                    'total_survey_fisher_information_Da': Information_Da,
                    'total_survey_fisher_information_H': Information_H}

    return Informations

def best_idx(FoM): return np.unravel_index( np.argmax(FoM), FoM.shape )

# def build_total_survey_information_metrics(config_survey_update, forecasts_survey):

#     total_time = np.zeros([len(forecasts_survey[tracer+'_mag_max_eval']) for tracer in config_survey_update['tracers']])
#     total_efficiency = np.zeros([len(forecasts_survey[tracer+'_mag_max_eval']) for tracer in config_survey_update['tracers']])
#     Information_fnl = np.zeros([len(forecasts_survey[tracer+'_mag_max_eval']) for tracer in config_survey_update['tracers']])
#     Information_Da = np.zeros([len(forecasts_survey[tracer+'_mag_max_eval']) for tracer in config_survey_update['tracers']])
#     Information_H = np.zeros([len(forecasts_survey[tracer+'_mag_max_eval']) for tracer in config_survey_update['tracers']])
#     for i in range(total_time.shape[0]):
#         for j in range(total_time.shape[1]):
#             for k in range(total_time.shape[2]):
#                 for l in range(total_time.shape[3]):
#                     total_time[i,j,k,l] = np.sum([config_survey_update[tracer + '_' + 'calendar_time'][m] for tracer, m in zip(config_survey_update['tracers'], [i,j,k,l])])
#                     total_efficiency[i,j,k,l] = np.sum([config_survey_update[tracer + '_' + 'spec_density'][m] for tracer, m in zip(config_survey_update['tracers'], [i,j,k,l])])/np.sum([config_dark_update[tracer + '_' + 'target_density'][m] for tracer, m in zip(config_survey_update['tracers'], [i,j,k,l])])
#                     Information_fnl[i,j,k,l] = np.sum([1/(forecasts_survey[tracer  + '_sigma_fnl_eff'][m]**2) for tracer, m in zip(config_survey_update['tracers'], [i,j,k,l])])
#                     Information_Da[i,j,k,l] = np.sum([1/(forecasts_survey[tracer  + '_sigma_Da_eff'][m]**2) for tracer, m in zip(config_survey_update['tracers'], [i,j,k,l])])
#                     Information_H[i,j,k,l] = np.sum([1/(forecasts_survey[tracer  + '_sigma_H_eff'][m]**2) for tracer, m in zip(config_survey_update['tracers'], [i,j,k,l])])

#     Informations = {}
#     Informations['total_survey_time'] = total_time
#     Informations['total_survey_efficiency'] = total_efficiency
#     Informations['total_survey_fisher_information_fnl'] = Information_fnl
#     Informations['total_survey_fisher_information_Da'] = Information_Da
#     Informations['total_survey_fisher_information_H'] = Information_H
    
#     return Informations
    
            


#for i,tracer in enumerate(config_dark_update['tracers']):
    
#    plt.plot(res[tracer+'_list_zbin_fnl'][2], res[tracer+'_list_sigma_fnl'][2], '-o',  color = config_dark_update['color'][i], label = tracer)
#    plt.ylabel(r'$\sigma(f_{\rm NL}^{\rm loc})$', fontsize=20)
#    plt.xlabel('redshift', fontsize=20)
#    plt.xlim(0, 6)
#    plt.ylim(0, 20)
#z = np.linspace(0, 6, 30)
#y = [np.interp(z, res[tracer+'_list_zbin_fnl'][2], 1/(np.array(res[tracer+'_list_sigma_fnl'][2])**2), left=0.000000000000001, right=0.0000000000001) #for tracer in config_dark_update['tracers']]
#plt.plot(z, np.sum(y, axis=0)**(-0.5), 'orange', ls='--', label = r'combined per $\Delta z$')
#all = np.sum(np.sum(y, axis=0))
#plt.axhline(all, label = 'combined all redshifts')
#plt.legend()