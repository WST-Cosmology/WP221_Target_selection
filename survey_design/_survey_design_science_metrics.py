import numpy as np
import math, copy
import matplotlib.pyplot as plt
import _tracer_spectroscopic_efficiency as tracer_spectroscopic_efficiency
import _survey_design_telescope_metrics as survey_design_telescope_metrics
import sys
import itertools
import numpy as np
sys.path.append('../forecasts/')
import bias_model
import power_spectrum_information
import fisher_matrix_local_png
import fisher_matrix_bao_SuEisenstein
import fisher_matrix_rsd
import fisher_matrix_neutrino_mass

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

def run_forecast_one_tracer(z, zarray, nz, bz, S_survey, nspec_deg2, cosmo):

    list_zbin_fnl, list_sigma_fnl, zeff, sigma_fnl_eff = fisher_matrix_local_png.sigma_fnl_single_tracer(z,
                                                                                np.interp(z, zarray, nz),
                                                                                np.interp(z, zarray, bz),
                                                                                S_survey,nspec_deg2,Deltaz=0.3,
                                                                                p=1,mod='bbks',kmax=0.1,cosmo=cosmo)
    
    list_zbin_bao, list_sigma_Da,list_sigma_H, zeff, sigma_Da_eff, sigma_H_eff = fisher_matrix_bao_SuEisenstein.sigma_Da_H_single_tracer(z,
                                                                                                                np.interp(z, zarray, nz),
                                                                                                                np.interp(z, zarray, bz),
                                                                                                                S_survey,nspec_deg2,Deltaz=0.3,
                                                                                                                cosmo=cosmo)
    
    list_zbin_mnu, list_sigma_b, list_sigma_mnu, zeff, sigma_b_eff, sigma_mnu_eff=fisher_matrix_neutrino_mass.sigma_mnu_single_tracer(z,
                                                                                                                np.interp(z, zarray, nz),
                                                                                                                np.interp(z, zarray, bz),
                                                                                                                S_survey,nspec_deg2,Deltaz=0.3,
                                                                                                                kmax=0.1,cosmo=cosmo)

    list_zbin_rsd, list_sigma_bs8, list_sigma_fs8, zeff, sigma_bs8_eff, sigma_fs8_eff=fisher_matrix_rsd.sigma_rsd_single_tracer(z,
                                                                                                                np.interp(z, zarray, nz),
                                                                                                                np.interp(z, zarray, bz),
                                                                                                                S_survey,nspec_deg2,Deltaz=0.3,
                                                                                                                kmax=0.1,cosmo=cosmo)

    results = {}
    results['fnl'] = [list_zbin_fnl, list_sigma_fnl, zeff, sigma_fnl_eff]
    results['bao'] = [ list_zbin_bao, list_sigma_Da,list_sigma_H, zeff, sigma_Da_eff, sigma_H_eff]
    results['neutrinos'] = [list_zbin_mnu, list_sigma_b, list_sigma_mnu, zeff, sigma_b_eff, sigma_mnu_eff]
    results['rsd'] = [list_zbin_rsd, list_sigma_bs8, list_sigma_fs8, zeff, sigma_bs8_eff, sigma_fs8_eff]

    return results

def run_forecast_two_tracers(z, zarray, nz1, nz2, bz1, bz2, S_survey, nspec1_deg2, nspec2_deg2, cosmo):

    list_zbin_fnl, list_sigma_fnl, zeff_fnl, sigma_fnl_eff= fisher_matrix_local_png.sigma_fnl_two_tracers(z,
                                                                        np.interp(z, zarray, nz1), np.interp(z, zarray, nz2),
                                                                        np.interp(z, zarray, bz1), np.interp(z, zarray, bz2),
                                                                        S_survey, nspec1_deg2, nspec2_deg2, Deltaz=0.3,
                                                                        p=1,mod='bbks',kmax=0.1, cosmo=cosmo)
    
    list_zbin_mnu, list_sigma_ba, list_sigma_bb, list_sigma_mnu, zeff_mnu, sigma_ba_eff, sigma_bb_eff, sigma_mnu_eff = fisher_matrix_neutrino_mass.sigma_mnu_two_tracers(z,
                                                                        np.interp(z, zarray, nz1), np.interp(z, zarray, nz2),
                                                                        np.interp(z, zarray, bz1), np.interp(z, zarray, bz2),
                                                                        S_survey, nspec1_deg2, nspec2_deg2, Deltaz=0.3,
                                                                        kmax=0.1,
                                                                        Sigma_mnu_fid=0.06, dSigma=0.06,
                                                                        cosmo=cosmo, Nk=200,
                                                                        return_F=False)

    list_zbin_rsd, list_sigma_bAs8, list_sigma_bBs8, list_sigma_fs8, zeff_rsd, sigma_bAs8_eff, sigma_bBs8_eff, sigma_fs8_eff = fisher_matrix_rsd.sigma_rsd_two_tracers(
                                                                        z, np.interp(z, zarray, nz1), np.interp(z, zarray, nz2),
                                                                        np.interp(z, zarray, bz1), np.interp(z, zarray, bz2),
                                                                        S_survey, nspec1_deg2, nspec2_deg2, 
        Deltaz=0.3, 
                                                                        kmax=0.1,
                                                                        cosmo=cosmo, Nk=100, Nmu=50,
                                                                        return_F=False)

    results = {}
    results['fnl'] = [list_zbin_fnl, list_sigma_fnl, zeff_fnl, sigma_fnl_eff]
    results['neutrinos'] = [list_zbin_mnu, list_sigma_ba, list_sigma_bb, list_sigma_mnu, zeff_mnu, sigma_ba_eff, sigma_bb_eff, sigma_mnu_eff]
    results['rsd'] = [list_zbin_rsd, list_sigma_bAs8, list_sigma_bBs8, list_sigma_fs8, zeff_rsd, sigma_bAs8_eff, sigma_bBs8_eff, sigma_fs8_eff]

    return results

def compute_FnP_one_tracer(z, zarray, nz, bz, S_survey, nspec_deg2, cosmo):

    results = {}
    for k in [0.001, 0.1, 1]:
        list_zbin, list_n, list_b, list_Pm, list_nP, list_Vsur = power_spectrum_information.compute_nbP(z, np.interp(z, zarray, nz), np.interp(z, zarray, bz),
                                                                          S_survey, nspec_deg2, k=k, Deltaz=0.3, cosmo=cosmo)
        results[f'list_nP_tracer_k{k}'] = np.array(list_nP)
        list_Vtracer = np.array(list_Vsur)
        results[f'nP_tracer_eff_k{k}'] = np.average(results[f'list_nP_tracer_k{k}'], weights=list_Vtracer)
        results[f'F_tracer_k{k}'] = np.sum(list_Vtracer * (results[f'list_nP_tracer_k{k}']/(results[f'list_nP_tracer_k{k}']+1))**2) / 1e10

    return results

def compute_FnP_two_tracers(z, zarray, nz1, nz2, bz1, bz2, S_survey, nspec1_deg2, nspec2_deg2, cosmo):

    results = {}
    for k in [0.001, 0.1, 1]:
        list_zbin, list_nA, list_nB, list_bA, list_bB, list_Pm, list_nPA, list_nPB, list_Vsur = power_spectrum_information.compute_nbP_two_tracers(
                                                                        z, np.interp(z, zarray, nz1), np.interp(z, zarray, nz2),
                                                                        np.interp(z, zarray, bz1), np.interp(z, zarray, bz2),
                                                                        S_survey, nspec1_deg2, nspec2_deg2, 
                                                                        k=k, Deltaz=0.3,cosmo=cosmo)
        results[f'list_nP_tracer_k{k}'] = np.array(list_nPA) + np.array(list_nPB)
        list_Vtracer = np.array(list_Vsur)

        results[f'nP_tracer_eff_k{k}'] = np.average(results[f'list_nP_tracer_k{k}'], weights=list_Vtracer)
        results[f'F_tracer_k{k}'] = np.sum(list_Vtracer * (results[f'list_nP_tracer_k{k}']/(results[f'list_nP_tracer_k{k}']+1))**2) / 1e10
    
    return results


def Survey_design_science_metrics(config_survey_update, cosmo, redshift_eval_range = None, mag_max_eval_range = None, multi_mag_bin_approach = False):

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
        mask_redshift_eval_range = (z_centers >= redshift_eval_range[i][0])*(z_centers <= redshift_eval_range[i][1])

        print('Computing forecasts: Survey ', config_survey_update['survey_type'], ' --- tracer: ', tracer)

        forecasts[tracer+'_mag_max_eval'] = mag_centers_frame

        forecasts[tracer+'_list_zbin_fnl'] = []
        forecasts[tracer+'_zeff_fnl'] = []
        forecasts[tracer+'_list_sigma_fnl'] = []
        forecasts[tracer+'_sigma_fnl_eff'] = []
        
        forecasts[tracer+'_list_zbin_Da'] = []
        forecasts[tracer+'_zeff_Da'] = []
        forecasts[tracer+'_list_sigma_Da'] = []
        forecasts[tracer+'_sigma_Da_eff'] = []

        forecasts[tracer+'_list_zbin_H'] = []
        forecasts[tracer+'_zeff_H'] = []
        forecasts[tracer+'_list_sigma_H'] = []
        forecasts[tracer+'_sigma_H_eff'] = []

        forecasts[tracer+'_list_zbin_Mnu'] = []
        forecasts[tracer+'_zeff_Mnu'] = []
        forecasts[tracer+'_list_sigma_Mnu'] = []
        forecasts[tracer+'_sigma_Mnu_eff'] = []

        forecasts[tracer+'_list_zbin_rsd'] = []
        forecasts[tracer+'_zeff_rsd'] = []
        forecasts[tracer+'_list_sigma_rsd'] = []
        forecasts[tracer+'_sigma_rsd_eff'] = []

        zarray = z_centers[mask_redshift_eval_range]
        z = np.linspace(zarray[0], zarray[-1], 100)
        mag_separation = 25
        has_cross_mag_2bin = False
        index_max_mag_first_bin = 0
        mag_max_first_bin = 0

        for j, mag_max in enumerate(mag_centers_frame):

            if not multi_mag_bin_approach[i]:

                nz = nz_distrib_frame[j][mask_redshift_eval_range]
                nspec_deg2 = nspec_deg2_frame[j]
                bz = linear_bias(zarray, mag_max, tracer = tracer)
                # BA0 PS constraints on Da and H
                results = run_forecast_one_tracer(z, zarray, nz, bz, S_survey, nspec_deg2, cosmo)

                list_zbin_fnl, list_sigma_fnl, zeff_fnl, sigma_fnl_eff = results['fnl']
                list_zbin_bao, list_sigma_Da,list_sigma_H, zeff_bao, sigma_Da_eff, sigma_H_eff = results['bao']
                list_zbin_mnu, list_sigma_b, list_sigma_mnu, zeff_mnu, sigma_b_eff, sigma_mnu_eff = results['neutrinos']
                list_zbin_rsd, list_sigma_bs8, list_sigma_fs8, zeff_rsd, sigma_bs8_eff, sigma_fs8_eff = results['rsd']

                forecasts[tracer+'_list_zbin_fnl'].append(list_zbin_fnl)
                forecasts[tracer+'_zeff_fnl'].append(zeff_fnl)
                forecasts[tracer+'_list_sigma_fnl'].append(list_sigma_fnl)
                forecasts[tracer+'_sigma_fnl_eff'].append(sigma_fnl_eff)
    
                forecasts[tracer+'_list_zbin_Da'].append(list_zbin_bao)
                forecasts[tracer+'_zeff_Da'].append(zeff_bao)
                forecasts[tracer+'_list_sigma_Da'].append(list_sigma_Da)
                forecasts[tracer+'_sigma_Da_eff'].append(sigma_Da_eff)
    
                forecasts[tracer+'_list_zbin_H'].append(list_zbin_bao)
                forecasts[tracer+'_zeff_H'].append(zeff_bao)
                forecasts[tracer+'_list_sigma_H'].append(list_sigma_H)
                forecasts[tracer+'_sigma_H_eff'].append(sigma_H_eff)

                forecasts[tracer+'_list_zbin_Mnu'].append(list_zbin_mnu)
                forecasts[tracer+'_zeff_Mnu'].append(zeff_mnu)
                forecasts[tracer+'_list_sigma_Mnu'].append(list_sigma_mnu)
                forecasts[tracer+'_sigma_Mnu_eff'].append(sigma_mnu_eff)

                forecasts[tracer+'_list_zbin_rsd'].append(list_zbin_rsd)
                forecasts[tracer+'_zeff_rsd'].append(zeff_rsd)
                forecasts[tracer+'_list_sigma_rsd'].append(list_sigma_fs8)
                forecasts[tracer+'_sigma_rsd_eff'].append(sigma_fs8_eff)

            elif multi_mag_bin_approach[i]:

                if mag_max < mag_separation:
                    nz1 = nz_distrib_frame[j][mask_redshift_eval_range]
                    nspec1_deg2 = nspec_deg2_frame[j]
                    bz1 = linear_bias(zarray, mag_max, tracer = tracer)

                    results1 = run_forecast_one_tracer(z, zarray, nz1, bz1, S_survey, nspec1_deg2, cosmo)

                    list_zbin_fnl_1, list_sigma_fnl_1, zeff_fnl_1, sigma_fnl_eff_1 = results1['fnl']
                    list_zbin_bao_1, list_sigma_Da_1,list_sigma_H_1, zeff_bao_1, sigma_Da_eff_1, sigma_H_eff_1 = results1['bao']
                    list_zbin_mnu_1, list_sigma_b_1, list_sigma_mnu_1, zeff_mnu_1, sigma_b_eff_1, sigma_mnu_eff_1 = results1['neutrinos']
                    list_zbin_rsd_1, list_sigma_bs8_1, list_sigma_fs8_1, zeff_rsd_1, sigma_bs8_eff_1, sigma_fs8_eff_1 = results1['rsd']

                    sigma_fnl_eff_joint = sigma_fnl_eff_1
                    sigma_Da_eff_joint, sigma_H_eff_joint = sigma_Da_eff_1, sigma_H_eff_1
                    sigma_mnu_eff_joint = sigma_mnu_eff_1
                    sigma_fs8_eff_joint = sigma_fs8_eff_1
                    
                else: 

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

                    #alone constraints
                    results2 = run_forecast_one_tracer(z, zarray, nz2, bz2, S_survey, nspec2_deg2, cosmo)

                    list_zbin_fnl_2, list_sigma_fnl_2, zeff_fnl_2, sigma_fnl_eff_2 = results2['fnl']
                    list_zbin_bao_2, list_sigma_Da_2,list_sigma_H_2, zeff_bao_2, sigma_Da_eff_2, sigma_H_eff_2 = results2['bao']
                    list_zbin_mnu_2, list_sigma_b_2, list_sigma_mnu_2, zeff_mnu_2, sigma_b_eff_2, sigma_mnu_eff_2 = results2['neutrinos']
                    list_zbin_rsd_2, list_sigma_bs8_2, list_sigma_fs8_2, zeff_rsd_2, sigma_bs8_eff_2, sigma_fs8_eff_2 = results2['rsd']

                    #joint constraints
                    results12 = run_forecast_two_tracers(z, zarray, nz1, nz2, bz1, bz2, S_survey, nspec1_deg2, nspec2_deg2, cosmo)
                    list_zbin_fnl, list_sigma_fnl, zeff, sigma_fnl_eff = results12['fnl']
                    list_zbin_mnu, list_sigma_ba, list_sigma_bb, list_sigma_mnu, zeff, sigma_ba_eff, sigma_bb_eff, sigma_mnu_eff = results12['neutrinos'] 
                    list_zbin_rsd, list_sigma_bAs8, list_sigma_bBs8, list_sigma_fs8, zeff_rsd, sigma_bAs8_eff, sigma_bBs8_eff, sigma_fs8_eff = results12['rsd']

                    sigma_fnl_eff_joint = sigma_fnl_eff
                    sigma_Da_eff_joint, sigma_H_eff_joint = (1/(sigma_Da_eff_1**2)+1/(sigma_Da_eff_2**2))**(-0.5), (1/(sigma_H_eff_1**2)+1/(sigma_H_eff_2**2))**(-0.5)
                    sigma_mnu_eff_joint = sigma_mnu_eff
                    sigma_fs8_eff_joint = sigma_fs8_eff

                
                forecasts[tracer+'_list_zbin_fnl'].append(list_zbin_fnl_1)
                forecasts[tracer+'_zeff_fnl'].append(zeff_fnl_1)
                forecasts[tracer+'_list_sigma_fnl'].append(list_sigma_fnl_1)
                forecasts[tracer+'_sigma_fnl_eff'].append(sigma_fnl_eff_joint)

                forecasts[tracer+'_list_zbin_Da'].append(list_zbin_bao_1)
                forecasts[tracer+'_zeff_Da'].append(zeff_bao_1)
                forecasts[tracer+'_list_sigma_Da'].append(None)
                forecasts[tracer+'_sigma_Da_eff'].append(sigma_Da_eff_joint)
    
                forecasts[tracer+'_list_zbin_H'].append(list_zbin_bao_1)
                forecasts[tracer+'_zeff_H'].append(zeff_bao_1)
                forecasts[tracer+'_list_sigma_H'].append(None)
                forecasts[tracer+'_sigma_H_eff'].append(sigma_H_eff_joint)

                forecasts[tracer+'_list_zbin_Mnu'].append(list_zbin_mnu_1)
                forecasts[tracer+'_zeff_Mnu'].append(zeff_mnu_1)
                forecasts[tracer+'_list_sigma_Mnu'].append(None)
                forecasts[tracer+'_sigma_Mnu_eff'].append(sigma_mnu_eff_joint)
    
                forecasts[tracer+'_list_zbin_rsd'].append(list_zbin_rsd_1)
                forecasts[tracer+'_zeff_rsd'].append(zeff_rsd_1)
                forecasts[tracer+'_list_sigma_rsd'].append(None)
                forecasts[tracer+'_sigma_rsd_eff'].append(sigma_fs8_eff_joint)

    return forecasts

def Survey_design_nP_metrics(config_survey_update, cosmo, redshift_eval_range = None, mag_max_eval_range = None, multi_mag_bin_approach = False):

    S_survey = config_survey_update['S_survey']
    N_fibres = config_survey_update['N_fibres']
    S_FoV = config_survey_update['S_FoV']
    t_exp = config_survey_update['exposure_time'] 
    observational_fraction = config_survey_update['observation_fraction']
    config_survey_update = copy.deepcopy(config_survey_update)

    nP = {}
    
    for i, tracer in enumerate(config_survey_update['tracers']):

        mag_centers = config_survey_update[tracer + '_' + 'mag_centers']
        z_centers = config_survey_update[tracer + '_' + 'redshift_centers']
        mask_mag_max_eval_range = (mag_centers >= mag_max_eval_range[i][0])*(mag_centers <= mag_max_eval_range[i][1]) #mag_max to evaluate
        mag_centers_frame = mag_centers[mask_mag_max_eval_range]
        nz_distrib_frame = config_survey_update[tracer + '_' + 'spec_redshift_density'][mask_mag_max_eval_range]
        nspec_deg2_frame = config_survey_update[tracer + '_' + 'spec_density'][mask_mag_max_eval_range]
        mask_redshift_eval_range = (z_centers >= redshift_eval_range[i][0])*(z_centers <= redshift_eval_range[i][1])

        print('Computing forecasts: Survey ', config_survey_update['survey_type'], ' --- tracer: ', tracer)

        nP[tracer+'_mag_max_eval'] = mag_centers_frame
        nP[tracer+'_nP_eff_k0.001'] = []
        nP[tracer+'_nP_eff_k0.1'] = []
        nP[tracer+'_nP_eff_k1'] = []
        nP[tracer+'_F_k0.001'] = []
        nP[tracer+'_F_k0.1'] = []
        nP[tracer+'_F_k1'] = []

        zarray = z_centers[mask_redshift_eval_range]
        z = np.linspace(zarray[0], zarray[-1], 100)
        mag_separation = 25
        has_cross_mag_2bin = False
        index_max_mag_first_bin = 0
        mag_max_first_bin = 0


        for j, mag_max in enumerate(mag_centers_frame):

            if not multi_mag_bin_approach[i]:

                nz = nz_distrib_frame[j][mask_redshift_eval_range]
                nspec_deg2 = nspec_deg2_frame[j]
                bz = linear_bias(zarray, mag_max, tracer = tracer)

                res_nP = compute_FnP_one_tracer(z, zarray, nz, bz, S_survey, nspec_deg2, cosmo)

            elif multi_mag_bin_approach[i]:

                if mag_max < mag_separation:
                    nz1 = nz_distrib_frame[j][mask_redshift_eval_range]
                    nspec1_deg2 = nspec_deg2_frame[j]
                    bz1 = linear_bias(zarray, mag_max, tracer = tracer)

                    res_nP = compute_FnP_one_tracer(z, zarray, nz1, bz1, S_survey, nspec1_deg2, cosmo)
                    
                else: 

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

                    res_nP = compute_FnP_two_tracers(z, zarray, nz1, nz2, bz1, bz2, S_survey, nspec1_deg2, nspec2_deg2, cosmo)

            nP[tracer+'_nP_eff_k0.001'].append(res_nP['nP_tracer_eff_k0.001'])
            nP[tracer+'_nP_eff_k0.1'].append(res_nP['nP_tracer_eff_k0.1'])
            nP[tracer+'_nP_eff_k1'].append(res_nP['nP_tracer_eff_k1'])
            nP[tracer+'_F_k0.001'].append(res_nP['F_tracer_k0.001'])
            nP[tracer+'_F_k0.1'].append(res_nP['F_tracer_k0.1'])
            nP[tracer+'_F_k1'].append(res_nP['F_tracer_k1'])

    return nP

def best_idx(FoM): return np.unravel_index( np.argmax(FoM), FoM.shape )

def build_total_survey_information_metrics(config_survey_update, forecasts_survey, nP_metrics):
    tracers = config_survey_update['tracers']
    shape = tuple(len(forecasts_survey[tracer + '_mag_max_eval']) for tracer in tracers)

    total_time       = np.zeros(shape)
    total_efficiency = np.zeros(shape)
    Information_fnl  = np.zeros(shape)
    Information_Da   = np.zeros(shape)
    Information_H    = np.zeros(shape)
    Information_Mnu = np.zeros(shape)
    Information_rsd = np.zeros(shape)

    Information_FnP_metrics_k0001 = np.zeros(shape)
    Information_FnP_metrics_k01 = np.zeros(shape)
    Information_FnP_metrics_k1 = np.zeros(shape)
    

    # Preload fixed config values
    S_survey             = config_survey_update['S_survey']
    N_fibres             = config_survey_update['N_fibres']
    S_FoV                = config_survey_update['S_FoV']
    t_exp                = config_survey_update['exposure_time']
    observational_fraction = config_survey_update['observation_fraction']
    t_full_sky_one_exp   = t_exp * (S_survey / S_FoV)
    n_fibres             = N_fibres / S_FoV

    for idx in itertools.product(*[range(s) for s in shape]):
        tracer_idx = list(zip(tracers, idx))  # [(tracer, i), ...]

        #total_time[idx] = np.sum([config_survey_update[f"{tracer}_fibre_time"][m] for tracer, m in zip(tracers, idx)])
        number_of_pointings = np.sum([
            config_survey_update[tracer + '_target_pointings'][m]
            for tracer, m in tracer_idx
        ])

        comp, passes = survey_design_telescope_metrics.completeness_vs_pass(number_of_pointings, n_fibres)
        total_time[idx] = (
            t_full_sky_one_exp
            * np.array(passes)[np.array(comp) > 0.95][0]
            / (observational_fraction * 365.25 * 24 * 3600)
        )

        total_efficiency[idx] = (
            np.sum([config_survey_update[tracer + '_spec_density'][m]   for tracer, m in tracer_idx])
            / np.sum([config_survey_update[tracer + '_target_density'][m] for tracer, m in tracer_idx])
        )

        Information_fnl[idx] = np.sum([
            1 / forecasts_survey[tracer + '_sigma_fnl_eff'][m] ** 2
            for tracer, m in tracer_idx
        ])
        Information_Da[idx] = np.sum([
            1 / forecasts_survey[tracer + '_sigma_Da_eff'][m] ** 2
            for tracer, m in tracer_idx
        ])
        Information_H[idx] = np.sum([
            1 / forecasts_survey[tracer + '_sigma_H_eff'][m] ** 2
            for tracer, m in tracer_idx
        ])
        Information_Mnu[idx] = np.sum([
            1.0 / forecasts_survey[f"{tracer}_sigma_Mnu_eff"][m]**2
            for tracer, m in zip(tracers, idx)])

        Information_rsd[idx] = np.sum([
            1.0 / forecasts_survey[f"{tracer}_sigma_rsd_eff"][m]**2
            for tracer, m in zip(tracers, idx)])

        Information_FnP_metrics_k0001[idx] = np.sum([
            nP_metrics[f"{tracer}_F_k0.001"][m]
            for tracer, m in zip(tracers, idx)])
        Information_FnP_metrics_k01[idx] = np.sum([
            nP_metrics[f"{tracer}_F_k0.1"][m]
            for tracer, m in zip(tracers, idx)])
        Information_FnP_metrics_k1[idx] = np.sum([
            nP_metrics[f"{tracer}_F_k1"][m]
            for tracer, m in zip(tracers, idx)])
    return {
        'total_survey_time':                 total_time,
        'total_survey_efficiency':           total_efficiency,
        'total_survey_fisher_information_fnl': Information_fnl,
        'total_survey_fisher_information_Da':  Information_Da,
        'total_survey_fisher_information_H':   Information_H,
        'total_survey_fisher_information_Mnu': Information_Mnu,
        'total_survey_fisher_information_rsd': Information_rsd,
        'total_survey_information_FnP_k0.001': Information_FnP_metrics_k0001,
        'total_survey_information_FnP_k0.1': Information_FnP_metrics_k01,
        'total_survey_information_FnP_k1': Information_FnP_metrics_k1,
       # 'total_survey_fisher_nP_ratio_k0.1': ,
       # 'total_survey_fisher_nP_ratio_k1': ,
    }