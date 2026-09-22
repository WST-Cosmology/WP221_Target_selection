import numpy as np
import matplotlib.pyplot as plt

import pyccl as ccl
import sys
sys.path.append('../forecasts/')
import fisher_matrix_bao_SuEisenstein

import _tracer_spectroscopic_efficiency
import _survey_design_telescope_metrics
import _surveys
import _survey_design_science_metrics

import pickle
def save_pickle(dat, filename, **kwargs):
    file = open(filename,'wb')
    pickle.dump(dat, file)
    file.close()
def load(filename, **kwargs):
    with open(filename, 'rb') as fin:
        return pickle.load(fin, **kwargs)

h=0.677
deltac=1.686
H0=100*h
c_ls=300*10**3
nlim=10000
n_s=0.968
cosmo = ccl.Cosmology(Omega_c=0.27, Omega_b=0.045, h=h, A_s=2.1e-9, n_s=n_s,transfer_function='boltzmann_camb')

path = '../target_selection/photom_redshift_distribution/'

config_survey_bright_desilike = {'survey_type': 'Bright_desilike',
                 'N_fibres': 30000,
                  'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 180, 
                 'observation_fraction': 0.8 * 0.5 * 0.2,
                 'tracer_N_zm_file' : [path+f'LS_BG_BRIGHT_WST.npz'],
                 'tracers' : ['BG_bright'],
                 'limiting_mag_band': ['r'],
                 'color' : ['gold']}

config_survey_bright_desilike_bao = {'survey_type': 'Bright_desilike_bao',
                 'N_fibres': 30000,
                  'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 180, 
                 'observation_fraction': 0.8 * 0.5 * 0.2,
                 'tracer_N_zm_file' : [path+f'LS_BG_BRIGHT_WST.npz'],
                 'tracers' : ['BG_bright'],
                 'limiting_mag_band': ['r'],
                 'color' : ['gold']}

config_survey_grey_desilike = {'survey_type': 'Grey_desilike',
                 'N_fibres': 30000,
                      'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 1000, 
                 'observation_fraction': 0.8 * 0.5 * 0.35,
                 'tracer_N_zm_file' : [path+f'LS_BG_FAINT_WST.npz', 
                                       path+f'LS_LRG_WST.npz',
                                       path+f'LS_ELG_WST.npz'],
                 'tracers' : ['BG_faint', 'LRG', 'ELG'],
                 'limiting_mag_band': ['r', 'z', 'g'],
                 'color' : ['darkorange', 'brown','peru']}

config_survey_grey_desilike_bao = {'survey_type': 'Grey_desilike_bao',
                 'N_fibres': 30000,
                      'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 1000, 
                 'observation_fraction': 0.8 * 0.5 * 0.35,
                 'tracer_N_zm_file' : [path+f'LS_BG_FAINT_WST.npz', 
                                       path+f'LS_LRG_WST.npz',
                                       path+f'LS_ELG_WST.npz'],
                 'tracers' : ['BG_faint', 'LRG', 'ELG'],
                 'limiting_mag_band': ['r', 'z', 'g'],
                 'color' : ['darkorange', 'brown','peru']}

config_survey_grey_magmax = {'survey_type': 'Grey_MagMax',
                 'N_fibres': 30000,
                      'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 1000, 
                 'observation_fraction': 0.8 * 0.5 * 0.35,
                 'tracer_N_zm_file' : [path+f'COSMOS_H_MagLim_WST.npz',
                                      path+f'COSMOS_H_MagLim_WST.npz',
                                      path+f'COSMOS_H_MagLim_WST.npz',
                                      path+f'COSMOS_H_MagLim_WST.npz'],
                 'tracers' : ['MagMax',
                              'MagMax_lowz', 'MagMax_midz', 'MagMax_highz'],
                 'limiting_mag_band': ['H', 'H', 'H', 'H'],
                 'color' : ['cyan',
                            'b', 'darkblue', 'dodgerblue']}

config_survey_dark_qso_only = {'survey_type': 'Dark_qso_only',
                 'N_fibres': 30000,
                      'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 1000, 
                 'observation_fraction': 0.8 * 0.5 * 0.45,
                 'tracer_N_zm_file' : [None],
                 'tracers' : ['QSO_qlf'],
                 'limiting_mag_band': ['r'],
                    'color' : ['k']}

config_survey_dark_qso_only_bao = {'survey_type': 'Dark_qso_only_bao',
                 'N_fibres': 30000,
                      'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 1000, 
                 'observation_fraction': 0.8 * 0.5 * 0.45,
                 'tracer_N_zm_file' : [None],
                 'tracers' : ['QSO_qlf'],
                 'limiting_mag_band': ['r'],
                    'color' : ['k']}

config_survey_dark_lbg_only = {'survey_type': 'Dark_lbg_only',
                 'N_fibres': 30000,
                      'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 1000, 
                 'observation_fraction': 0.8 * 0.5 * 0.45,
                 'tracer_N_zm_file' : [path+f'COSMOS_LBG_udropout_highz.npz',
                                       path+f'COSMOS_LBG_gdropout.npz', 
                                       path+f'COSMOS_LBG_rdropout.npz'],
                 'tracers' : ['LBGu', 'LBGg', 'LBGr'],
                 'limiting_mag_band': ['r', 'i', 'z'],
                    'color' : ['m','g','r']}

mag_max_eval_range = {'Bright_desilike': [[19, 21]],
                      'Grey_desilike'  : [[21, 22], [20, 23], [23,25]],
                      'Bright_desilike_bao': [[19, 21]],
                      'Grey_desilike_bao'  : [[21, 22], [20, 23], [23,25]],
                      'Grey_MagMax': [[19, 22], 
                                      [19, 22], [19, 22], [19, 22]], 
                      'Dark_lbg_only'  : [[24.2, 26], [24.2, 26], [24.2, 26]],
                     'Dark_qso_only'  : [[23, 25]],
                     'Dark_qso_only_bao'  : [[23, 25]]}

redshift_eval_range = {'Bright_desilike': [[0, 1.5]],
                      'Grey_desilike'  :  [[0, 2], [0, 2], [0, 2]],
                      'Bright_desilike_bao': [[0, 0.4]],
                      'Grey_desilike_bao'  :  [[0.4, 0.7], [0.7, 1.2], [1.2, 1.6]],
                      'Grey_MagMax': [[0, 1.5], [0., 0.5], [0.5, 1.], [1., 1.5]], 
                      'Dark_lbg_only'  :  [[2, 4.5], [2.5, 5.5], [4, 6]],
                      'Dark_qso_only'  : [[0.8, 3.1]],
                      'Dark_qso_only_bao'  : [[1.6, 3.1]]}

multi_mag_bin_approach = {'Bright_desilike': [False],
                          'Grey_desilike'  :  [False, False, False, False],
                          'Bright_desilike_bao': [False],
                          'Grey_desilike_bao'  :  [False, False, False, False],
                          'Grey_MagMax': [False,False,False,False], 
                          'Dark_lbg_only'  : [True, True, True],
                         'Dark_qso_only'  : [False],
                         'Dark_qso_only_bao'  : [False]}

config_surveys = [config_survey_grey_desilike_bao]


for i, config_survey in enumerate(config_surveys):

    survey = config_survey['survey_type']

    #if survey != 'Dark': continue

    config_survey_update = _survey_design_telescope_metrics.Survey_design_telescope_metrics(config_survey,mag_max_eval_range=mag_max_eval_range[survey], max_mag=None)   

    per_tracer_forecasts = _survey_design_science_metrics.Survey_design_science_metrics(config_survey_update, cosmo, 
                                                                               redshift_eval_range =redshift_eval_range[survey], 
                                                                               mag_max_eval_range=mag_max_eval_range[survey],
                                                                                multi_mag_bin_approach=multi_mag_bin_approach[survey],
                                                                                       which_param=['bao','rsd'])

    #FnP_tracer = _survey_design_science_metrics.Survey_design_nP_metrics(config_survey_update, cosmo, 
    #                                                                           redshift_eval_range =redshift_eval_range[survey], 
    #                                                                           mag_max_eval_range=mag_max_eval_range[survey],
    #                                                                            multi_mag_bin_approach=multi_mag_bin_approach[survey])

    total_Informations = _survey_design_science_metrics.build_total_survey_information_metrics(config_survey_update, per_tracer_forecasts, which_param=['bao','rsd'])

    file_to_save = {}
    file_to_save['config_survey'] = config_survey_update
    file_to_save['per_tracer_forecasts'] = per_tracer_forecasts
    file_to_save['total_survey_Informations'] = total_Informations

    save_pickle(file_to_save, './telescope_and_science_metrics/'+'survey_design_' + survey + '.pkl')

    