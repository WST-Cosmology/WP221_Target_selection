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

config_survey_bright = {'survey_type': 'Bright',
                 'N_fibres': 30000,
                  'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 180, 
                 'observation_fraction': 0.8 * 0.5 * 0.2,
                 'tracer_N_zm_file' : [path+f'LS_BG_BRIGHT_WST.npz'],
                 'tracers' : ['BG_bright'],
                 'color' : ['gold']}
config_survey_grey = {'survey_type': 'Grey',
                 'N_fibres': 30000,
                      'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 1000, 
                 'observation_fraction': 0.8 * 0.5 * 0.35,
                 'tracer_N_zm_file' : [path+f'LS_BG_FAINT_WST.npz', 
                                       path+f'LS_LRG_WST.npz',
                                       path+f'LS_ELG_WST.npz'],
                 'tracers' : ['BG_faint', 'LRG', 'ELG'],
                 'color' : ['darkorange', 'brown','peru']}

config_survey_dark = {'survey_type': 'Dark',
                 'N_fibres': 30000,
                      'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 1000, 
                 'observation_fraction': 0.8 * 0.5 * 0.45,
                 'tracer_N_zm_file' : [#path+f'COSMOS_QSO_WST_QSO_no_H.npz', 
                                       path+f'COSMOS_LBG_udropout_highz.npz',
                                       path+f'COSMOS_LBG_gdropout.npz', 
                                       path+f'COSMOS_LBG_rdropout.npz'],
                 'tracers' : [#'QSO', 
                              'LBGu', 'LBGg', 'LBGr'],
                    'color' : [#'k', 
                               'm','g','r']}

mag_max_eval_range = {'Bright': [[19, 21]],
                      'Grey'  : [[21, 22], [20, 23], [23,25]],
                      'Dark'  : [
                          #[23, 24.5], 
                          [24.2, 26], [24.2, 26], [24.2, 26]]}

redshift_eval_range = {'Bright': [[0, 1.5]],
                      'Grey'  :  [[0, 2], [0, 2], [0, 2]],
                      'Dark'  :  [[0, 4], 
                                  [2, 4.5], [2.5, 5], [4, 6]]}

multi_mag_bin_approach = {'Bright': [False],
                          'Grey'  :  [False, False,False],
                          'Dark'  : [#False, 
                                     True, True, True]}

config_surveys = [config_survey_bright, 
                  config_survey_grey, 
                  config_survey_dark
                 ]

for i, config_survey in enumerate(config_surveys):

    survey = config_survey['survey_type']

    config_survey_update = _survey_design_telescope_metrics.Survey_design_telescope_metrics(config_survey,mag_max_eval_range=mag_max_eval_range[survey], max_mag=None)   

    per_tracer_forecasts = _survey_design_science_metrics.Survey_design_science_metrics(config_survey_update, cosmo, 
                                                                               redshift_eval_range =redshift_eval_range[survey], 
                                                                               mag_max_eval_range=mag_max_eval_range[survey],
                                                                                multi_mag_bin_approach=multi_mag_bin_approach[survey])

    FnP_tracer = _survey_design_science_metrics.Survey_design_nP_metrics(config_survey_update, cosmo, 
                                                                               redshift_eval_range =redshift_eval_range[survey], 
                                                                               mag_max_eval_range=mag_max_eval_range[survey],
                                                                                multi_mag_bin_approach=multi_mag_bin_approach[survey])

    total_Informations = _survey_design_science_metrics.build_total_survey_information_metrics(config_survey_update, per_tracer_forecasts, FnP_tracer)

    file_to_save = {}
    file_to_save['config_survey'] = config_survey_update
    file_to_save['per_tracer_forecasts'] = per_tracer_forecasts
    file_to_save['per_tracer_FnP'] = FnP_tracer
    file_to_save['total_survey_Informations'] = total_Informations

    save_pickle(file_to_save, './telescope_and_science_metrics/'+'survey_design_' + survey + '.pkl')

    