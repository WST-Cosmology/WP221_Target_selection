#bright time: BG_bright
#grey time: BG_Faint, LRG, ELG
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
                 'tracer_N_zm_file' : [path+f'LS_BG_FAINT_WST.npz', path+f'LS_LRG_WST.npz', path+f'LS_ELG_WST.npz'],
                 'tracers' : ['BG_faint', 'LRG', 'ELG'],
                 'color' : ['darkorange', 'brown','peru']}

config_survey_dark = {'survey_type': 'Dark',
                 'N_fibres': 30000,
                      'S_FoV': 3,
                 'S_survey': 18000,
                 'exposure_time': 1000, 
                 'observation_fraction': 0.8 * 0.5 * 0.45,
                 'tracer_N_zm_file' : [path+f'COSMOS_QSO_WST_QSO_no_H.npz', path+f'COSMOS_LBG_udropout_highz.npz',
                                       path+f'COSMOS_LBG_gdropout.npz', path+f'COSMOS_LBG_rdropout.npz'],
                 'tracers' : ['QSO', 'LBGu', 'LBGg', 'LBGr'],
                    'color' : ['k', 'm','g','r']}