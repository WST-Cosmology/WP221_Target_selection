data_path = '/global/cfs/cdirs/desi/users/cpayerne/data_WP221_Target_selection/'

get_XMMLSS_CLAUDS_HSC = 'curl -O -L "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/clauds/desprez/PublicRelease/XMMLSS-SExtractor-Lephare.fits"'

get_COSMOS_CLAUDS_HSC = 'curl -O -L "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/clauds/desprez/PublicRelease/COSMOS-SExtractor-Lephare.fits"'

get_COSMOS_wise_data 'curl -O "https://irsa.ipac.caltech.edu/workspace/TMP_AVox9i_30065/Gator/irsa/24256/wise_allwise.allwise_p3as_psd_24256.tbl"'

extract_WISE_cosmos ='ra>148 and ra<153 and dec>0.5 and dec<4.5'

