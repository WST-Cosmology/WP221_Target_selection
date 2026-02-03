data_path = '/global/cfs/cdirs/desi/users/cpayerne/data_WP221_Target_selection/'
photometric_catalogs = data_path + 'photometric_catalogs/'
tracer_catalogs = data_path + 'tracer_catalogs/'

get_XMMLSS_CLAUDS_HSC = 'curl -O -L "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/clauds/desprez/PublicRelease/XMMLSS-SExtractor-Lephare.fits"'

get_COSMOS_CLAUDS_HSC = 'curl -O -L "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/clauds/desprez/PublicRelease/COSMOS-SExtractor-Lephare.fits"'

get_COSMOS_wise_data = 'curl -O "https://irsa.ipac.caltech.edu/workspace/TMP_AVox9i_30065/Gator/irsa/24256/wise_allwise.allwise_p3as_psd_24256.tbl"'

extract_WISE_cosmos ='ra>148 and ra<153 and dec>0.5 and dec<4.5'

Surface_cosmos_deg2 = 11.26 #deg2
Surface_xmm_deg2 = 9.56 #deg2


#Les étoiles sont bien incluses dans les catalogues, et il y a plusieurs flag pour les identifier, tout depend du type de catalogue (SExtractor ou HSCpipe):
#HSCpipe
#- isCompact —> flag tous les objets PSF-likes (Existe pour toutes les bandes de HSC)
#- isStarTemp —> flag les objets mieux fittés avec une SED d’étoiles que les SED de galaxies
#- isStar —>  la combinaison “&” des deux flags

#SExtractor
#- OBJ_TYPE==2 —> flag les étoiles
#- COMPACT==1 —> flag les objets compacts

