#!/home/u17/jiachuanxu/python2_virtualenv/bin/python2.6

import sys
sys.path.append('/home/u17/jiachuanxu/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool

file_source_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin_norm")
file_lens_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_clustering_fine_bin_norm")
data_file = os.path.join(dirname, "datav/WFIRST_shear_shear_opti")
#cov_file = os.path.join(dirname, "cov/WFIRST+LSST_SN10_opti_shear_shear_inv")
cov_file = os.path.join(dirname, "cov/WFIRST_WL_SN10_shear_shear_inv")
chain_file = "/extra/jiachuanxu/WFIRST_forecasts/chains/like_WFIRST_WL_SN10_opti_shear_shear_cos_opti"

initcosmo("halofit")
# initbins(Ncl, lmin,    lmax, lmax_shear, Rmin_bias, Ntomo_source, Ntomo_lens)
initbins( 20, 30.0, 4000.0,     4000.0,      21.0,           10,         10)
initpriors("photo_opti","shear_opti","none","none")
initsurvey("WFIRST_WL")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initclusters()
initia("none","none")

# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("shear_shear")
initdatainv(cov_file ,data_file)

#sample_params=sample_LCDM_only()
sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

#sample_main(sample_params,10000,560,1,chain_file, blind=False, pool=MPIPool())
sample_main(sample_params,1000,560,1,chain_file+"_5000", blind=False, pool=MPIPool())

