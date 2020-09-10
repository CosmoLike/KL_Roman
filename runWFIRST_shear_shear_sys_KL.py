#!/home/u17/jiachuanxu/python2_virtualenv/bin/python2.6

import sys
sys.path.append('/home/u17/jiachuanxu/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool

Ntomo=30
Ncl=10
sigmae=0.08
Pco=""
strat="KL"

file_source_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_grism_norm")
file_lens_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_clustering_fine_bin_norm")
#data_file = os.path.join(dirname, "datav/WFIRST_KL_shear_shear_opti_SPS002")
data_file = os.path.join(dirname, "datav/WFIRST_KL_shear_shear_opti_Ntomo30_Ncl10_zdistgrism_dmo")
#cov_file = os.path.join(dirname, "cov/WFIRST_Tully_Fisher_SN10_sigmae0.05_shear_shear_inv")
cov_file = os.path.join(dirname, "cov/WFIRST_KL_Ntomo30_Ncl10_sigmae0.08_grismzdist_shear_shear_inv")
#chain_file = "/extra/jiachuanxu/WFIRST_forecasts/chains/like_WFIRST_KL_SN10_opti_shear_shear_sys_sigmae0.05"
bary_file = "datav/WFIRST_KL_shear_shear_opti_Ntomo30_Ncl10_zdistgrism_PCs"
chain_file = "/xdisk/timeifler/mig2020/extra/jiachuanxu/WFIRST_forecasts/chains/like_WFIRST_KL_shear_shear_cos_Ntomo30_Ncl10_sigmae0.08_grismzdist_bary"

initcosmo("halofit")
# initbins(Ncl, lmin,    lmax, lmax_shear, Rmin_bias, Ntomo_source, Ntomo_lens)
initbins( Ncl, 30.0,    4000.0,     4000.0,      21.0,           Ntomo,         Ntomo)
initpriors_IA_bary("photo_KL","shear_KL","none","none", False, 3.0, 1.2, 3.8, 2.0, True, 16.0, 5.0, 0.8)
initsurvey("WFIRST_KL")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initclusters()
initia("none","none")

# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("shear_shear")
initdatainvbary(cov_file ,data_file, bary_file)

#sample_params=sample_LCDM_only()
#sample_params= sample_cosmology_only()
sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
sample_params += ['bary_%d'%i for i in xrange(3)]
#print "Dim of param space: ", len(sample_params)
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

sample_main(sample_params,5000,700,1,chain_file+"_5000", blind=False, pool=MPIPool(), KL=True)

