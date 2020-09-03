#!/home/u17/jiachuanxu/python2_virtualenv/bin/python2.6

import sys
sys.path.append('/home/u17/jiachuanxu/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool

file_source_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin_norm")
file_lens_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_clustering_fine_bin_norm")
data_file = os.path.join(dirname, "datav/WFIRST_shear_shear_opti_Ntomo10_Ncl20_sigmae0.37_dmo_ia") # vanilla WL
#data_file = os.path.join(dirname, "datav/WFIRST_shear_shear_opti")
#data_file = os.path.join(dirname, "datav/WFIRST_shear_shear_opti_Ntomo10_Ncl20_sigmae0.37_dmo_ia")
#cov_file = os.path.join(dirname, "cov/WFIRST+LSST_SN10_opti_shear_shear_inv")
cov_file = os.path.join(dirname, "cov/WFIRST_WL_Ntomo10_Ncl20_sigmae0.37_baseline_shear_shear_inv")
#cov_file = os.path.join(dirname, "cov/WFIRST_WL_Ntomo10_Ncl20_sigmae0.37_IA_shear_shear_inv")
bary_file = os.path.join(dirname, "datav/WFIRST_shear_shear_opti_Ntomo10_Ncl20_sigmae0.37_ia")
chain_file = "/xdisk/timeifler/mig2020/extra/jiachuanxu/WFIRST_forecasts/chains/like_WFIRST_WL_SN10_opti_shear_shear_sys_opti_IA_bary"

initcosmo("halofit")
# initbins(Ncl, lmin,    lmax, lmax_shear, Rmin_bias, Ntomo_source, Ntomo_lens)
initbins( 20, 30.0, 4000.0,     4000.0,      21.0,           10,         10)
#initpriors("photo_opti","shear_opti","none","none")
initpriors_IA_bary("photo_opti","shear_opti","none","none",
    3.0,1.2,3.8,2.0,# prior std for IA: A_ia, beta_ia, eta_ia, eta_highz, if gaussian
    16.0,5.0,0.8)# prior std for baryon PCs: Q1, Q2, Q3, if gaussian
#initpriors_IA("photo_opti","shear_opti","none","none",
#    3.0,1.2,3.8,2.0)# prior std for IA: A_ia, beta_ia, eta_ia, eta_highz, if gaussian
initsurvey("WFIRST_WL")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initclusters()
initia("NLA_HF","GAMA")
#initia("none","GAMA")

# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("shear_shear")
initdatainvbary(cov_file ,data_file, bary_file)
#initdatainv(cov_file, data_file)

#sample_params=sample_LCDM_only()
#sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance_IA(get_N_tomo_shear())
sample_params = sample_cosmology_shear_nuisance_IA_bary(get_N_tomo_shear()) # 35 sampled params
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

#sample_main(sample_params,10000,560,1,chain_file, blind=False, pool=MPIPool())
sample_main(sample_params,5000,700,1,chain_file+"_5000", blind=False, pool=MPIPool())
"""
print "Sampled Parameters: ",sample_params
starting_point = InputCosmologyParams.fiducial().convert_to_vector_filter(sample_params)
starting_point += InputNuisanceParams().fiducial().convert_to_vector_filter(sample_params)
starting_point += InputNuisanceParamsGRS().fiducial().convert_to_vector_filter(sample_params)
std = InputCosmologyParams.fiducial_sigma().convert_to_vector_filter(sample_params)
std += InputNuisanceParams().fiducial_sigma().convert_to_vector_filter(sample_params)
std += InputNuisanceParamsGRS().fiducial_sigma().convert_to_vector_filter(sample_params)
likelihood = LikelihoodFunctionWrapper(sample_params)
print "FIDUCIAL COSMOLOGY: ",starting_point
deviated_point = [starting_point[i]+0.1*std[i] for i in range(len(starting_point))]
print "likelihood at fiducial point: %.4f"%(likelihood(starting_point))
print "likelihood around fiducial point: %.4f"%(likelihood(deviated_point))

test_point = [0.249146116691,0.759908403828,0.87177318073,-0.872807021823,0.109406808791,0.0603944633065,0.830207263835]
test_point += InputNuisanceParams().fiducial().convert_to_vector_filter(sample_params)
test_point += InputNuisanceParamsGRS().fiducial().convert_to_vector_filter(sample_params)
print "likelihood at test point: %.4f"%(likelihood(test_point))
"""
