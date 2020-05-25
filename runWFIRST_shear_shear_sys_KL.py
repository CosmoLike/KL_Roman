#!/home/u17/jiachuanxu/python2_virtualenv/bin/python2.6

import sys
sys.path.append('/home/u17/jiachuanxu/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool

Ntomo=10
Ncl=10
sigmae=0.05
Pco = ""# '', '_DEu95CPL_', '_DEl95CPL_'
strat = "KL" # 'KL', 'WL'

file_source_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_grism_norm")
file_lens_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_clustering_fine_bin_norm")
#data_file = os.path.join(dirname, "datav/WFIRST_KL_shear_shear_opti_SPS002")
data_file = os.path.join(dirname, "datav/WFIRST_%s_shear_shear_opti_grism_Ntomo%2d_Ncl%2d_sigmae%.2f"%(strat,Ntomo,Ncl,sigmae))
#cov_file = os.path.join(dirname, "cov/WFIRST_Tully_Fisher_SN10_sigmae0.05_shear_shear_inv")
cov_file = os.path.join(dirname, "cov/WFIRST_%s_Ntomo%2d_Ncl%2d_sigmae%.2f_shear_shear_inv"%(strat,Ntomo,Ncl,sigmae))
#cov_file = os.path.join(dirname, "cov/WFIRST_KL_Ntomo10_Ncl10_sigmae0.05_KLnorm_shear_shear_inv")
#chain_file = "/extra/jiachuanxu/WFIRST_forecasts/chains/like_WFIRST_KL_SN10_opti_shear_shear_sys_sigmae0.05"
chain_file = "/xdisk/timeifler/mig2020/extra/jiachuanxu/WFIRST_forecasts/chains/like_WFIRST_%s_shear_shear_cos_Ntomo%2d_Ncl%2d_sigmae%.2f_grismzdist%s"%(strat,Ntomo, Ncl, sigmae,Pco)

initcosmo("halofit")
# initbins(Ncl, lmin,    lmax, lmax_shear, Rmin_bias, Ntomo_source, Ntomo_lens)
initbins(  Ncl, 30.0,  4000.0,     4000.0,      21.0,        Ntomo,      Ntomo)
initpriors_KL("photo_opti","shear_opti","none","none")
initsurvey("WFIRST_KL")
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
#print "Dim of param space: ", len(sample_params)
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

sample_main(sample_params,1000,560,1,chain_file+"_1000", blind=False, pool=MPIPool(), KL=True)

