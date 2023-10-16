import sys
sys.path.append('/home/u17/jiachuanxu/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool

#########################################################
dirname = "/home/u17/jiachuanxu/CosmoLike/KL_WFIRST"
outdirname = "/xdisk/timeifler/jiachuanxu/DESI2KL"
Ntomo_src, Ntomo_lens = 10, 10
Ncl = 15
ell_min, ell_max, ell_max_shear = 20.0, 3000.0, 3000.0
Rmin_bias = 21.0
strat = "LSST_Y10"
nz_src_file = "zdistris/src_LSSTY10"
nz_lens_file = "zdistris/lens_LSSTY10"
data_vector_file = "datav/LSST_Y10_shear_shear_Ntomo%d_Ncl%d_dmo"
invcovmat_file = "invcov/LSST_Y10_ssss_invcov_Ncl%d_Ntomo%d"
baryon_PCS_file = "datav/LSST_Y10_shear_shear_Ntomo%d_Ncl%d_9sim.pca"
#chain_output_file = "chains/LSST_Y10_ss_Ncl%d_Ntomo%d"
chain_output_file = "chains/LSST_Y10_LCDM_ss_Ncl%d_Ntomo%d"
external_probe = "none"
DE_FLAG = False
NPCs_used = 2
############################################################
file_source_z = os.path.join(dirname, nz_src_file)
file_lens_z = os.path.join(dirname, nz_lens_file)
data_file = os.path.join(dirname, data_vector_file%(Ntomo_src, Ncl))
cov_file = os.path.join(outdirname, invcovmat_file%(Ncl, Ntomo_src))
bary_file = os.path.join(dirname, baryon_PCS_file%(Ntomo_src, Ncl))
chain_file = os.path.join(outdirname, chain_output_file%(Ncl, Ntomo_src))

initcosmo("halofit")
initbins(Ncl,ell_min,ell_max,ell_max_shear,Rmin_bias,Ntomo_src,Ntomo_lens)
initpriors_IA_bary("photo_LSST_Y10","shear_LSST_Y10","none",external_probe,
    True, 3.0,1.2,3.8,2.0,
    True, 40.0,10.0,0.8)
initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initclusters()
initia("NLA_HF","GAMA")
initprobes("shear_shear")
initdatainvbary(cov_file ,data_file, bary_file)

#sample_params=sample_LCDM_only()
#sample_params= sample_cosmology_only()
sample_params = sample_cosmology_shear_nuisance_IA_bary(get_N_tomo_shear(), 
	NPCs=NPCs_used, DE=DE_FLAG)
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

sample_main(sample_params,8000,400,1,chain_file+"_8000", blind=False, pool=MPIPool())
