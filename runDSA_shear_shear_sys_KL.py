import sys
sys.path.append('/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool
# from argparse import ArgumentParser

# parser = ArgumentParser()

# parser.add_argument('iSelection', type=int,
#                     help='Index of target selection scheme')
# parser.add_argument('iSN', type=int,
#                     help='Index of shape noise scheme')
# args = parser.parse_args()
#########################################################
dirname = "/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST"
outdirname = "/home/u15/yhhuang/cosmology/dsa"
Ntomo_src, Ntomo_lens = 4, 10
Ncl = 15
ell_min, ell_max, ell_max_shear = 20.0, 3000.0, 3000.0
Rmin_bias = 21.0
strat = "DSA_allsky"
#external_prior = "Planck15_BAO_H070p6_JLA_w0wa" # default: "none"
external_prior = "none"

#sigmae_list = np.array([0.02, 0.04, 0.06, 0.10, 0.20, 0.30])*np.sqrt(2)
#Nsrc_list = np.array([0.4761, 0.1629, 0.1553, 0.0881, 0.1006, 0.0740])
nz_src_files = "zdistris/zdistri_DSA_allsky"
nz_lens_file = "zdistris/lens_LSSTY1"
data_vector_file = "datav/DSA_allsky_shear_shear_Ntomo%d_Ncl%d_dmo_Om40_test"
invcovmat_file = "invcov/DSA_allsky_ssss_invcov_Ncl%d_Ntomo%d"
chain_output_file = "chains/DSA_allsky_LCDM_ss_Ncl%d_Ntomo%d_Om40_test"
DE_FLAG = False
############################################################
file_source_z = os.path.join(dirname, nz_src_files)
file_lens_z = os.path.join(dirname, nz_lens_file)
data_file = os.path.join(dirname, data_vector_file%(Ntomo_src, Ncl))
cov_file = os.path.join(outdirname, invcovmat_file%(Ncl, Ntomo_src))
chain_file = os.path.join(outdirname, chain_output_file%(Ncl, Ntomo_src))

initcosmo("halofit")
initbins(Ncl,ell_min,ell_max,ell_max_shear,Rmin_bias,Ntomo_src,Ntomo_lens)
#initpriors_KL("photo_opti","shear_opti","none","none")
initpriors_IA_bary("spec_DSA_allsky", "shear_DSA_allsky", "none", external_prior,
    False, 3.0, 1.2, 3.8, 2.0, 
    False, 16, 1.9, 0.7)
initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initclusters()
initia("none","GAMA")
initprobes("shear_shear")
initdatainv(cov_file ,data_file)

#sample_params=sample_LCDM_only()
#sample_params= sample_cosmology_only()
sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear(), DE=DE_FLAG)
# Fix Q3, not constraining that
#sample_params += ['bary_%d'%i for i in xrange(2)]
#print "Dim of param space: ", len(sample_params)
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

sample_main(sample_params,5000,400,1,chain_file+"_5000", blind=False, pool=MPIPool(), KL=True)

