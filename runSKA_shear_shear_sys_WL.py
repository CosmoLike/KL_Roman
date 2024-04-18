import sys
sys.path.append('/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
# from schwimmbad import MPIPool
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
Ntomo_src, Ntomo_lens = 10, 10
Ncl = 15
ell_min, ell_max, ell_max_shear = 20.0, 3000.0, 3000.0
Rmin_bias = 21.0
strat = "SKA_WL"

## external prior, e.g. "Planck15_BAO_H070p6_JLA_w0wa"
## default is "none"
external_prior = "none"    

#sigmae_list = np.array([0.02, 0.04, 0.06, 0.10, 0.20, 0.30])*np.sqrt(2)
#Nsrc_list = np.array([0.4761, 0.1629, 0.1553, 0.0881, 0.1006, 0.0740])
nz_src_files = "zdistris/zdistri_SKA"
nz_lens_file = "zdistris/lens_LSSTY1"
data_vector_file = "datav/SKA_WL_shear_shear_Ntomo%d_Ncl%d_dmo"
invcovmat_file = "invcov/SKA_WL_ssss_invcov_Ncl%d_Ntomo%d"
chain_output_file = "chains/SKA_WL_LCDM_ss_Ncl%d_Ntomo%d"

## flag
DE_FLAG = False
KL_FLAG = False     # true if perform KL forecast
one = False         # one component
photoz_flag = True  # enable different sigma_photoz senario

## mcmc setting
nsteps = 1000
nwalkers = 400
nthreads = 1

############################################################
file_source_z = os.path.join(dirname, nz_src_files)
file_lens_z = os.path.join(dirname, nz_lens_file)
data_file = os.path.join(dirname, data_vector_file%(Ntomo_src, Ncl))
cov_file = os.path.join(outdirname, invcovmat_file%(Ncl, Ntomo_src))
chain_file = os.path.join(outdirname, chain_output_file%(Ncl, Ntomo_src))

initcosmo("halofit")
initbins(Ncl,ell_min,ell_max,ell_max_shear,Rmin_bias,Ntomo_src,Ntomo_lens)

initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initpriors_IA_bary("spec_SKA_WL", "shear_SKA_WL", "none", external_prior,
    True, 3.0, 1.2, 3.8, 2.0, 
    False, 16, 1.9, 0.7)
initclusters()
initia("NLA_HF","GAMA")
initprobes("shear_shear")
initdatainv(cov_file ,data_file)

# sample_params=sample_LCDM_only()
# sample_params= sample_cosmology_onlxy()

# Fix Q3, not constraining that
#sample_params += ['bary_%d'%i for i in xrange(2)]
#print "Dim of param space: ", len(sample_params)
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

# only sample two parameters
# sample_params = ['omega_m','sigma_8']

# The `sample_main` function is being called with several parameters:
# sample_main(sample_params, nsteps, nwalkers, nthreads, chain_file+"_%d"%(nsteps), blind=False, pool=MPIPool(), KL=KL_FLAG, one=one, photoz_flag=photoz_flag)
# sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear(), DE=DE_FLAG)
# sample_main(sample_params,5000,400,1,chain_file+"_5000", blind=False, pool=MPIPool(), KL=True)

