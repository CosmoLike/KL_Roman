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
Ntomo_src, Ntomo_lens = 4, 4
Ncl = 15
ell_min, ell_max, ell_max_shear = 20.0, 3000.0, 3000.0
Rmin_bias = 21.0
strat = "SKA_KL"

## external prior, e.g. "Planck15_BAO_H070p6_JLA_w0wa"
## default is "none"
external_prior = "none"    

#sigmae_list = np.array([0.02, 0.04, 0.06, 0.10, 0.20, 0.30])*np.sqrt(2)
#Nsrc_list = np.array([0.4761, 0.1629, 0.1553, 0.0881, 0.1006, 0.0740])
nz_src_files = "zdistris/zdistri_trecs_KL"
nz_lens_file = "zdistris/lens_LSSTY1"
data_vector_file = "datav/SKA_KL_shear_shear_Ntomo%d_Ncl%d_dmo_OneComp"
invcovmat_file = "invcov/SKA_KL_ssss_invcov_Ncl%d_Ntomo%d_OneComp"
# chain_output_file = "chains/SKA_KL_OmS8_ss_Ncl%d_Ntomo%d_OneComp"
chain_output_file = "chains/SKA_KL_LCDM_ss_Ncl%d_Ntomo%d_OneComp"

## flag
DE_FLAG = False
KL_FLAG = True         # true if perform KL forecast
one = True              # one component
photoz_flag = False     # enable different sigma_photoz senario

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
initpriors_IA_bary("spec_SKA_KL", "shear_SKA_KL", "none", external_prior,
    False, 3.0, 1.2, 3.8, 2.0, 
    False, 16, 1.9, 0.7)
initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initclusters()
initia("none","GAMA")
initprobes("shear_shear")
initdatainv(cov_file ,data_file)

# only sample two parameters
sample_params = ['omega_m','sigma_8']
# sample_params=sample_LCDM_only()
# sample_params= sample_cosmology_only()

# The `sample_main` function is being called with several parameters:
sample_main(sample_params, nsteps, nwalkers, nthreads, chain_file+"_%d"%(nsteps), blind=False, pool=MPIPool(), KL=KL_FLAG, one=one, photoz_flag=photoz_flag)

