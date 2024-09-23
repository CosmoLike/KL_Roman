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
strat = "SKA_WL" 

## flag
DE_FLAG = False     # dynamical dark energy
KL_FLAG = False     # true if perform KL forecast
one = False         # one component
photoz_flag = False # enable different sigma_photoz senario

## sample parameters
# sample_params = ['omega_m','sigma_8']
# sample_params = sample_cosmology_only()
sample_params = sample_LCDM_only()

## directory and file names
nz_src_files = "zdistris/zdistri_trecs_WL"
nz_lens_file = "zdistris/lens_LSSTY1"
data_vector_file = "datav/%s_shear_shear_Ntomo%d_Ncl%d_dmo"
invcovmat_file = "invcov/%s_ssss_invcov_Ncl%d_Ntomo%d"
# chain_output_file = "chains/%s_OmS8_ss_Ncl%d_Ntomo%d"
chain_output_file = "chains/%s_LCDM_ss_Ncl%d_Ntomo%d"

## external prior, e.g. "Planck15_BAO_H070p6_JLA_w0wa"
## default is "none"
external_prior = "none"   

## mcmc setting
nsteps = 1000
nwalkers = 400
nthreads = 1

############################################################
file_source_z = os.path.join(dirname, nz_src_files)
file_lens_z = os.path.join(dirname, nz_lens_file)
data_file = os.path.join(dirname, data_vector_file%(strat, Ntomo_src, Ncl))
cov_file = os.path.join(outdirname, invcovmat_file%(strat, Ncl, Ntomo_src))
chain_file = os.path.join(outdirname, chain_output_file%(strat, Ncl, Ntomo_src))

initcosmo("halofit")
initbins(Ncl,ell_min,ell_max,ell_max_shear,Rmin_bias,Ntomo_src,Ntomo_lens)
initpriors_IA_bary("spec_%s"%(strat), "shear_%s"%(strat), "none", external_prior,
    True, 3.0, 1.2, 3.8, 2.0, 
    False, 16, 1.9, 0.7)
initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initclusters()
initia("NLA_HF","GAMA")
initprobes("shear_shear")
initdatainv(cov_file ,data_file)

sample_main(sample_params, nsteps, nwalkers, nthreads, chain_file+"_%d"%(nsteps), 
    pool=MPIPool(), blind=False, KL=KL_FLAG, one=one, photoz_flag=photoz_flag)

