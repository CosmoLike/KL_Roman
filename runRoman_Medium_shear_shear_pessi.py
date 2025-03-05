import sys
sys.path.append('/home/u17/jiachuanxu/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-nsteps', type=int, default=2000, help='MCMC steps')
parser.add_argument('-nwalkers', type=int, default=400, help='N walkers')
args = parser.parse_args()
print "MCMC sampler steps = %d, walkers = %d"%(args.nsteps, args.nwalkers)

ell_max = 4000.0
neff = 41.3
area = 2415
print "ell max = %.0f; neff = %.0f; area = %.0f;"%(ell_max, neff, area)

#########################################################
dirname = "/home/u17/jiachuanxu/CosmoLike/KL_WFIRST"
outdirname = "/xdisk/timeifler/jiachuanxu/RomanPIT"
Ntomo_src, Ntomo_lens = 10, 10
Ncl = 15
ell_min, ell_max_shear = 20.0, 4000.0
Rmin_bias = 21.0
strat = "Roman_Medium"
#print "Survey strat = ", strat
nz_src_file = "n_eff40.nz"
nz_lens_file = "lens_LSSTY1"
data_vector_file = "Roman_WL_403_shear_shear_Ntomo10_Ncl15_Haley_dndz_dmo"
invcovmat_file = "Roman_Medium_ssss_invcov_Ncl15_Ntomo10_Haley_dndz"
baryon_PCS_file = "Roman_Medium_shear_shear_Ntomo10_Ncl15_9sim.pca"
#external_probe = "better_DESI_BAO"
external_probe = "none"
chain_output_file = "Roman_Medium_ss_Ncl15_Ntomo10_pessi_IAbary"

#cosmo_model = "LCDM_split"
#cosmo_model = "s8split_only"
#cosmo_model = "OmS8"
#cosmo_model = "w0wa"
cosmo_model = "w0waCDM"
#runmode = "halofit_split"
runmode = "halofit"
############################################################
file_source_z = os.path.join(dirname, "zdistris", nz_src_file)
file_lens_z = os.path.join(dirname, "zdistris", nz_lens_file)
data_file = os.path.join(dirname, "datav", data_vector_file)
cov_file = os.path.join(outdirname, "invcov", invcovmat_file)
bary_file = os.path.join(dirname, "datav", baryon_PCS_file)
chain_file = os.path.join(outdirname, "chains", chain_output_file)

initcosmo(runmode)
initbins(Ncl,ell_min,ell_max,ell_max_shear,Rmin_bias,Ntomo_src,Ntomo_lens)
initia("none","GAMA")
initpriors_IA_bary("photo_pessi","shear_pessi","none",external_probe,
    True, 3.0,1.2,3.8,2.0,
    True, 80,10.0,0.8)
initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
#initclusters()
initprobes("shear_shear")
initdatainvbary(cov_file, data_file, bary_file)

sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear(), 
    MG=False, NPCs=1, cosmology=cosmo_model, source_photo_z=False, 
    shear_calibration=False, IA=True)
#sample_params = ['omega_m', ] + sample_params 

### test likelihood evaluation bias
#if samp_bary:
#    parval = [0.3156, 0.831, 0.831, 0.0]
#else:
#    parval = [0.3156, 0.831, 0.831]
#test_logpost = test_likelihood(sample_params, parval)
#print "test likelihood:", test_logpost

### run mcmc chains
sample_main(sample_params,args.nsteps,args.nwalkers,1,chain_file+"_%d"%args.nsteps, blind=False, pool=MPIPool())
