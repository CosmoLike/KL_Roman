import sys
sys.path.append('/home/u17/jiachuanxu/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('iSelection', type=int,
                    help='Index of target selection scheme')
parser.add_argument('iSN', type=int,
                    help='Index of shape noise scheme')
args = parser.parse_args()
#########################################################
dirname = "/home/u17/jiachuanxu/CosmoLike/KL_WFIRST"
outdirname = "/xdisk/timeifler/jiachuanxu/DESI2KL"
Ntomo_src, Ntomo_lens = 4, 10
Ncl = 15
ell_min, ell_max, ell_max_shear = 20.0, 3000.0, 3000.0
Rmin_bias = 21.0
strat = "DESI2_KL_%d%d"%(args.iSelection, args.iSN)
nz_src_files = [
	"zdistris/zdistri_DESI2_KL_LS_DR9_sample1_v2",
    "zdistris/zdistri_DESI2_KL_LS_DR9_sample2_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Any_sample1_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Any_sample2_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Bright_sample1_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Bright_sample2_v2",
]
nz_lens_file = "zdistris/zdistri_WFIRST_LSST_clustering_fine_bin_norm"
data_vector_file = "datav/DESI2_KL_%d0_shear_shear_Ntomo%d_Ncl%d_dmo"
invcovmat_file = "invcov/DESI2_KL_v3_%d%d_ssss_invcov_Ncl%d_Ntomo%d"
baryon_PCS_file = "datav/DESI2_KL_%d%d_shear_shear_Ntomo%d_Ncl%d_9sim.pca"
#chain_output_file = "chains/DESI2_KL_v3_PlanckBAOJLA_%d%d_ss_Ncl%d_Ntomo%d"
chain_output_file = "chains/DESI2_KL_v3_s8split_only_zlow015__%d%d_ss_Ncl%d_Ntomo%d"
#chain_output_file = "chains_test/DESI2_KL_s8split_only_zlow015_%d%d_ss_Ncl%d_Ntomo%d"
#external_prior = "Planck15_BAO_H070p6_JLA_w0wa" # default: "none"
external_prior = "none"
NPCs_used = 2
#cosmo_model = "LCDM_split"
cosmo_model = "s8split_only"
runmode = "halofit_split"
############################################################
file_source_z = os.path.join(dirname, nz_src_files[args.iSelection])
file_lens_z = os.path.join(dirname, nz_lens_file)
data_file = os.path.join(dirname, 
    data_vector_file%(args.iSelection, Ntomo_src, Ncl))
cov_file = os.path.join(outdirname, 
    invcovmat_file%(args.iSelection, args.iSN, Ncl, Ntomo_src))
bary_file = os.path.join(dirname, 
    baryon_PCS_file%(args.iSelection, args.iSN, Ntomo_src, Ncl))
chain_file = os.path.join(outdirname, 
    chain_output_file%(args.iSelection, args.iSN, Ncl, Ntomo_src))

initcosmo(runmode)
initbins(Ncl,ell_min,ell_max,ell_max_shear,Rmin_bias,Ntomo_src,Ntomo_lens)
initpriors_IA_bary("spec_DESI2", "shear_KL_DESI2", "none", external_prior,
    False, 3.0, 1.2, 3.8, 2.0, 
    True, 20.0, 6.0, 2.0)
initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
#initclusters()
initia("none","none")
initprobes("shear_shear")
initdatainvbary(cov_file ,data_file, bary_file)

sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear(), 
    MG=False, NPCs=NPCs_used, cosmology=cosmo_model, source_photo_z=True, 
    shear_calibration=True, IA=False)

sample_main(sample_params,5000,400,1,chain_file+"_5000", blind=False, pool=MPIPool(), KL=True)

