#!/home/u17/jiachuanxu/python2_virtualenv/bin/python2.6

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

#sigmae_list = np.array([0.02, 0.04, 0.06, 0.10, 0.20, 0.30])*np.sqrt(2)
#Nsrc_list = np.array([0.4761, 0.1629, 0.1553, 0.0881, 0.1006, 0.0740])
nz_src_files = [
	"zdistris/zdistri_DESI2_KL_LS_DR9_sample1_v2",
    "zdistris/zdistri_DESI2_KL_LS_DR9_sample2_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Any_sample1_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Any_sample2_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Bright_sample1_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Bright_sample2_v2",
]
nz_lens_file = "zdistris/zdistri_WFIRST_LSST_clustering_fine_bin_norm"

data_vector_file = "datav/DESI2_KL_v2_%d_ss_Ncl%d_Ntomo%d_LCDM+Eagle"
invcovmat_file = "invcov/DESI2_KL_v2_%d%d_ssss_invcov_Ncl%d_Ntomo%d"
baryon_PCS_file = "datav/DESI2_KL_v2_%d_ss_Ncl%d_Ntomo%d.pca"
chain_output_file = "chains/DESI2_KL_v2_%d%d_ss_Ncl%d_Ntomo%d"

############################################################
file_source_z = os.path.join(dirname, nz_src_files[args.iSelection])
file_lens_z = os.path.join(dirname, nz_lens_file)
data_file = os.path.join(dirname, 
    data_vector_file%(args.iSelection, Ncl, Ntomo_src))
cov_file = os.path.join(dirname, 
    invcovmat_file%(args.iSelection, args.iSN, Ncl, Ntomo_src))
bary_file = os.path.join(dirname, 
    baryon_PCS_file%(args.iSelection, Ncl, Ntomo_src))
chain_file = os.path.join(outdirname, 
    chain_output_file%(args.iSelection, args.iSN, Ncl, Ntomo_src))

initcosmo("halofit")
initbins(Ncl,ell_min,ell_max,ell_max_shear,Rmin_bias,Ntomo_src,Ntomo_lens)
#initpriors_KL("photo_opti","shear_opti","none","none")
initpriors_IA_bary("spec_DESI2", "shear_KL_DESI2", "none", "none", 
    False, 3.0, 1.2, 3.8, 2.0, 
    True, 16.0, 5.0, 0.8)
initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initclusters()
initia("none","none")
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

sample_main(sample_params,7000,700,1,chain_file+"_5000", blind=False, pool=MPIPool(), KL=True)

