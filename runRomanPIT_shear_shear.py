import sys
sys.path.append('/home/u17/jiachuanxu/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('i_depth', type=int, help='Index of survey depth scenario')
parser.add_argument('i_ellmax', type=int, help='Index of ell max')
parser.add_argument('baryPCA', type=int, help='Marg. over PC1')
parser.add_argument('-nsteps', type=int, default=2000, help='MCMC steps')
parser.add_argument('-nwalkers', type=int, default=400, help='N walkers')
args = parser.parse_args()

ellmax_list = [1000., 2000., 3000., 4000.]
neff_list = [20., 25., 30., 35, 40.]
Q1_priors = { # use Illustris Q1 value for the std. of the Q1 prior
    '00' : 7,
    '01' : 23,
    '02' : 42,
    '03' : 69,
    '10' : 8,
    '11' : 26,
    '12' : 48,
    '13' : 84,
    '20' : 9,
    '21' : 28,
    '22' : 54,
    '23' : 98,
    '30' : 9,
    '31' : 30,
    '32' : 59,
    '33' : 111,
    '40' : 10,
    '41' : 32,
    '42' : 64,
    '43' : 124,
}
ellmax = ellmax_list[args.i_ellmax]
neff = neff_list[args.i_depth]
Q1_std = Q1_priors["%d%d"%(args.i_depth, args.i_ellmax)]
#print "ell max = %.0f; neff = %.0f; Q1 std = %.1f"%(ellmax, neff, Q1_std)

#########################################################
dirname = "/home/u17/jiachuanxu/CosmoLike/KL_WFIRST"
outdirname = "/xdisk/timeifler/jiachuanxu/RomanPIT"
Ntomo_src, Ntomo_lens = 10, 10
Ncl = 15
ell_min, ell_max_shear = 20.0, 3000.0
Rmin_bias = 21.0
strat = "Roman_WL_%d%d"%(args.i_depth, args.i_ellmax)
#print "Survey strat = ", strat
nz_src_file = "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin_norm_neff%.0f"%(neff)
nz_lens_file = "zdistris/lens_LSSTY1"
data_vector_file = "datav/Roman_WL_%d%d_shear_shear_Ntomo10_Ncl15_dmo"%(args.i_depth, args.i_ellmax)
invcovmat_file = "invcov/Roman_WL_%d%d_ssss_invcov_Ncl15_Ntomo10"%(args.i_depth, args.i_ellmax)
baryon_PCS_file = "datav/Roman_WL_%d%d_shear_shear_Ntomo10_Ncl15_9sim.pca"%(args.i_depth, args.i_ellmax)
#chain_output_file = "chains/LSST_Y1_ss_Ncl%d_Ntomo%d"
chain_output_file = "chains/Roman_WL_%d%d_zlow1_ss_Ncl15_Ntomo10"%(args.i_depth, args.i_ellmax)
external_probe = "none"
NPCs_used = args.baryPCA
if NPCs_used>0:
    samp_bary = True
else:
    samp_bary = False
#cosmo_model = "LCDM_split"
cosmo_model = "s8split_only"
runmode = "halofit_split"
############################################################
file_source_z = os.path.join(dirname, nz_src_file)
file_lens_z = os.path.join(dirname, nz_lens_file)
data_file = os.path.join(dirname, data_vector_file)
cov_file = os.path.join(outdirname, invcovmat_file)
bary_file = os.path.join(dirname, baryon_PCS_file)
chain_file = os.path.join(outdirname, chain_output_file)

initcosmo(runmode)
initbins(Ncl,ell_min,ell_max,ell_max_shear,Rmin_bias,Ntomo_src,Ntomo_lens)
initpriors_IA_bary("photo_opti","shear_opti","none",external_probe,
    False, 3.0,1.2,3.8,2.0,
    samp_bary, Q1_std,10.0,0.8)
initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
#initclusters()
initia("NLA_HF","GAMA")
initprobes("shear_shear")
initdatainvbary(cov_file, data_file, bary_file)

sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear(), 
    MG=False, NPCs=NPCs_used, cosmology=cosmo_model, source_photo_z=False, 
    shear_calibration=False, IA=False)

sample_main(sample_params,args.nsteps,args.nwalkers,1,chain_file+"_%d"%args.nsteps, blind=False, pool=MPIPool())
