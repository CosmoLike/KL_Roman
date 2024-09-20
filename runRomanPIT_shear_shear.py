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
    '00' : 3.9,
    '01' : 7.7,
    '02' : 10.0,
    '03' : 12.0,
    '10' : 4.5,
    '11' : 9.3,
    '12' : 12.3,
    '13' : 15.2,
    '20' : 5.1,
    '21' : 10.8,
    '22' : 14.6,
    '23' : 18.6,
    '30' : 5.7,
    '31' : 12.2,
    '32' : 16.9,
    '33' : 22.3,
    '40' : 6.2,
    '41' : 13.7,
    '42' : 19.2,
    '43' : 26.1,
}
ell_max = ellmax_list[args.i_ellmax]
neff = neff_list[args.i_depth]
Q1_std = Q1_priors["%d%d"%(args.i_depth, args.i_ellmax)]
print "ell max = %.0f; neff = %.0f; Q1 std = %.1f"%(ell_max, neff, Q1_std)

#########################################################
dirname = "/home/u17/jiachuanxu/CosmoLike/KL_WFIRST"
outdirname = "/xdisk/timeifler/jiachuanxu/RomanPIT"
Ntomo_src, Ntomo_lens = 10, 10
Ncl = 15
ell_min, ell_max_shear = 20.0, ellmax_list[args.i_ellmax]
Rmin_bias = 21.0
strat = "Roman_WL_%d%d"%(args.i_depth, args.i_ellmax)
#print "Survey strat = ", strat
nz_src_file = "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin_norm_neff%.0f"%(neff)
nz_lens_file = "zdistris/lens_LSSTY1"
data_vector_file = "datav/Roman_WL_%d%d_shear_shear_Ntomo10_Ncl15_dmo"%(args.i_depth, args.i_ellmax)
invcovmat_file = "invcov/Roman_WL_%d%d_ssss_invcov_Ncl15_Ntomo10"%(args.i_depth, args.i_ellmax)
baryon_PCS_file = "datav/Roman_WL_%d%d_shear_shear_Ntomo10_Ncl15_9sim.pca"%(args.i_depth, args.i_ellmax)
#chain_output_file = "chains/LSST_Y1_ss_Ncl%d_Ntomo%d"
external_probe = "none"
NPCs_used = args.baryPCA
if NPCs_used>0:
    samp_bary = True
else:
    samp_bary = False
print "NPCs used = %d", NPCs_used
chain_output_file = "chains/Roman_WL_%d%d_zlow1_ss_Ncl15_Ntomo10_PC%d"%(args.i_depth, args.i_ellmax, args.baryPCA)

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
initia("none","GAMA")
initpriors_IA_bary("photo_opti","shear_opti","none",external_probe,
    False, 3.0,1.2,3.8,2.0,
    samp_bary, Q1_std,10.0,0.8)
initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
#initclusters()
initprobes("shear_shear")
initdatainvbary(cov_file, data_file, bary_file)

sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear(), 
    MG=False, NPCs=NPCs_used, cosmology=cosmo_model, source_photo_z=False, 
    shear_calibration=False, IA=False)

### test likelihood evaluation bias
if samp_bary:
    parval = [0.851, 0.831, 0.2]
else:
    parval = [0.851, 0.831]
test_logpost = test_likelihood(sample_params, parval)
print "test likelihood:", test_logpost

### run mcmc chains
sample_main(sample_params,args.nsteps,args.nwalkers,1,chain_file+"_%d"%args.nsteps, blind=False, pool=MPIPool())
