import sys
sys.path.append('/home/u17/jiachuanxu/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('i_depth', type=int, help='Index of survey depth scenario')
parser.add_argument('i_area', type=int, help='Index of survey area scenario')
parser.add_argument('i_ellmax', type=int, help='Index of ell max')
parser.add_argument('baryPCA', type=int, help='Marg. over PC1')
parser.add_argument('-nsteps', type=int, default=2000, help='MCMC steps')
parser.add_argument('-nwalkers', type=int, default=400, help='N walkers')
args = parser.parse_args()
assert (args.i_depth>=0) and (args.i_depth<5)
assert (args.i_area>=0) and (args.i_area<5)
assert (args.i_ellmax>=0) and (args.i_ellmax<4)
assert (args.baryPCA>=0) and (args.baryPCA<2)
print "MCMC sampler steps = %d, walkers = %d"%(args.nsteps, args.nwalkers)
ellmax_list = [1000., 2000., 3000., 4000.]
neff_list = [20., 25., 30., 35, 40.]
area_list = [1000.0, 4000.0, 7000.0, 10000.0, 13000.0]
# Q1_priors = { # use Illustris Q1 value for the std. of the Q1 prior
#     '00' : 3.9,
#     '01' : 7.7,
#     '02' : 10.0,
#     '03' : 12.0,
#     '10' : 4.5,
#     '11' : 9.3,
#     '12' : 12.3,
#     '13' : 15.2,
#     '20' : 5.1,
#     '21' : 10.8,
#     '22' : 14.6,
#     '23' : 18.6,
#     '30' : 5.7,
#     '31' : 12.2,
#     '32' : 16.9,
#     '33' : 22.3,
#     '40' : 6.2,
#     '41' : 13.7,
#     '42' : 19.2,
#     '43' : 26.1,
# }
ell_max = ellmax_list[args.i_ellmax]
neff = neff_list[args.i_depth]
# Q1_std = Q1_priors["%d%d"%(args.i_depth, args.i_ellmax)]
area = area_list[args.i_area]
#print "ell max = %.0f; neff = %.0f; Q1 std = %.1f"%(ell_max, neff, Q1_std)
print "ell max = %.0f; neff = %.0f; area = %.0f; bary = %d"%(ell_max, neff, area, args.baryPCA)

#########################################################
dirname = "/home/u17/jiachuanxu/CosmoLike/KL_WFIRST"
outdirname = "/xdisk/timeifler/jiachuanxu/RomanPIT"
Ntomo_src, Ntomo_lens = 10, 10
Ncl = 15
ell_min, ell_max_shear = 20.0, ellmax_list[args.i_ellmax]
Rmin_bias = 21.0
strat = "Roman_WL_%d%d%d"%(args.i_depth, args.i_area, args.i_ellmax)
#print "Survey strat = ", strat
#nz_src_file = "zdistri_WFIRST_LSST_lensing_fine_bin_norm_neff%.0f"%(neff)
nz_src_file = "n_eff%.0f.nz"%(neff)
nz_lens_file = "lens_LSSTY1"
# data_vector_file = "Roman_WL_%d%d_shear_shear_Ntomo10_Ncl15_dmo"%(args.i_depth, args.i_ellmax)
data_vector_file = "Roman_WL_%d0%d_shear_shear_Ntomo10_Ncl15_Haley_dndz_dmo"%(args.i_depth, args.i_ellmax)
invcovmat_file = "Roman_WL_%d%d%d_ssss_invcov_Ncl15_Ntomo10_Haley_dndz_v2"%(args.i_depth, args.i_area, args.i_ellmax)
baryon_PCS_file = "Roman_WL_%d%d_shear_shear_Ntomo10_Ncl15_9sim.pca"%(args.i_depth, args.i_ellmax)
#chain_output_file = "chains/LSST_Y1_ss_Ncl%d_Ntomo%d"
external_probe = "better_DESI_BAO"
#NPCs_used = args.baryPCA
if args.baryPCA==1:
    samp_bary = True
    # use PC1 marginalized inverse covariance
    # invcovmat_file = invcovmat_file.replace("invcov", "invcov_barymarg_2std")
    invcovmat_file = invcovmat_file.replace("invcov_Ncl15", 
        "invcov_barymarg_illustris_Ncl15")
    chain_output_file = "Roman_WL_%d%d%d_ss_Ncl15_Ntomo10_barymarg_illustris_Haley_dndz_v2"%(args.i_depth, args.i_area, args.i_ellmax)
    print "Baryonic effects marginalized! (illustris)"
elif args.baryPCA==0:
    samp_bary = False
    chain_output_file = "Roman_WL_%d%d%d_ss_Ncl15_Ntomo10_fixbary_Haley_dndz_v2"%(args.i_depth, args.i_area, args.i_ellmax)
    print "No baryonic effects marginalized!"
# elif args.baryPCA==2:
#     samp_bary = True
#     # use PC1 marginalized inverse covariance
#     invcovmat_file = invcovmat_file.replace("invcov", "invcov_barymarg_5std")
#     chain_output_file = "Roman_WL_%d%d%d_zlow1_ss_Ncl15_Ntomo10_margbary_5std_s8sl"%(args.i_depth, args.i_area, args.i_ellmax)
#     print "Baryonic effects marginalized! (5xstd)"

#cosmo_model = "LCDM_split"
#cosmo_model = "s8split_only"
cosmo_model = "OmS8"
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
initpriors_IA_bary("photo_opti","shear_opti","none",external_probe,
    False, 3.0,1.2,3.8,2.0,
    False, 12,10.0,0.8)
initsurvey(strat)
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
#initclusters()
initprobes("shear_shear")
initdatainvbary(cov_file, data_file, bary_file)

sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear(), 
    MG=False, NPCs=0, cosmology=cosmo_model, source_photo_z=False, 
    shear_calibration=False, IA=False)
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
