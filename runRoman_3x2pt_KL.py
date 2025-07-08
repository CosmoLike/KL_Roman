import sys, os
sys.path.append('/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST')

from cosmolike_libs_opti import *
from schwimmbad import MPIPool
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-nsteps", type=int, default=2000, help="MCMC steps")
parser.add_argument("-nwalkers", type=int, default=400, help="N walkers")
args = parser.parse_args()
print "MCMC sampler steps = %d, walkers = %d"%(args.nsteps, args.nwalkers)

def parse_ini_file(fname):
    params = {}
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue
            if ":" in line:
                key, value = map(str.strip, line.split(":", 1))
                try:
                    params[key] = int(value)
                except ValueError:
                    try:
                        params[key] = float(value)
                    except ValueError:
                        params[key] = value
    return params

# path
dirname = "/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST/"
outdirname = "/home/u15/yhhuang/cosmology/CosmoLike/3Dx2D/"

# file fmt
data_vector_file = "datav/Roman_KL_3x2pt_Ntomo%d_Ncl%d_dmo"
invcovmat_file = "invcov/Roman_KL_3x2pt_invcov_Ncl%d_Ntomo%d"
chain_output_file = "chains/Roman_KL_3x2pt_Ncl%d_Ntomo%d"

# external prior
external_prior = "none"

# flag
DE_FLAG = False
KL_FLAG = True

# read parameters from ini file
params = parse_ini_file("params.ini")
Ntomo_src = params.get("Ntomo_source")
Ntomo_lens = params.get("Ntomo_lens")
Ncl = params.get("Ncl")
ell_min = params.get("lmin")
ell_max = params.get("lmax")
ell_max_shear = params.get("lmax_shear")
Rmin_bias = params.get("Rmin_bias")
nz_src_files = params.get("shear_REDSHIFT_FILE")
nz_lens_files = params.get("clustering_REDSHIFT_FILE")

###########################################################
file_source_z = os.path.join(dirname, nz_src_files)
file_lens_z = os.path.join(dirname, nz_lens_files)
data_file = os.path.join(dirname, data_vector_file%(Ntomo_src, Ncl))
cov_file = os.path.join(outdirname, invcovmat_file%(Ncl, Ntomo_src))
chain_file = os.path.join(outdirname, chain_output_file%(Ncl, Ntomo_src))

initcosmo("halofit")
initbins(Ncl, ell_min, ell_max, ell_max_shear, Rmin_bias, Ntomo_src, Ntomo_lens)
initpriors_IA_bary(
    "spec_Roman", "shear_Roman", "none", "none",
    False, 3.0, 1.2, 3.8, 2.0,
    False, 16, 1.9, 0.7
)
initsurvey("Roman_KL")
initgalaxies(file_source_z, file_lens_z, "gaussian", "gaussian", "SN10")
initclusters()
initia("none", "GAMA")
initprobes("3x2pt")
initdatainv(cov_file, data_file)

# sample_params = sample_cosmology_only()
sample_params = sample_cosmology_clustering_KL(Ntomo_lens)

sample_main(sample_params, args.nsteps, args.nwalkers, 1, chain_file+'_%d'%args.nsteps, 
            blind=False, pool=MPIPool(), KL=KL_FLAG)