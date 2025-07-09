#!/usr/bin/python
import sys, os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import LogNorm
from numpy import linalg as LA

def parse_ini_file(fname):
    params = {}
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if ':' in line:
                key, value = map(str.strip, line.split(':', 1))
                try:
                    params[key] = int(value)
                except ValueError:
                    try:
                        params[key] = float(value)
                    except ValueError:
                        params[key] = value
    return params

if len(sys.argv) != 2:
    print('Use: python [ini_file]')
    sys.exit()

params = parse_ini_file(sys.argv[1])
Ncl = params.get('Ncl')
Ntomo = params.get('Ntomo_source')
Area = params.get('survey_area')
Nsrc = params.get('n_gal')
ShapeNoise = 0.05
ellmax = params.get('lmax')
probes = params.get('probes')


infile_fmt = 'cov/Roman_KL_%s_cov_Ncl%d_Ntomo%d'
outfile_fmt = 'invcov/Roman_KL_%s_invcov_Ncl%d_Ntomo%d'
DATA_DIR = '../3Dx2D/'

# the numbers below can be computed knowing the data vector settings, 
# e.g. 10 tomographic source bins results in 55 shear-shear power spectra. 
# Or they can be read off when running the covariance code, i.e. type 'compute_covariance_fourier 100000' 
# and look for the output mentioning number of ggl bins accepted and/or 
# number of cluster weak lensing bins accepted. 
# The default numbers below most likely don't correspond to your binning choices.
nggl    = 21                        # number of ggl power spectra
nlens   = int(Ntomo)                # number of lens bins 
nshear  = int((Ntomo+1)*Ntomo/2)    # number of shear tomographic power spectra
nlenscl = 0                         # number of cluster redshift bins 
ncl     = int(Ncl)                  # number of ell-bins
nclgcl  = 0                         # number of cluster ell-bins
ngcl    = 0                         # number of cluster-source galaxy power spectra
nrich   = 0                         # number of richness bins

# set data vector dimensions
nbuff = 0
if probes == '3x2pt':
    probes_list = ['ssss', 'llll', 'lsls', 'llls', 'lsss', 'llss']
elif probes == 'ggl_cl':
    probes_list = ['llll', 'lsls', 'llls']
    nbuff = int(nshear * ncl)
    nshear = 0
elif probes == 'shear_shear':
    probes_list = ['ssss']
    nlens = 0
    nggl = 0
else:
    print(f'Unsupported probes: {probes}')
    exit(-1)

ndata   = int((nshear+nggl+nlens)*ncl + nlenscl*nrich + nrich*ngcl*nclgcl)
n2pt    = int((nshear+nggl+nlens)*ncl)
ncluster = int(nlenscl*nrich)
n2ptcl = int(n2pt+ncluster)
nclusterN_WL = int(ncluster+nrich*ngcl*nclgcl)
print(f'Data Vector Size = {ndata}')

# use all-1 mask for new
mask = np.ones(ndata)

# figure format
plot_corrmat = True # True for correlation matrix, False for covariance matrix
plot_separate = False   # True for separate plots, False for combined plot
fontsize = 10

fig_title = "$n_\mathrm{src}=$%.2f arcmin$^{-2}$, $\ell_\mathrm{max}$=%.0f"%(Nsrc, ellmax)
if plot_corrmat:
    fig_filename = "../3Dx2D/figure/Roman_KL_%s_corr_Ncl%d_Ntomo%d.png"%(probes, Ncl, Ntomo)
else:
    fig_filename = "../3Dx2D/figure/Roman_KL_%s_cov_Ncl%d_Ntomo%d.png"%(probes, Ncl, Ntomo)

if plot_separate and probes == '3x2pt':
    fig_list = [
        "../3Dx2D/figure/Roman_KL_shear_cov_Ncl%d_Ntomo%d.png"%(Ncl, Ntomo),
        "../3Dx2D/figure/Roman_KL_ggl_cov_Ncl%d_Ntomo%d.png"%(Ncl, Ntomo),
        "../3Dx2D/figure/Roman_KL_lens_cov_Ncl%d_Ntomo%d.png"%(Ncl, Ntomo)
    ]
    start_list = [0, nshear*ncl, (nshear+nggl)*ncl]
    end_list = [nshear*ncl, (nshear+nggl)*ncl, n2pt]

cov = np.zeros((ndata, ndata))
for probe in probes_list:
    infile = DATA_DIR + infile_fmt%(probe, Ncl, Ntomo)
    if not os.path.exists(infile):
        print(f'File {infile} does not exist! Continue')
        continue
    
    covfile = np.genfromtxt(infile)
    for l in range(0, covfile.shape[0]):
        bin1, bin2 = int(covfile[l,0]), int(covfile[l,1])
        # normal format: g/ng
        if covfile.shape[1] == 10:
            cov_g, cov_ng = covfile[l,8], covfile[l,9]
        # break-down version
        elif covfile.shape[1] == 12:
            cov_g = covfile[l,8] + covfile[l,9] + covfile[l,10]
            cov_ng = covfile[l,11]
        # covariance only
        elif covfile.shape[1] == 3:
            cov_g = 0; cov_ng = covfile[l,2]
        else:
            print(f'Wrong covariance matrix file format! {covfile.shape[1]} columns')
            exit(-1)
        
        if nbuff != 0:
            bin1 = bin1 - nbuff
            bin2 = bin2 - nbuff
        if bin1 == bin2:
            factor = 1
        else:
            factor = mask[bin1] * mask[bin2]
        
        cov[bin1, bin2] += cov_g * factor
        cov[bin2, bin1] += cov_g * factor
# correlation matrix
cor = cov / np.outer(np.diagonal(cov)**0.5, np.diagonal(cov)**0.5)
a = np.sort(LA.eigvals(cor[:,:]))
print("Eigenvalues range of corrmat (%.3e, %.3e)"%(np.min(a), np.max(a)))
print("Negative eigenvalues of corrmat:{}".format(a[a<0]))

############### invert 2pt covariance #################
a = np.sort(LA.eigvals(cov[0:n2pt,0:n2pt]))
print("Eigenvalues range of 2pt cov: (%.3e, %.3e)"%(np.min(a), np.max(a)))

inv = LA.inv(cov[0:n2pt,0:n2pt])
outfile = DATA_DIR + outfile_fmt%(probes, Ncl, Ntomo)
f = open(outfile, "w")
for i in range(0,n2pt):
    inv[i,i]=inv[i,i]*mask[i]
    for j in range(0,n2pt):
        f.write("%d %d %e\n" %( i,j, inv[i,j]))
f.close()

# ############### invert shear covariance #################
# inv = LA.inv(cov[0:nshear*ncl,0:nshear*ncl])
# a = np.sort(LA.eigvals(cov[0:nshear*ncl,0:nshear*ncl]))
# print("Eigvals range of covmat (%.3e, %.3e)"%(np.min(a), np.max(a)))
# f = open(outname, 'w')
# for i in range(0,nshear*ncl):
#   for j in range(0,nshear*ncl):
#       f.write("%d %d %e\n" %(i, j, inv[i,j]*mask[i]*mask[j]))
# f.close()

# ############### invert clustering covariance #################
# inv = LA.inv(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl])
# a = np.sort(LA.eigvals(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl]))
# print("min+max eigenvalues clustering cov:")
# print(np.min(a), np.max(a))
# outfile = "cov/"+outname[k]+"_pos_pos_inv"
# f = open(outfile, "w")
# for i in range(0,nlens*ncl):
#   inv[i,i]=inv[i,i]*mask[(nshear+nggl)*ncl+i]
#   for j in range(0,nlens*ncl):
#       f.write("%d %d %e\n" %(i,j, inv[i,j]))
# f.close()

# ############### invert full2pt+clusterN+clusterWL covariance #################
# precond = 1.e-7
# for i in range(0,ncluster):
#   cov[n2pt+i,:]*= precond
#   cov[:,n2pt+i]*= precond
# inv = LA.inv(cov)
# a = np.sort(LA.eigvals(cov))
# print "min+max eigenvalues of full 2ptclusterN+clusterWL pre-conditioned matrix:"
# print np.min(a), np.max(a)
# if (np.min(a)<0):
#   print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile[k])
# for i in range(0,ncluster):
#   inv[n2pt+i,:]*= precond
#   inv[:,n2pt+i]*= precond

# outfile = "cov/"+outname[k]+"_3x2pt_clusterN_clusterWL_inv"
# f = open(outfile, "w")
# for i in range(0,ndata):
#   inv[i,i]=inv[i,i]*mask[i]
#   for j in range(0,ndata):
#   f.write("%d %d %e\n" %( i,j, inv[i,j]))
# f.close()



# ############### invert clusterN+clusterWL covariance #################
# inv = LA.inv(cov[n2pt:n2pt+nclusterN_WL,n2pt:n2pt+nclusterN_WL])
# a = np.sort(LA.eigvals(cov[n2pt:n2pt+nclusterN_WL,n2pt:n2pt+nclusterN_WL]))
# print "min+max eigenvalues of clusterN_WL pre-conditioned matrix:"
# print np.min(a), np.max(a)
# if (np.min(a)<0):
#   print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile[k])

# for i in range(0,ncluster):
#   inv[i,:]*= precond
#   inv[:,i]*= precond

# outfile = "cov/"+outname[k]+"_clusterN_clusterWL_inv"
# f = open(outfile, "w")
# for i in range(0,nclusterN_WL):
#   inv[i,i]=inv[i,i]*mask[n2pt+i]
#   for j in range(0,nclusterN_WL):
#       f.write("%d %d %e\n" %( i,j, inv[i,j]))
# f.close()

### Plot correlation matrix or covariance matrix

labels = [
r'$C^{\kappa \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
r'$C^{g \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
r'$C^{gg}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
# r'$N^{c}\left(r,z_{\mathrm{c}_i}\right)$',
# r'$C^{c \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
]
ticks = np.zeros(4)
ticks[1] = nshear*ncl
ticks[2] = (nshear+nggl)*ncl
ticks[3] = n2pt
# ticks[4] = n2pt+ncluster
# ticks[5] = ndata
tickx = 0.5*(ticks[:-1]+ticks[1:])

fig = plt.figure(figsize=(10,10))
ax = plt.gca()
if plot_corrmat:
    cor = cov/np.outer(np.diagonal(cov)**0.5, np.diagonal(cov)**0.5)
    cor = cor*np.outer(mask, mask)

    im = ax.imshow(cor, interpolation='nearest', origin='lower', vmin=-1, vmax=1,
                    cmap='seismic')
else:
    im = ax.imshow(cov, interpolation='nearest', origin='lower', cmap='seismic',
                    norm=LogNorm(vmin=1e-25, vmax=2e-10))

ax.set_title(fig_title, fontsize=16)
plt.xticks(tickx, labels, fontsize=fontsize)
plt.yticks(tickx-0.5, labels, fontsize=fontsize)
plt.colorbar(im)
plt.savefig(fig_filename, format="png", dpi=150, bbox_inches='tight')
plt.close()

if plot_separate:
    print('Plot separate covariance matrices')
    for i, (fname, start, end) in enumerate(zip(fig_list, start_list, end_list)):
        fig = plt.figure(figsize=(10,10))
        ax = plt.gca()
        im = ax.imshow(cov[start:end, start:end], interpolation='nearest', origin='lower',
                        cmap='seismic', norm=LogNorm(vmin=1e-25, vmax=2e-10))
        ax.set_title(fig_title, fontsize=16)
        plt.xlabel(labels[i], fontsize=fontsize)
        plt.ylabel(labels[i], fontsize=fontsize)
        plt.colorbar(im)
        plt.savefig(fname, format="png", dpi=150, bbox_inches='tight')
        plt.close()






