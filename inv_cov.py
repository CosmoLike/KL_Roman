#!/usr/bin/python
import sys, os
import math, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import linalg as LA
import numpy as np


# covariance matrix
# DATA_DIR = '/xdisk/timeifler/jiachuanxu/DESI2KL/'
# #infile_fmt = "DESI2_KL_v2_%d%d_ssss_cov_Ncl%d_Ntomo%d"
# #outfile_fmt = "DESI2_KL_v2_%d%d_ssss_invcov_Ncl%d_Ntomo%d"
# infile_fmt = "LSST_Y%d_ssss_cov_Ncl15_Ntomo10"
# outfile_fmt = "LSST_Y%d_ssss_invcov_Ncl15_Ntomo10"
DATA_DIR = '../../dsa/'
'''
### DSA-2000
infile_fmt = "DSA_allsky_ssss_cov_Ncl15_Ntomo1_OneComp"
outfile_fmt = "DSA_allsky_ssss_invcov_Ncl15_Ntomo1_OneComp"
Ncl = 15
Area_list = [30000]
N_area = 1
N_selection = 1
N_tomo_list = [1]
Nsrc_list = np.array([0.0459])*3600
N_shape_noise = 1
SN_list = [0.05]
'''
### SKA
infile_fmt = "SKA_WL_ssss_cov_Ncl15_Ntomo1"
outfile_fmt = "SKA_WL_ssss_invcov_Ncl15_Ntomo1"
# infile_fmt = "SKA_KL_ssss_cov_Ncl15_Ntomo1_OneComp"
# outfile_fmt = "SKA_KL_ssss_invcov_Ncl15_Ntomo1_OneComp"
Ncl = 15
Area_list = [30000]
N_area = 1
N_selection = 1
N_tomo_list = [1]
Nsrc_list = np.array([4.49])*3600
N_shape_noise = 1
SN_list = [0.3]		# shape noise: 0.3 (WL) or 0.05 (KL)
plot_corrmat = True

for iArea in range(N_area):
	#Area = Area_list[iArea]
	for jSelect in range(N_selection):
		N_tomo = N_tomo_list[jSelect]
		Nsrc = Nsrc_list[jSelect]
		#year = years[jSelect]
		Area = Area_list[jSelect]
		for kSN in range(N_shape_noise):
			SN = SN_list[kSN]
			fig_title = "$\Omega_s=$%d deg$^{2}$, $n_\mathrm{src}=$%d deg$^{-2}$, $\sigma_\epsilon^\mathrm{rms}$=%.2f"%(Area, Nsrc, SN)
			#infile = DATA_DIR+"cov/"+infile_fmt%(jSelect,kSN,Ncl,N_tomo)
			infile = DATA_DIR+"cov/"+infile_fmt#%(year)
			if not os.path.exists(infile):
				print(f'File {infile} does not exist! Continue')
				continue
			#outname = DATA_DIR+"invcov/"+outfile_fmt%(jSelect,kSN,Ncl,N_tomo)
			outname = DATA_DIR+"invcov/"+outfile_fmt#%(year)
			### Set data vector dimensions
			nggl = 0 	# number of ggl power spectra
			ngcl = 0	# number of cluster-source galaxy power spectra
			nlens = 0 	# number of lens bins 
			nlenscl= 0 	# number of cluster redshift bins 
			nshear = int((N_tomo+1)*N_tomo/2) # # of shear tomo power spectra
			ncl=int(Ncl)		# number of ell-bins
			nclgcl=0	# number of cluster ell-bins
			nrich=0		# number of richness bins
			ndata = int((nshear+nggl+nlens)*ncl+nlenscl*nrich+nrich*ngcl*nclgcl)
			n2pt = int((nshear+nggl+nlens)*ncl)
			ncluster = int(nlenscl*nrich)
			n2ptcl=int(n2pt+ncluster)
			nclusterN_WL=int(ncluster+nrich*ngcl*nclgcl)
			print(f'Data Vector Size = {ndata}')
			### Just use all-1 mask for new
			mask = np.ones(ndata)
			### Read covariance matrix
			covfile = np.genfromtxt(infile)
			_ndata = int(np.max(covfile[:,0])+1)
			assert _ndata == ndata, f'Inconsistent Ndata {ndata} v.s. {_ndata}'
			cov = np.zeros((ndata,ndata))
			for l in range(0,covfile.shape[0]):
				bin1, bin2 = int(covfile[l,0]), int(covfile[l,1])
				cov_g, cov_ng = covfile[l,8], covfile[l,9]
				cov[bin1, bin2] = cov_g+cov_ng
				cov[bin2, bin1] = cov_g+cov_ng
			### And calculate correlation matrix
			cor = cov/np.outer(np.diagonal(cov)**0.5, np.diagonal(cov)**0.5)
			# for i in range(0,ndata):
			#     for j in range(0,ndata):
			#     	if (cov[i,i]*cov[j,j] >0):
			#        		cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])
			a = np.sort(LA.eigvals(cor[:,:]))
			print("Eigenvalues range of corrmat (%.3e, %.3e)"%(np.min(a), np.max(a)))
			print("Negative eigenvalues of corrmat:{}".format(a[a<0]))
			# ############### invert shear covariance #################
			inv = LA.inv(cov[0:nshear*ncl,0:nshear*ncl])
			a = np.sort(LA.eigvals(cov[0:nshear*ncl,0:nshear*ncl]))
			print("Eigvals range of covmat (%.3e, %.3e)"%(np.min(a), np.max(a)))
			f = open(outname, 'w')
			for i in range(0,nshear*ncl):
				for j in range(0,nshear*ncl):
					f.write("%d %d %e\n" %(i, j, inv[i,j]*mask[i]*mask[j]))
			f.close()
			'''	
			# ############### invert clustering covariance #################
			inv = LA.inv(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl])
			a = np.sort(LA.eigvals(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl]))
			print "min+max eigenvalues clustering cov:"
			print np.min(a), np.max(a)
			outfile = "cov/"+outname[k]+"_pos_pos_inv"
			f = open(outfile, "w")
			for i in range(0,nlens*ncl):
				inv[i,i]=inv[i,i]*mask[(nshear+nggl)*ncl+i]
				for j in range(0,nlens*ncl):
					f.write("%d %d %e\n" %(i,j, inv[i,j]))
			f.close()

			# ############### invert 2pt covariance #################
			a = np.sort(LA.eigvals(cov[0:n2pt,0:n2pt]))
			print "min+max eigenvalues 2pt cov:"
			print np.min(a), np.max(a)
			inv = LA.inv(cov[0:n2pt,0:n2pt])
			outfile = "cov/"+outname[k]+"_3x2pt_inv" 
			f = open(outfile, "w")
			for i in range(0,n2pt):
				inv[i,i]=inv[i,i]*mask[i]
				for j in range(0,n2pt):
					f.write("%d %d %e\n" %( i,j, inv[i,j]))
			f.close()



			# # ############### invert full2pt+clusterN+clusterWL covariance #################
			precond = 1.e-7
			for i in range(0,ncluster):
			  cov[n2pt+i,:]*= precond
			  cov[:,n2pt+i]*= precond
			inv = LA.inv(cov)
			a = np.sort(LA.eigvals(cov))
			print "min+max eigenvalues of full 2ptclusterN+clusterWL pre-conditioned matrix:"
			print np.min(a), np.max(a)
			if (np.min(a)<0):
			  print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile[k])
			for i in range(0,ncluster):
			  inv[n2pt+i,:]*= precond
			  inv[:,n2pt+i]*= precond

			outfile = "cov/"+outname[k]+"_3x2pt_clusterN_clusterWL_inv"
			f = open(outfile, "w")
			for i in range(0,ndata):
			  inv[i,i]=inv[i,i]*mask[i]
			  for j in range(0,ndata):
			    f.write("%d %d %e\n" %( i,j, inv[i,j]))
			f.close()



			# # ############### invert clusterN+clusterWL covariance #################
			inv = LA.inv(cov[n2pt:n2pt+nclusterN_WL,n2pt:n2pt+nclusterN_WL])
			a = np.sort(LA.eigvals(cov[n2pt:n2pt+nclusterN_WL,n2pt:n2pt+nclusterN_WL]))
			print "min+max eigenvalues of clusterN_WL pre-conditioned matrix:"
			print np.min(a), np.max(a)
			if (np.min(a)<0):
			  print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile[k])

			for i in range(0,ncluster):
			  inv[i,:]*= precond
			  inv[:,i]*= precond

			outfile = "cov/"+outname[k]+"_clusterN_clusterWL_inv"
			f = open(outfile, "w")
			for i in range(0,nclusterN_WL):
			  	inv[i,i]=inv[i,i]*mask[n2pt+i]
			  	for j in range(0,nclusterN_WL):
			  		f.write("%d %d %e\n" %( i,j, inv[i,j]))
			f.close()
			'''
			### Plot correlation matrix or covariance matrix
			if plot_corrmat:
				labels = [
	r'$C^{\kappa \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
	# r'$C^{g \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
	# r'$C^{gg}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
	# r'$N^{c}\left(r,z_{\mathrm{c}_i}\right)$',
	# r'$C^{c \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
				]
				ticks = np.zeros(2)
				ticks[1] = nshear*ncl
				# ticks[2] = (nshear+nggl)*ncl
				# ticks[3] = n2pt
				# ticks[4] = n2pt+ncluster
				# ticks[5] = ndata
				tickx = 0.5*(ticks[:-1]+ticks[1:])

				cor = cov/np.outer(np.diagonal(cov)**0.5, np.diagonal(cov)**0.5)
				cor = cor*np.outer(mask, mask)
				# cor = np.zeros((ndata,ndata))
				# for i in range(0,ndata):
				# 	for j in range(0,ndata):
				# 		if (cov[i,i]*cov[j,j] >0):
				# 			cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])

				fs= 10
				# for i in range(0,5):
			  	# 	tickx[i] = 0.5*(ticks[i]+ticks[i+1])
			  	# 	plt.plot([ticks[i]-0.5,ticks[i]-0.5],[-.5,ndata-0.5],linestyle ='-',color = 'k')
			  	# 	plt.plot([-.5,ndata-0.5],[ticks[i]-0.5,ticks[i]-0.5],linestyle ='-',color = 'k')

				plt.subplot(1, 1, 1)
				ax = plt.gca()
				# im = ax.imshow(np.log10(cov[:,:]), interpolation="nearest",vmin=-25, vmax=-10)
				im = ax.imshow(cor, interpolation='nearest', origin='lower', vmin=-1, vmax=1,
								cmap='seismic')
				ax.set_title(fig_title, fontsize=16)
				plt.xticks(tickx, labels,fontsize=fs)
				plt.yticks(tickx-0.5, labels,fontsize=fs)
				#plt.tick_params(axis = 'x',length = 0, pad = 15)
				#plt.tick_params(axis = 'y',length = 0, pad = 5)

				plt.colorbar(im)
				#plt.show()
				figfilename_fmt = DATA_DIR+"figure/"+outfile_fmt+".png"
				#plt.savefig(figfilename_fmt%(jSelect,kSN,Ncl,N_tomo), format="png")
				# plt.savefig(figfilename_fmt%(year), format="png")
				plt.savefig(figfilename_fmt, format="png")
