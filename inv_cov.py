#!/usr/bin/python
import sys, os
import math, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import linalg as LA
import numpy as np

Ncl = 15
#which_survey = "DESI2" # "LSST"
which_survey = "RomanPIT"
#DATA_DIR = '/xdisk/timeifler/jiachuanxu/DESI2KL/'
DATA_DIR = "/xdisk/timeifler/jiachuanxu/RomanPIT/"

if which_survey=="DESI2":
	# covariance matrix for DESI-II KL
	infile_fmt = "DESI2_KL_v3_%d%d_ssss_cov_Ncl%d_Ntomo%d"
	outfile_fmt = "DESI2_KL_v3_%d%d_ssss_invcov_Ncl%d_Ntomo%d"
	Area_list = [14000, 14000, 14000, 14000, 14000, 14000];
	N_tomo_list = [4, 4, 4, 4, 4, 4];N_selection = 6
	#Nsrc_list = np.array([0.4761, 0.1629, 0.1553, 0.0881, 0.1006, 0.0740])*3600 # old TS
	Nsrc_list = np.array([0.5256, 0.1717, 0.1789, 0.0953, 0.1189, 0.0811])*3600 # new TS
	#SN_list = [0.02*1.4142, 0.04*1.4142, 0.06*1.4142, 0.10*1.4142, 0.20*1.4142, 0.30*1.4142]
	SN_list = [0.04, 0.05, 0.06, 0.07, 0.09, 0.11];N_shape_noise = 6
elif which_survey=="LSST":
	# covariance matrix for LSST
	infile_fmt = "LSST_Y%d_ssss_cov_Ncl15_Ntomo10"
	outfile_fmt = "LSST_Y%d_ssss_invcov_Ncl15_Ntomo10"
	Area_list = [12300, 14300]
	N_selection = 2
	years = [1, 10]
	N_tomo_list = [10, 10]
	Nsrc_list = np.array([11.112, 27.737])
	SN_list = [0.26,];N_shape_noise = 1
elif which_survey=="RomanPIT":
	# covariance matrix for Roman PIT
	infile_fmt = "Roman_WL_%d%d_ssss_cov_Ncl15_Ntomo10"
	outfile_fmt = "Roman_WL_%d%d_ssss_invcov_Ncl15_Ntomo10"
	Area = 1000
	N_depth = 5 # 5 number density bins
	N_tomo = 10
	Nsrc_list = np.array([20, 25, 30, 35, 40])
	# different ell_max
	ellmax_list = [1000, 2000, 3000, 4000];N_ellmax = 4
else:
	print(f'{which_survey} not supported')
	exit(1)

plot_corrmat = True

# for jSelect in range(N_selection):
# 	N_tomo = N_tomo_list[jSelect]
# 	Nsrc = Nsrc_list[jSelect]
# 	Area = Area_list[jSelect]
# 	for kSN in range(N_shape_noise):
# 		SN = SN_list[kSN]
for jSelect in range(N_depth):
	Nsrc = Nsrc_list[jSelect]
	for kSN in range(N_ellmax):
		ellmax = ellmax_list[kSN]
		fig_title = \
		"$n_\mathrm{src}=$%.0f arcmin$^{-2}$, $\ell_\mathrm{max}$=%.0f"%(Nsrc, ellmax)
		if which_survey=='DESI2':
			infile = DATA_DIR+"cov/"+infile_fmt%(jSelect,kSN,Ncl,N_tomo)
			outname= DATA_DIR+"invcov/"+outfile_fmt%(jSelect,kSN,Ncl,N_tomo)
		elif which_survey=='LSST':
			infile = DATA_DIR+"cov/"+infile_fmt%(years[jSelect])
			outname= DATA_DIR+"invcov/"+outfile_fmt%(years[jSelect])
		elif which_survey=="RomanPIT":
			infile = DATA_DIR+"cov/"+infile_fmt%(jSelect,kSN)
			outname= DATA_DIR+"invcov/"+outfile_fmt%(jSelect,kSN)

		if not os.path.exists(infile):
			print(f'File {infile} does not exist! Continue')
			continue
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
			# traditional cosmocov format, g/ng
			if covfile.shape[1]==10:
				cov_g, cov_ng = covfile[l,8], covfile[l,9]
			# break-down version
			elif covfile.shape[1]==12:
				cov_g  = covfile[l,8]+covfile[l,9]+covfile[l,10]
				cov_ng = covfile[l,11]
			# covariance only
			elif covfile.shape[1]==3:
				cov_g = 0; cov_ng = covfile[l,2]
			else:
				print(f'Wrong covariance matrix file format! {covfile.shape[1]} columns')
				exit(-1)
			if bin1==bin2:
				factor = 1
			else:
				factor = mask[bin1] * mask[bin2]
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

			plt.subplot(1, 1, 1, figsize=(10,10))
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
			if which_survey=="DESI":
				figfilename_fmt = "test_imgs/"+outfile_fmt%(jSelect,kSN,Ncl,N_tomo)+".png"
			elif which_survey=="LSST":
				figfilename_fmt = "test_imgs/"+outfile_fmt%(years[jSelect])+".png"
			elif which_survey=="RomanPIT":
				figfilename_fmt = "test_imgs/"+outfile_fmt%(jSelect,kSN)+".png"
			plt.savefig(figfilename_fmt, format="png")
			#plt.savefig(figfilename_fmt%(year), format="png")
			plt.close()
