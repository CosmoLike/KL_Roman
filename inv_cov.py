#!/usr/bin/python
import sys
import math, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import linalg as LA
import numpy as np

# covariance matrix
DATA_DIR = '/xdisk/timeifler/jiachuanxu/DESI2KL/'
infile_fmt = "DESI2_KL_%d%d%d_ssss_cov_Ncl%d_Ntomo%d"
Ncl = 15
N_area = 3
N_selection = 4
N_tomo_list = [4, 4, 4, 3]
N_shape_noise = 4
plot_corrmat = True

for i in range(N_area):
	for j in range(N_selection):
		N_tomo = N_tomo_list[j]
		for k in range(N_shape_noise):
			infile = DATA_DIR+"cov/"+infile_fmt%(i,j,k,Ncl,N_tomo)
			if not os.path.exists(infile):
				print(f'File {infile} does not exist! Continue')
				continue
			outname = DATA_DIR+"invcov/"+infile_fmt%(i,j,k,Ncl,N_tomo)

			# data vector
			data= ['/home/u17/jiachuanxu/CosmoLike/KL_WFIRST/datav/WFIRST_shear_shear_opti_Ntomo10_Ncl20_sigmae0.37_dmo_ia']

			nggl = 0 	# number of ggl power spectra
			ngcl = 0	# number of cluster-source galaxy power spectra
			nlens = 0 	# number of lens bins 
			nlenscl= 0 	# number of cluster redshift bins 
			nshear = (N_tomo+1)*N_tomo/2 # number of shear tomographic power spectra
			ncl=Ncl		# number of ell-bins
			nclgcl=0	# number of cluster ell-bins
			nrich=0		# number of richness bins


			ndata = (nshear+nggl+nlens)*ncl+nlenscl*nrich+nrich*ngcl*nclgcl
			n2pt = (nshear+nggl+nlens)*ncl 
			ncluster = nlenscl*nrich 
			n2ptcl=n2pt+ncluster
			nclusterN_WL=ncluster+nrich*ngcl*nclgcl
			print(f'Data Vector Size = {ndata}')
			mask = np.ones(ndata)
			#datafile= np.genfromtxt()
			#mask = np.zeros(ndata)
			#for i in range(0,datafile.shape[0]):
			#		if (datafile[i,1] >1.0e-15): 
			#			mask[i]=1.0
			#	for m in mask:
			#		m = 1.0
  	
  			covfile = np.genfromtxt(infile)
  			_ndata = int(np.max(covfile[:,0])+1)
  			assert _ndata == ndata, f'Inconsistent Ndata {ndata} v.s. {_ndata}'
			cov = np.zeros((ndata,ndata))

			for l in range(0,covfile.shape[0]):
				bin1, bin2 = int(covfile[l,0]), int(covfile[l,1])
				cov_g, cov_ng = covfile[l,8], covfile[l,9]
			  	cov[bin1, bin2] = cov_g+cov_ng
			  	cov[bin2, bin1] = cov_g+cov_ng
	 
	
			cor = cov/np.outer(LA.diag(cov)**0.5, LA.diag(cov)**0.5)

			# for i in range(0,ndata):
			#     for j in range(0,ndata):
			#     	if (cov[i,i]*cov[j,j] >0):
			#        		cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])

			a = np.sort(LA.eigvals(cor[:,:]))
			print("Eigenvalues range of corrmat (%.3e, %.3e)"%(np.min(a), np.max(a)))
			print("Negative eigenvalues of corrmat:{}".format(a[a<0]))
			#for i in range(0,a.shape[0]):
			#	if (a[i]< 0.0): print a[i]


			# ############### invert shear covariance #################
			inv = LA.inv(cov[0:nshear*ncl,0:nshear*ncl])
			a = np.sort(LA.eigvals(cov[0:nshear*ncl,0:nshear*ncl]))
			print("Eigvals range of covmat (%.3e, %.3e)"%(np.min(a), np.max(a)))
			f = open(outname, "w")
			for i in range(0,nshear*ncl):
			  	for j in range(0,nshear*ncl):
			  		f.write("%d %d %e\n" %(i,j, inv[i,j]*mask[i]*mask[j]))
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
			if plot_corrmat:
				labels = (
	r'$C^{\kappa \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
	# r'$C^{g \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
	# r'$C^{gg}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
	# r'$N^{c}\left(r,z_{\mathrm{c}_i}\right)$',
	# r'$C^{c \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',
				)
				ticks = np.zeros(2)
				ticks[1] = nshear*ncl
				# ticks[2] = (nshear+nggl)*ncl
				# ticks[3] = n2pt
				# ticks[4] = n2pt+ncluster
				# ticks[5] = ndata
				tickx = 0.5*(ticks[:-1]+ticks[1:])

				cor = cov/np.outer(LA.diag(cov)**0.5, LA.diag(cov)**0.5)
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
				im = ax.imshow(cor, interpolation='nearest', origin='lower', vmin=-1, vmax=1)
				plt.xticks(tickx, labels,fontsize=fs)
				plt.yticks(tickx-0.5, labels,fontsize=fs)
				#plt.tick_params(axis = 'x',length = 0, pad = 15)
				#plt.tick_params(axis = 'y',length = 0, pad = 5)

				plt.colorbar(im)
				#plt.show()
				plt.savefig("test_imgs/"+outname+".png", format="png")
				
				# plt.figure()
				# #plt.imshow(np.log10(cov[0:1500,2000:]), interpolation="nearest",vmin=-25, vmax=-10)
				# plt.imshow(np.log10(cov[:,:]), interpolation="nearest",vmin=-25, vmax=-10)
				# #plt.imshow(cor[n2ptcl:n2ptcl+200,300:nshear*ncl], interpolation="nearest",vmax=0.5)
				# plt.colorbar()
				# plt.show()

	




