#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include <fftw3.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>


#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/emu17/P_cb/emu.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift_spline.c"
#include "../cosmolike_core/theory/halo.c"
#include "../cosmolike_core/theory/HOD.c"
#include "../cosmolike_core/theory/pt.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/cluster.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "../cosmolike_core/theory/covariances_3D.c"
#include "../cosmolike_core/theory/covariances_fourier.c"
#include "../cosmolike_core/theory/covariances_cluster.c"
#include "init_WFIRST_forecasts.c"

double cov_G_shear_shear_tomo_one(double l, double delta_l, int z1, int z2, int z3, int z4);
void run_cov_N_N (char *OUTFILE, char *PATH, int nzc1, int nzc2,int start);
void run_cov_cgl_N (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int nzc2, int start);
void run_cov_cgl_cgl (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int start);
void run_cov_cgl_cgl_all (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster);
void run_cov_shear_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_shear_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);
void run_cov_ggl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_ggl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);
void run_cov_cl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_cl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);

void run_cov_ggl_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_shear_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_shear_shear_one(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);

double cov_G_shear_shear_tomo_one(double l, double delta_l, int z1, int z2, int z3, int z4){
  double C13, C14, C23, C24, N13 =0, N14=0, N23=0, N24=0;
  double fsky = survey.area/41253.0;
  /* one component power spectra */
  C13 = 0.25*C_shear_tomo_nointerp(l,z1,z3);C24 = 0.25*C_shear_tomo_nointerp(l,z2,z4);
  C14 = 0.25*C_shear_tomo_nointerp(l,z1,z4);C23 = 0.25*C_shear_tomo_nointerp(l,z2,z3);

  /* this is one component shot noise*/
  if (z1 == z3){N13= pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);}
  if (z1 == z4){N14= pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);}
  if (z2 == z3){N23= pow(survey.sigma_e,2.0)/(2.0*nsource(z2)*survey.n_gal_conversion_factor);}
  if (z2 == z4){N24=pow(survey.sigma_e,2.0)/(2.0*nsource(z2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+ C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23+N13*N24+N14*N23)/((2.*l+1.)*delta_l*fsky);
}

void run_cov_N_N (char *OUTFILE, char *PATH, int nzc1, int nzc2,int start)
{
  int nN1, nN2,i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      i += Cluster.N200_Nbin*nzc1+nN1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;

      cov =cov_N_N(nzc1,nN1, nzc2, nN2);
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,0.0,0.0, nzc1, nN1, nzc2, nN2,cov,0.0);
    }
  }
  fclose(F1);
}

void run_cov_cgl_N (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int nzc2, int start)
{
  int nN1, nN2, nl1, nzc1, nzs1,i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc1 = ZC(N1);
  nzs1 = ZSC(N1);
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
     for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
       i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
       i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
       j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
       j += Cluster.N200_Nbin*nzc2+nN2;

       cov =cov_cgl_N(ell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2);
       fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell_Cluster[nl1], 0., nzc1, nzs1, nzc2, nN2,cov,0.);
     }
   }
 }
 fclose(F1);
}

void run_cov_cgl_cgl (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int start)
{
  int nN1, nN2, nl1, nzc1, nzs1, nl2, nzc2, nzs2,i,j;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc1 = ZC(N1);
  nzs1 = ZSC(N1);
  nzc2 = ZC(N2);
  nzs2 = ZSC(N2);
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
      for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
        for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
          i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
          i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
          j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
          j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;

          c_g = 0;
          c_ng = cov_NG_cgl_cgl(ell_Cluster[nl1],ell_Cluster[nl2],nzc1,nN1, nzs1, nzc2, nN2,nzs2);
          if (nl2 == nl1){c_g =cov_G_cgl_cgl(ell_Cluster[nl1],dell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2,nzs2);}
          fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell_Cluster[nl1],ell_Cluster[nl2], nzc1, nzs1, nzc2, nzs2,c_g, c_ng);
        }
      }
    }
  }
  fclose(F1);
}

void run_cov_cgl_cgl_all (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster)
{
  int nN1, nN2, nl1, nzc1, nzs1, nl2, nzc2, nzs2,i,j, N1,N2;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s",PATH,OUTFILE);
  F1 =fopen(filename,"w");
  for (N1 = 0; N1 < tomo.cgl_Npowerspectra; N1 ++){
    for (N2 = 0; N2 < tomo.cgl_Npowerspectra; N2 ++){
      nzc1 = ZC(N1);
      nzs1 = ZSC(N1);
      nzc2 = ZC(N2);
      nzs2 = ZSC(N2);
      for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
        for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
          for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
            for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
              i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
              i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
              j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
              j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;         
              c_g = 0;
              c_ng = cov_NG_cgl_cgl(ell_Cluster[nl1],ell_Cluster[nl2],nzc1,nN1, nzs1, nzc2, nN2,nzs2);
              if (nl2 == nl1){c_g =cov_G_cgl_cgl(ell_Cluster[nl1],dell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2,nzs2);}
              fprintf(F1,"%d %d %e %e %d %d %d  %d %d %d  %e %e\n",i,j,ell_Cluster[nl1],ell_Cluster[nl2], nzc1, nN1,nzs1, nzc2, nN2,nzs2,c_g, c_ng);
            }
          }
        }
      }
    }
  }
  fclose(F1);
}

void run_cov_shear_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start)
{
  int nz1,nz2, nN2, nl1, nzc1, i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nz1 = Z1(N1);
  nz2 = Z2(N1);
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      cov = 0.;
      i = like.Ncl*N1+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;

      if (ell[nl1] < like.lmax_shear){cov =cov_shear_N(ell[nl1],nz1,nz2, nzc2, nN2);}
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., nz1, nz2, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_shear_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN1, nN2, nzs1, nzs2,nl2, nzc2, nzs3,i,j;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzs1 = Z1(N1);
  nzs2 = Z2(N1);
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*N1+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;          
        c_g = 0.;
        c_ng = 0.;
        if (ell[nl1] < like.lmax_shear){
          c_ng = cov_NG_shear_cgl(ell[nl1],ell_Cluster[nl2],nzs1, nzs2, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.001){ 
            c_g =cov_G_shear_cgl(ell[nl1],dell_Cluster[nl2],nzs1,nzs2, nzc2, nN2,nzs3);
          }
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], nzs1, nzs2, nzc2, nzs3,c_g, c_ng);
      }
    }
  }
  fclose(F1);
}

void run_cov_ggl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start)
{
  int zl,zs, nN2, nl1, nzc1, i,j;
  double cov,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(N1);
  zs = ZS(N1);
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
      cov = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight){
        cov =cov_ggl_N(ell[nl1],zl,zs, nzc2, nN2);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., zl, zs, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_ggl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN2, zl, zs, nzs1, nl2, nzc2, nzs3,i,j;
  double c_g, c_ng,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(N1);
  zs = ZS(N1);
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*(tomo.shear_Npowerspectra+N1)+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;

        c_g = 0; c_ng = 0.;
        weight = test_kmax(ell[nl1],zl);
        if (weight){
          c_ng = cov_NG_ggl_cgl(ell[nl1],ell_Cluster[nl2],zl,zs, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){c_g =cov_G_ggl_cgl(ell[nl1],dell_Cluster[nl2],zl,zs, nzc2, nN2,nzs3);}
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], zl, zs, nzc2, nzs3,c_g, c_ng);
      }
    }
  }
  fclose(F1);
}

void run_cov_cl_N (char *OUTFILE, char *PATH, double *ell, double *dell,int N1, int nzc2, int start)
{
  int zl,zs, nN2, nl1, nzc1, i,j;
  double cov,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
      cov = 0.;
      weight = test_kmax(ell[nl1],N1);
      if (weight){
        cov =cov_cl_N(ell[nl1],N1,N1,nzc2,nN2);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., N1, N1, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_cl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN2,nzc2, nzs3,i,j,nl2;
  double c_g, c_ng,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+N1)+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;

        c_g = 0; c_ng = 0.;
        weight = test_kmax(ell[nl1],N1);
        if (weight){
          c_ng = cov_NG_cl_cgl(ell[nl1],ell_Cluster[nl2],N1,N1, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){
            c_g =cov_G_cl_cgl(ell[nl1],dell_Cluster[nl2],N1,N1, nzc2, nN2,nzs3);
          }
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], N1,N1, nzc2, nzs3,c_g, c_ng);
        //printf("%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], N1,N1, nzc2, nzs3,c_g, c_ng);
      }
    }
  }
  fclose(F1);
}

void run_cov_ggl_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl,zs,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w"); 
  zl = ZL(n1); zs = ZS(n1);
  printf("\nN_ggl = %d (%d, %d)\n", n1,zl,zs);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight && ell[nl2] < like.lmax_shear){
        if (test_zoverlap(zl,z3)*test_zoverlap(zl,z4)){
          c_ng = cov_NG_gl_shear_tomo(ell[nl1],ell[nl2],zl,zs,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_gl_shear_tomo(ell[nl1],dell[nl1],zl,zs,z3,z4);
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],zl,zs,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],z1);
      if (weight && ell[nl2] < like.lmax_shear){
        if (test_zoverlap(z1,z3)*test_zoverlap(z1,z4)){
          c_ng = cov_NG_cl_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_cl_shear_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        }
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,zl,zs,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (z1 == zl){
        weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],zl);
        if (weight){
          c_ng = cov_NG_cl_gl_tomo(ell[nl1],ell[nl2],z1,z2,zl,zs);
          if (nl1 == nl2){
            c_g =  cov_G_cl_gl_tomo(ell[nl1],dell[nl1],z1,z2,zl,zs);
          }
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2, ell[nl1],ell[nl2],z1,z2,zl,zs,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  z3 = n2; z4 = n2;
  printf("N_cl_2 = %d\n", n2);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (z1 == z3){
        weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],z3);
        printf("ell[nl1]=%le ell[nl2]=%le test_kmax %d %d\n",ell[nl1],ell[nl2],test_kmax(ell[nl1],z1),test_kmax(ell[nl2],z3));
        if (weight) {
          c_ng = cov_NG_cl_cl_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_cl_cl_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl1,zl2,zs1,zs2,nl1,nl2, weight;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  zl1 = ZL(n1); zs1 = ZS(n1);
  printf("\nN_tomo_1 = %d (%d, %d)\n", n1,zl1,zs1);
  zl2 = ZL(n2); zs2 = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl2,zs2);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],zl1)*test_kmax(ell[nl2],zl2);
      if (weight && zl1 == zl2) {
        c_ng = cov_NG_gl_gl_tomo(ell[nl1],ell[nl2],zl1,zs1,zl2,zs2);
      }
      if (nl1 == nl2){
        c_g =  cov_G_gl_gl_tomo(ell[nl1],dell[nl1],zl1,zs1,zl2,zs2);
      }
      if (weight ==0 && n2 != n1){
        c_g = 0;
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2, ell[nl1],ell[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);   
    }
  }
  fclose(F1);
}


void run_cov_shear_shear(char *OUTFILE, char *PATH, double *ell, double *dell,int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  z1 = Z1(n1); z2 = Z2(n1);
  printf("N_shear = %d\n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n",n2,z3,z4);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (ell[nl1] < like.lmax_shear && ell[nl2] < like.lmax_shear){
        c_ng = cov_NG_shear_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
      }
      if (nl1 == nl2){
        c_g =  cov_G_shear_shear_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        if (ell[nl1] > like.lmax_shear && n1!=n2){c_g = 0.;} 
      }         
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",like.Ncl*n1+nl1,like.Ncl*(n2)+nl2,ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
      //printf("%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*n1+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_shear_shear_one(char *OUTFILE, char *PATH, double *ell, double *dell,int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  z1 = Z1(n1); z2 = Z2(n1);
  printf("N_shear = %d\n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n",n2,z3,z4);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      // if (ell[nl1] < like.lmax_shear && ell[nl2] < like.lmax_shear){
      //   c_ng = cov_NG_shear_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
      // }
      if (nl1 == nl2){
        c_g =  cov_G_shear_shear_tomo_one(ell[nl1],dell[nl1],z1,z2,z3,z4);
        if (ell[nl1] > like.lmax_shear && n1!=n2){c_g = 0.;} 
      }         
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",like.Ncl*n1+nl1,like.Ncl*(n2)+nl2,ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
      //printf("%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*n1+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}


int main(int argc, char** argv)
{
  int i,l,m,n,o,s,p,nl1,t,k;
  char OUTFILE[400],filename[400];
  
  // Setting Scenarios (survey area, source density, lens density)
  // As a comparison, LSST (12300 deg2) assumes 15 bins from ell=20 to 3000
  // Roman HLIS (5000 deg2) assumes 20 bins from ell=30 to 4000
  int N_scenarios_area = 1;
  double survey_area[1] = {30000.0};
  char survey_names[1][100] = {"SKA_WL"}; // DSA_allsky, SKA_WL, etc

  // IA model
  char ia_model[100] = "NLA_HF";          // NLA_HF (WL) or none (KL)

  // 1 if single component
  int one = 0;
  
  // Six sets of target selection criteria, each with different n(z)
  int N_scenarios_selection = 1;
  // Start with 4 source tomo bins 
  //int Ntomo_source[6] = {4, 4, 4, 4, 4, 4};
  int Ntomo_source[1] = {4};
  char dndz[1][100] = {"zdistris/zdistri_trecs_WL"};
  printf("%d target selection scenarios\n", N_scenarios_selection);

  // Six shape noise scenarios
  // Note that we do not include correlation between shape noise and target 
  // selection here.
  int N_scenarios_shape_noise = 1;
  printf("%d shape noise scenarios\n", N_scenarios_shape_noise);
  // Lens galaxies not used, set to random value
  float lens_density = 66.0;
  // Lens galaxies not used, set to random value
  int Ntomo_lens = 4;
  double Rmin_bias = 21.0; // not used 
  // 15 ell bins in Fourier space, from 20 to 3000
  int Nell = 15;
  double ell_min = 20.0;
  double ell_max = 3000.0;
  double ell_max_shear = 3000.0;
  // Now count how many scenarios
  //int N_scenarios = N_scenarios_selection * N_scenarios_shape_noise;
  int N_scenarios = 1;
  //double scenario_table[1][3]={ {2000.0, 8.0, 66.0} };

  int hit=atoi(argv[1]);
  Ntable.N_a=20;

  k=1;
  //set l-bins for shear, ggl, clustering, clusterWL
  double logdl=(log(ell_max)-log(ell_min))/Nell;
  double *ell, *dell, *ell_Cluster, *dell_Cluster;
  ell=create_double_vector(0,Nell-1);
  dell=create_double_vector(0,Nell-1);
  ell_Cluster=create_double_vector(0,Cluster.lbin-1);
  dell_Cluster=create_double_vector(0,Cluster.lbin-1);
  int j=0;
  printf("Ell array:\n");
  for(i=0;i<Nell;i++){
    ell[i]=exp(log(ell_min)+(i+0.5)*logdl);
    dell[i]=exp(log(ell_min)+(i+1)*logdl)-exp(log(ell_min)+(i*logdl));
    if(ell[i]<ell_max_shear) printf("%le\n",ell[i]);
    if(ell[i]>ell_max_shear){
      ell_Cluster[j]=ell[i];
      dell_Cluster[j]=dell[i];
      printf("%le %le\n",ell[i],ell_Cluster[j]);
      j++;
    }
  } 


  for(t=0;t<N_scenarios;t++){
    // int temp = t;
    // int i_selection = temp/N_scenarios_shape_noise;
    // temp -= i_selection * N_scenarios_shape_noise;
    // int i_shape_noise = temp;
    // temp -= i_shape_noise;
    // assert(temp==0);
    int i_selection = t;
    int i_SN = 0;

    //RUN MODE setup
    init_cosmo_runmode("halofit");
    init_binning_fourier(Nell, ell_min, ell_max, ell_max_shear, Rmin_bias, 
    Ntomo_source[i_selection], Ntomo_lens);
    char _surveyname[10];
    //sprintf(_surveyname, "DESI2_KL_%d%d", i_selection, i_shape_noise);
    sprintf(_surveyname, survey_names[i_selection]);
    char _photoz_prior[100];
    char _shearm_prior[100];
    sprintf(_photoz_prior, "spec_%s", survey_names[i_selection]);
    sprintf(_shearm_prior, "shear_%s", survey_names[i_selection]);
    init_priors_IA_bary(_photoz_prior, _shearm_prior,"none","none",
      // IA_flag, A, beta, eta, etaZ
      false, 3.0, 1.2, 3.8, 2.0, 
      // bary_flag, Q1, Q2, Q3
      false, 16, 1.9, 0.7);
    init_survey(_surveyname);
    // init survey name, area, n_gal, shape noise, magnitude limit, K-correction
    //sprintf(survey.name, "%s_%d%d", "DESI2_KL_v2", i_selection, i_shape_noise);
    //survey.area = survey_area;
    //survey.n_gal = source_density[i_selection];
    //survey.sigma_e = shape_noise_rms[i_shape_noise];
    // init source and lens n(z) and photo-z 
    init_galaxies(dndz[i_selection], 
      "zdistris/lens_LSSTY1", 
      "gaussian", "gaussian", "SN10");// the last arg is lens sample
    init_clusters();                  // not used if we don't have clusters
    init_IA(ia_model, "GAMA");
    init_probes("shear_shear");
    // sprintf(covparams.outdir, 
    //   "/xdisk/timeifler/jiachuanxu/DESI2KL/covpara_v2/");
    sprintf(covparams.outdir,
        "/home/u15/yhhuang/cosmology/dsa/cov/");

    printf("----------------------------------\n");  
    printf("area: %.2f n_source: %.2f n_lens: %.2f\n",
      survey.area,survey.n_gal,survey.n_lens);
    printf("----------------------------------\n");
    /******************************* START ************************************/
    /********************** cosmic shear - cosmic shear  **********************/
    if(like.shear_shear==1){
      // one component
      if (one == 1){
        sprintf(OUTFILE, "%s_ssss_cov_Ncl%d_Ntomo%d_OneComp",survey.name,like.Ncl,tomo.shear_Nbin);
        for (l=0; l<tomo.shear_Npowerspectra; l++){
          for (m=l; m<tomo.shear_Npowerspectra; m++){
            if (k==hit){
              printf("catch k=%d", hit);
              sprintf(filename, "%s%s_%d", covparams.outdir, OUTFILE, k);
              if (fopen(filename, "r") != NULL){
                printf("File %s already exist! Foreced not to overwrite!\n", filename);
                exit(1);
              }
              else{
                run_cov_shear_shear_one(OUTFILE, covparams.outdir, ell, dell, l, m, k);
                printf("Exit normally!\n");
                return 0;
              }
            }
            k = k+1;
          }
        }
      }
      // normal case
      else {
        sprintf(OUTFILE, "%s_ssss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
        for (l=0; l<tomo.shear_Npowerspectra; l++){
          for (m=l;m<tomo.shear_Npowerspectra; m++){
            if(k==hit){
              printf("catch k=%d\n", hit); 
              sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
              if (fopen(filename, "r") != NULL){
                printf("File %s already exist! Forced not to overwrite!\n",
                        filename);
                exit(1);
              }
              else {
                run_cov_shear_shear(OUTFILE,covparams.outdir,ell,dell,l,m,k);
                printf("Exit normally!\n");return 0;
              }
            }
            k=k+1;
          }
        }
      }
    }
    /**********************       ggl    -       ggl     **********************/
    if(like.shear_pos==1){
      sprintf(OUTFILE, "%s_lsls_cov_Ncl%d_Ntomo%d",
        survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.ggl_Npowerspectra; l++){
        for (m=l;m<tomo.ggl_Npowerspectra; m++){
          if(k==hit){
			      printf("catch k=%d\n", hit);
            sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
              run_cov_ggl(OUTFILE,covparams.outdir,ell,dell,l,m,k);
              printf("Exit normally!\n");return 0;
            }
          }
          k=k+1;
        }
      }
    }
    /**********************   clustering - clustering    **********************/
    if(like.pos_pos==1){
      sprintf(OUTFILE,"%s_llll_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.clustering_Npowerspectra; l++){ //auto bins only for now!
        for (m=l;m<tomo.clustering_Npowerspectra; m++){
          if(k==hit){ 
            printf("catch k=%d\n", hit); 
			      sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
              run_cov_clustering(OUTFILE,covparams.outdir,ell,dell,l,m,k);
              printf("Exit normally!\n");return 0;
            }
          }
          k=k+1;
        }
      }
    }
    /**********************   clustering - cosmic shear  **********************/
    if((like.pos_pos==1) && (like.shear_shear==1)){
      sprintf(OUTFILE,"%s_llss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.clustering_Npowerspectra; l++){
        for (m=0;m<tomo.shear_Npowerspectra; m++){
          if(k==hit){
			      printf("catch k=%d\n", hit);
            sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
              run_cov_clustering_shear(OUTFILE,covparams.outdir,ell,dell,l,m,k);
              printf("Exit normally!\n");return 0;
            }
          }
          k=k+1;
        }
      }
    }
    /**********************   clustering -       ggl     **********************/
    if((like.pos_pos==1) && (like.shear_pos==1)){
      sprintf(OUTFILE,"%s_llls_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.clustering_Npowerspectra; l++){
        for (m=0;m<tomo.ggl_Npowerspectra; m++){
          if(k==hit){
			      printf("catch k=%d\n", hit);
            sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
              run_cov_clustering_ggl(OUTFILE,covparams.outdir,ell,dell,l,m,k);
              printf("Exit normally!\n");return 0;
            }
          }
          k=k+1;
        }
      }
    }
    /**********************       ggl    - cosmic shear  **********************/
    if((like.shear_pos==1) && (like.shear_shear)){
      sprintf(OUTFILE,"%s_lsss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.ggl_Npowerspectra; l++){
        for (m=0;m<tomo.shear_Npowerspectra; m++){
          if(k==hit){
			      printf("catch k=%d\n", hit);
            sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
              run_cov_ggl_shear(OUTFILE,covparams.outdir,ell,dell,l,m,k);
              printf("Exit normally!\n");return 0;
            }
          }
          k=k+1;
        }
      }
    }
    //****************************** 
    //******cluster covariance****** 
    //******************************

    /******************** cluster counts - cluster counts *********************/
    if(like.clusterN==1){
      sprintf(OUTFILE,"%s_nn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.cluster_Nbin; l++){
        for (m=0;m<tomo.cluster_Nbin; m++){
          if(k==hit){
			      printf("catch k=%d\n", hit);
            sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
              run_cov_N_N (OUTFILE,covparams.outdir,l,m,k);
              printf("Exit normally!\n");return 0;
            }
          }
          k=k+1;
        }
      }
    }
    /******************* cluster lensing - cluster lensing ********************/
    if(like.clusterWL==1){
      sprintf(OUTFILE,"%s_cscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.cgl_Npowerspectra; l++){
        for (m=0;m<tomo.cgl_Npowerspectra; m++){
          if(k==hit){
			      printf("catch k=%d\n", hit);
            sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
              run_cov_cgl_cgl (OUTFILE,covparams.outdir,ell_Cluster,dell_Cluster,l,m,k);
              printf("Exit normally!\n");return 0;
            }
          }
          k=k+1;
        } 
      }
    }
    /******************** cluster lensing - cluster counts ********************/
    if((like.clusterWL==1) && (like.clusterN==1)){
      sprintf(OUTFILE,"%s_csn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.cgl_Npowerspectra; l++){
        for (m=0;m<tomo.cluster_Nbin; m++){
          if(k==hit){
			      printf("catch k=%d\n", hit);
            sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
              run_cov_cgl_N (OUTFILE,covparams.outdir,ell_Cluster,dell_Cluster,l,m,k);
              printf("Exit normally!\n");return 0;
            }
          }
          k=k+1;
        }
      }
    }
    /********************** cosmic shear - cluster counts *********************/
    if((like.shear_shear==1) && (like.clusterN==1)){
      printf(OUTFILE,"%s_ssn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.shear_Npowerspectra; l++){
        for (m=0;m<tomo.cluster_Nbin; m++){
          if(k==hit){
			      printf("catch k=%d\n", hit);
            sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
             run_cov_shear_N (OUTFILE,covparams.outdir,ell,dell,l,m,k);
             printf("Exit normally!\n");return 0;
           }
         }
         k=k+1;
       }
     }
   }
   /********************** cosmic shear - cluster lensing *********************/
   if((like.shear_shear==1) && (like.clusterWL==1)){
     sprintf(OUTFILE,"%s_sscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
     for (l=0;l<tomo.shear_Npowerspectra; l++){
        for (m=0;m<tomo.cgl_Npowerspectra; m++){
            if(k==hit){
			        printf("catch k=%d\n", hit);
              sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
              if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
              else {
                run_cov_shear_cgl (OUTFILE,covparams.outdir,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
                printf("Exit normally!\n");return 0;
              }
            }
            k=k+1;
        }
      }
    }
    /**********************     ggl - cluster counts     **********************/
    if((like.shear_pos==1) && (like.clusterN==1)){
      sprintf(OUTFILE,"%s_lsn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.ggl_Npowerspectra; l++){
        for (m=0;m<tomo.cluster_Nbin; m++){
          if(k==hit){
			      printf("catch k=%d\n", hit);
            sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
             run_cov_ggl_N (OUTFILE,covparams.outdir,ell,dell,l,m,k);
             printf("Exit normally!\n");return 0;
            }
          }
          k=k+1;
        }
      }
    }
    /**********************     ggl    - cluster lensing **********************/
    if((like.shear_pos==1) && (like.clusterWL==1)){
      sprintf(OUTFILE,"%s_lscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.ggl_Npowerspectra; l++){
        for (m=0;m<tomo.cgl_Npowerspectra; m++){
            if(k==hit){
			        printf("catch k=%d\n", hit);
              sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
              if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
              else {
                run_cov_ggl_cgl (OUTFILE,covparams.outdir,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
                printf("Exit normally!\n");return 0;
              }
            }       
            k=k+1;
        }
      }
    }
    /********************** clustering - cluster counts  **********************/
    if((like.pos_pos==1) && (like.clusterN==1)){
      sprintf(OUTFILE,"%s_lln_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.clustering_Npowerspectra; l++){
        for (m=0;m<tomo.cluster_Nbin; m++){
          if(k==hit){
			      printf("catch k=%d\n", hit);
            sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
            if (fopen(filename, "r") != NULL){
              printf("File %s already exist! Forced not to overwrite!\n",
                      filename);
              exit(1);
            }
            else {
             run_cov_cl_N (OUTFILE,covparams.outdir,ell,dell,l,m,k);
             printf("Exit normally!\n");return 0;
            }
          } 
          k=k+1;
        }
      }
    }
    /********************** clustering - cluster lensing **********************/
    if((like.pos_pos==1) && (like.clusterWL==1)){
      sprintf(OUTFILE,"%s_llcs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
      for (l=0;l<tomo.clustering_Npowerspectra; l++){
        for (m=0;m<tomo.cgl_Npowerspectra; m++){
            if(k==hit){
			        printf("catch k=%d\n", hit);
              sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
              if (fopen(filename, "r") != NULL){
				        printf("File %s already exist! Forced not to overwrite!\n", 
					       filename);
				        exit(1);
			        }
              else {
                run_cov_cl_cgl (OUTFILE,covparams.outdir,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
                printf("Exit normally!\n");return 0;
              }
            }
            k=k+1;
        }
      }
    }
    /******************************** END *************************************/
  }
  printf("number of cov blocks for parallelization: %d\n",k-1); 
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0; 
}

