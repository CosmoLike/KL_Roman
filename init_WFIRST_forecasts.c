#include <stdbool.h>
double invcov_read(int READ, int ci, int cj);
double data_read(int READ, int ci);
double bary_read(int READ, int PC, int cj);
void init_data_inv(char *INV_FILE, char *DATA_FILE);
void init_data_inv_bary(char *INV_FILE, char *DATA_FILE, char *BARY_FILE);
void init_priors(char *cosmoPrior1, char *cosmoPrior2, char *cosmoPrior3, char *cosmoPrior4);
void init_priors_IA(char *Prior1, char *Prior2, char *Prior3, char *Prior4, 
	double A_ia_Prior, double beta_ia_Prior, double eta_ia_Prior, double etaZ_ia_Prior);
void init_priors_IA_bary(char *Prior1, char *Prior2, char *Prior3, char *Prior4, 
	bool IA_flag, double A_ia_Prior, double beta_ia_Prior, double eta_ia_Prior, double etaZ_ia_Prior, 
	bool BARY_flag, double Q1_Prior, double Q2_Prior, double Q3_Prior);
void init_priors_KL(char *Prior1, char *Prior2, char *Prior3, char *Prior4);
void init_survey(char *surveyname);
void init_galaxies(char *SOURCE_ZFILE, char *LENS_ZFILE, char *lensphotoz, char *sourcephotoz, char *galsample);
void init_cosmo_runmode(char *runmode);
void init_cosmo_runmode_DEu95CPL(char *runmode);
void init_cosmo_runmode_DEl95CPL(char *runmode);
void init_binning_fourier(int Ncl, double lmin, double lmax, double lmax_shear, double Rmin_bias, int Ntomo_source, int Ntomo_lens);
void init_probes(char *probes);
void init_clusters();

void set_galaxies_source();
void set_galaxies_SN10();
void set_galaxies_redmagic();
void set_clusters_WFIRST(); //set parameters for LSST/WFIRST forecasts
void init_lens_sample(char *lensphotoz, char *galsample);
void init_source_sample(char *sourcephotoz);


void set_wlphotoz_WFIRST_opti();
void set_clphotoz_WFIRST_opti();
void set_wlphotoz_WFIRST_KL();
void set_clphotoz_WFIRST_KL();
void set_wlphotoz_WFIRST_pessi();
void set_clphotoz_WFIRST_pessi();
void set_wlphotoz_DESI2_KL();
void set_clphotoz_DESI2_KL();
void set_wlphotoz_LSST_Y1();
void set_clphotoz_LSST_Y1();
void set_wlphotoz_LSST_Y10();
void set_clphotoz_LSST_Y10();
void set_wlphotoz_DSA_allsky();
void set_clphotoz_DSA_allsky();

void set_shear_priors_WFIRST_KL();
void set_shear_priors_DESI2_KL();
void set_shear_priors_WFIRST_opti();
void set_shear_priors_WFIRST_pessi();
void set_shear_priors_LSST_Y1();
void set_shear_priors_LSST_Y10();
void set_shear_priors_DSA_allsky();

void set_survey_parameters_to_WFIRST_WL();
void set_survey_parameters_to_WFIRST_KL();
void set_survey_parameters_to_DESI2_KL(char *surveyname);
void set_survey_parameters_to_LSST_Y1();
void set_survey_parameters_to_LSST_Y10();
void set_survey_parameters_to_DSA_allsky();

void init_clusterMobs();
void set_equal_tomo_bins();
void init_IA(char *model,char *lumfct);

void set_galaxies_DES_Y1();

double log_L_PlanckBAOJLA_w0wa();

int count_rows(char* filename,const char delimiter){
  FILE *file = fopen (filename, "r" );
  char line [1000];
  if (file != NULL) {
    fgets(line,sizeof line,file);
    fclose(file);
  }
  else{
    printf("count_rows: file %s not found.\nEXIT\n",filename);
    exit(1);
  }
  int count = 1;
  char *p;

  p = line;
  while (*p != '\0')
  {
    if (*p == delimiter){
        while (*p == delimiter){p++;}
        count++;
    }
      p++;
    }
   return count;
}

double invcov_read(int READ, int ci, int cj)
{
  int i,j,intspace;
  static double **inv =0;

  if(READ==0 || inv == 0){
    inv   = create_double_matrix(0, like.Ndata-1, 0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.INV_FILE,"r");
    for (i=0;i<like.Ndata; i++){
      for (j=0;j<like.Ndata; j++){
       fscanf(F,"%d %d %le\n",&intspace,&intspace,&inv[i][j]);  
     }
   }
   fclose(F);
   printf("FINISHED READING COVARIANCE\n");
 }    
 return inv[ci][cj];
}


double data_read(int READ, int ci)
{
  int i,intspace;
  static double *data = 0;
  
  if(READ==0 || data ==0){
    data  = create_double_vector(0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.DATA_FILE,"r");
    for (i=0;i<like.Ndata; i++){  
      fscanf(F,"%d %le\n",&intspace,&data[i]);
    }
    fclose(F);
    printf("FINISHED READING DATA VECTOR\n");
  }    
  return data[ci];
}

double bary_read(int READ, int PC, int cj)
{
  int i,j, N_PC=6;
  static double **bary =0;

  if(READ==0 || bary == 0){
    bary   = create_double_matrix(0, N_PC-1, 0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.BARY_FILE,"r");
    for (i=0;i<like.Ndata; i++){
      fscanf(F,"%le %le %le %le %le %le\n",&bary[0][i],&bary[1][i],&bary[2][i],&bary[3][i],&bary[4][i],&bary[5][i]);  
    }
    fclose(F);
    printf("FINISHED READING BARYON MATRIX\n");
  }    
  return bary[PC][cj];
}

void init_data_inv_bary(char *INV_FILE, char *DATA_FILE, char *BARY_FILE)
{
  double init;
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing data vector and covariance\n");
  printf("---------------------------------------\n");

  sprintf(like.INV_FILE,"%s",INV_FILE);
  printf("PATH TO INVCOV: %s\n",like.INV_FILE);
  sprintf(like.DATA_FILE,"%s",DATA_FILE);
  printf("PATH TO DATA: %s\n",like.DATA_FILE);
  sprintf(like.BARY_FILE,"%s",BARY_FILE);
  printf("PATH TO BARYONS: %s\n",like.BARY_FILE);
  init=data_read(0,1);
  init=bary_read(0,1,1);
  init=invcov_read(0,1,1);
}

void init_cosmo_runmode(char *runmode)
{
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Initializing Standard Runmode/Cosmology\n");
  printf("-------------------------------------------\n");
  set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP();
  set_prior_cosmology(); // if external data sets are included, this routine centers the external data set on top of the fiducial model for the simulated analysis.

  sprintf(pdeltaparams.runmode,"%s",runmode);
  printf("pdeltaparams.runmode =%s\n",pdeltaparams.runmode);
}

void init_cosmo_runmode_DEu95CPL(char *runmode)
{
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Initializing Standard Runmode/Cosmology\n");
  printf("-------------------------------------------\n");
  set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP_DEu95CPL();
  set_prior_cosmology(); // if external data sets are included, this routine centers the external data set on top of the fiducial model for the simulated analysis.

  sprintf(pdeltaparams.runmode,"%s",runmode);
  printf("pdeltaparams.runmode =%s\n",pdeltaparams.runmode);
}

void init_cosmo_runmode_DEl95CPL(char *runmode)
{
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Initializing Standard Runmode/Cosmology\n");
  printf("-------------------------------------------\n");
  set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP_DEl95CPL();
  set_prior_cosmology(); // if external data sets are included, this routine centers the external data set on top of the fiducial model for the simulated analysis.

  sprintf(pdeltaparams.runmode,"%s",runmode);
  printf("pdeltaparams.runmode =%s\n",pdeltaparams.runmode);
}


void init_binning_fourier(int Ncl, double lmin, double lmax, double lmax_shear, double Rmin_bias, int Ntomo_source, int Ntomo_lens)
{
  printf("-------------------------------------------\n");
  printf("Initializing Binning\n");
  printf("-------------------------------------------\n");
  
  like.Rmin_bias=Rmin_bias;
  like.Ncl=Ncl;
  like.lmin= lmin; //std=20
  like.lmax= lmax; //15,000
  like.lmax_shear = lmax_shear; //5000
  tomo.shear_Nbin=Ntomo_source;
  tomo.clustering_Nbin=Ntomo_lens;
  //compute cluster ell bins acc to 2PCF l-bins
  double ell;
  int i,k=0;
  double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
  for(i=0;i<like.Ncl;i++){
    ell=exp(log(like.lmin)+(i+0.5)*logdl);
    if (ell > like.lmax_shear){
      if (k==0) Cluster.l_min = ell;
      k=k+1;
    }
  } 
  Cluster.lbin = k;
  Cluster.l_max = lmax; //clusters go to highly nonlin as std
  //printf("%le %le %d\n",Cluster.l_min,Cluster.l_max,Cluster.lbin);
  like.lmax_kappacmb = 2999.;
  
  printf("number of ell bins Ncl: %d\n",like.Ncl);
  printf("minimum ell: %le\n",like.lmin);
  printf("maximum ell: %le\n",like.lmax);
}


void init_priors(char *Prior1, char *Prior2, char *Prior3, char *Prior4)
{
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing priors for marginalization\n");
  printf("---------------------------------------\n");
    
  if(strcmp(Prior1,"photo_opti")==0){
    set_wlphotoz_WFIRST_opti();
    set_clphotoz_WFIRST_opti();
  }
  if(strcmp(Prior1,"photo_pessi")==0){
    set_wlphotoz_WFIRST_pessi();
    set_clphotoz_WFIRST_pessi();
  }
  if(strcmp(Prior2,"shear_opti")==0){
    set_shear_priors_WFIRST_opti();
  }
  if(strcmp(Prior2,"shear_pessi")==0){
    set_shear_priors_WFIRST_pessi();
  }
  if(strcmp(Prior3,"GRS")==0) like.GRS=1;
  // external probes
  sprintf(like.ext_data,"%s",Prior4);
  if(strcmp(Prior4, "Planck15_BAO_H070p6_JLA_w0wa")==0){
    printf("Set external prior to Planck+BAO+JLA (w0-wa)\n");
    like.Planck15_BAO_H070p6_JLA_w0wa=1;  
  } 
  // else if(strcmp(Prior4,"Planck18_BAO_Riess18_Pantheon_w0wa")==0){
  //   like.Planck18_BAO_Riess18_Pantheon_w0wa=1;  
  // } 
  // else if(strcmp(Prior4,"Planck18_BAO_w0wa")==0){
  //   like.Planck18_BAO_w0wa=1; 
  // } 
  // else if(strcmp(Prior4,"Planck18_w0")==0){
  //   like.Planck18_w0=1;          
  // } 
  else {
    printf("Error from like_fourier.c: Prior4 can only be Planck, Planck15_BAO_w0wa, Planck15_BAO_H070p6_JLA_w0wa, Planck18_BAO_Riess18_Pantheon_w0wa, Planck18_BAO_w0wa, or Planck18_w0."); //CH: no real error handling.
  }
}

void init_priors_IA_bary(char *Prior1, char *Prior2, char *Prior3, char *Prior4, 
	bool IA_flag, double A_ia_Prior, double beta_ia_Prior, double eta_ia_Prior, double etaZ_ia_Prior,
	bool BARY_flag, double Q1_Prior, double Q2_Prior, double Q3_Prior)
{
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing priors for marginalization\n");
  printf("---------------------------------------\n");
  // Initializing photo-z priors  
  if(strcmp(Prior1,"photo_opti")==0){
    set_wlphotoz_WFIRST_opti();
    set_clphotoz_WFIRST_opti();
  }
  if(strcmp(Prior1,"photo_pessi")==0){
    set_wlphotoz_WFIRST_pessi();
    set_clphotoz_WFIRST_pessi();
  }
  if(strcmp(Prior1,"photo_LSST_Y1")==0){
    set_wlphotoz_LSST_Y1();
    set_clphotoz_LSST_Y1();
  }
  if(strcmp(Prior1,"photo_LSST_Y10")==0){
    set_wlphotoz_LSST_Y10();
    set_clphotoz_LSST_Y10();
  }
  if(strcmp(Prior1,"spec_WFIRST")==0){
    set_wlphotoz_WFIRST_KL();
    set_clphotoz_WFIRST_KL();
  }
  if(strcmp(Prior1,"spec_DESI2")==0){
    set_wlphotoz_DESI2_KL();
    set_clphotoz_DESI2_KL();
  }
  if(strcmp(Prior1,"spec_DSA_allsky")==0){
    set_wlphotoz_DSA_allsky();
    set_clphotoz_DSA_allsky();
  }
  if(strcmp(Prior1,"spec_SKA_WL")==0){
    set_wlphotoz_SKA_WL();
    set_clphotoz_SKA_WL();
  }
  // Initializing shear calibration priors
  if(strcmp(Prior2,"shear_opti")==0){
    set_shear_priors_WFIRST_opti();
  }
  if(strcmp(Prior2,"shear_pessi")==0){
    set_shear_priors_WFIRST_pessi();
  }
  if(strcmp(Prior2,"shear_LSST_Y1")==0){
    set_shear_priors_LSST_Y1();
  }
  if(strcmp(Prior2,"shear_LSST_Y10")==0){
    set_shear_priors_LSST_Y10();
  }
  if(strcmp(Prior2,"shear_KL_WFIRST")==0){
    set_shear_priors_WFIRST_KL();
  }
  if(strcmp(Prior2,"shear_KL_DESI2")==0){
    set_shear_priors_DESI2_KL();
  }
  if(strcmp(Prior2,"shear_DSA_allsky")==0){
    set_shear_priors_DSA_allsky();
  }
  if(strcmp(Prior2,"shear_SKA_WL")==0){
    set_shear_priors_SKA_WL();
  }
  if(strcmp(Prior3,"GRS")==0) like.GRS=1;
  // external probes
  sprintf(like.ext_data,"%s",Prior4);
  if(strcmp(Prior4, "Planck15_BAO_H070p6_JLA_w0wa")==0){
    printf("Set external prior to Planck+BAO+JLA (w0-wa)\n");
    like.Planck15_BAO_H070p6_JLA_w0wa=1;  
  } 
  // else if(strcmp(Prior4,"Planck18_BAO_Riess18_Pantheon_w0wa")==0){
  //   like.Planck18_BAO_Riess18_Pantheon_w0wa=1;  
  // } 
  // else if(strcmp(Prior4,"Planck18_BAO_w0wa")==0){
  //   like.Planck18_BAO_w0wa=1; 
  // } 
  // else if(strcmp(Prior4,"Planck18_w0")==0){
  //   like.Planck18_w0=1;          
  // } 
  else {
    printf("Error from like_fourier.c: Prior4 can only be Planck, Planck15_BAO_w0wa, Planck15_BAO_H070p6_JLA_w0wa, Planck18_BAO_Riess18_Pantheon_w0wa, Planck18_BAO_w0wa, or Planck18_w0."); //CH: no real error handling.
  }
  if(IA_flag){
    // Initializing IA priors
    prior.A_ia[0]=5.92; 
    prior.A_ia[1]=A_ia_Prior; 
    
    prior.beta_ia[0]=1.1; 
    prior.beta_ia[1]=beta_ia_Prior; 
    
    prior.eta_ia[0]=-0.47; 
    prior.eta_ia[1]=eta_ia_Prior; 
    
    prior.eta_ia_highz[0]=0.0; 
    prior.eta_ia_highz[1]=etaZ_ia_Prior; 
    like.IA=1;
  }
  // Initializing baryon effects priors
  if(BARY_flag){
    prior.bary_Q1[0]=0.0; 
    prior.bary_Q1[1]=Q1_Prior; 

    prior.bary_Q2[0]=0.0; 
    prior.bary_Q2[1]=Q2_Prior; 

    prior.bary_Q3[0]=0.0; 
    prior.bary_Q3[1]=Q3_Prior; 
    like.baryons=1;
  }
  // Confirmation
  printf("\n");
  printf("---------------------------------------\n");
  printf("Shear Calibration Prior\n");
  printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[0][0],prior.shear_calibration_m[0][1]);

  printf("\n");
  printf("---------------------------------------\n");
  printf("Photo-z priors Weak Lensing\n");
  printf("Delta_z=%le, Sigma (Delta_z)=%le\n",prior.bias_zphot_shear[0][0],prior.bias_zphot_shear[0][1]);
  printf("Sigma_z=%le, Sigma (Sigma_z)=%le\n",prior.sigma_zphot_shear[0][0],prior.sigma_zphot_shear[0][1]); 

  printf("\n");
  printf("---------------------------------------\n");
  printf("Photo-z priors Clustering\n");
  printf("Delta_z=%le, Sigma (Delta_z)=%le\n",prior.bias_zphot_clustering[0][0],prior.bias_zphot_clustering[0][1]);
  printf("Sigma_z=%le, Sigma (Sigma_z)=%le\n",prior.sigma_zphot_clustering[0][0],prior.sigma_zphot_clustering[0][1]);
  if(IA_flag){
    printf("\n");
    printf("---------------------------------------\n");
    printf("IA Priors\n");
    printf("A_IA=%le, A_IA_Prior=%le\n",prior.A_ia[0],prior.A_ia[1]);
    printf("beta_ia=%le, betaIA_Prior=%le\n",prior.beta_ia[0],prior.beta_ia[1]);
    printf("eta_ia=%le, etaIA_Prior=%le\n",prior.eta_ia[0],prior.eta_ia[1]);
    printf("eta_ia_highz=%le, etaZIA_Prior=%le\n",prior.eta_ia_highz[0],prior.eta_ia_highz[1]);
  }
  if(BARY_flag){
    printf("\n");
    printf("---------------------------------------\n");
    printf("Baryon Priors\n");
    printf("Q1=%le, Sigma (Q1)=%le\n",prior.bary_Q1[0],prior.bary_Q1[1]);
    printf("Q2=%le, Sigma (Q2)=%le\n",prior.bary_Q2[0],prior.bary_Q2[1]);
    printf("Q3=%le, Sigma (Q3)=%le\n",prior.bary_Q3[0],prior.bary_Q3[1]);
  }
}

void init_priors_IA(char *Prior1, char *Prior2, char *Prior3, char *Prior4, 
	double A_ia_Prior, double beta_ia_Prior, double eta_ia_Prior, double etaZ_ia_Prior)
{
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing priors for marginalization\n");
  printf("---------------------------------------\n");
  // Initializing photo-z priors  
  if(strcmp(Prior1,"photo_opti")==0){
    set_wlphotoz_WFIRST_opti();
    set_clphotoz_WFIRST_opti();
  }
  if(strcmp(Prior1,"photo_pessi")==0){
    set_wlphotoz_WFIRST_pessi();
    set_clphotoz_WFIRST_pessi();
  }
  if(strcmp(Prior1,"spec_WFIRST")==0){
    set_wlphotoz_WFIRST_KL();
    set_clphotoz_WFIRST_KL();
  }
  if(strcmp(Prior1,"spec_DESI2")==0){
    set_wlphotoz_DESI2_KL();
    set_clphotoz_DESI2_KL();
  }
  if(strcmp(Prior1,"spec_DSA_allsky")==0){
    set_wlphotoz_DSA_allsky();
    set_clphotoz_DSA_allsky();
  }
  if(strcmp(Prior1,"spec_SKA_WL")==0){
    set_wlphotoz_SKA_WL();
    set_clphotoz_SKA_WL();
  }
  // Initializing shear calibration priors
  if(strcmp(Prior2,"shear_opti")==0){
    set_shear_priors_WFIRST_opti();
  }
  if(strcmp(Prior2,"shear_pessi")==0){
    set_shear_priors_WFIRST_pessi();
  }
  if(strcmp(Prior2,"shear_KL_WFIRST")==0){
    set_shear_priors_WFIRST_KL();
  }
  if(strcmp(Prior2,"shear_KL_DESI2")==0){
    set_shear_priors_DESI2_KL();
  }
  if(strcmp(Prior2,"shear_DSA_allsky")==0){
    set_shear_priors_DSA_allsky();
  }
  if(strcmp(Prior2,"shear_SKA_WL")==0){
    set_shear_priors_SKA_WL();
  }
  if(strcmp(Prior3,"GRS")==0) like.GRS=1;
  // external probes
  sprintf(like.ext_data,"%s",Prior4);
  if(strcmp(Prior4, "Planck15_BAO_H070p6_JLA_w0wa")==0){
    printf("Set external prior to Planck+BAO+JLA (w0-wa)\n");
    like.Planck15_BAO_H070p6_JLA_w0wa=1;  
  } 
  // else if(strcmp(Prior4,"Planck18_BAO_Riess18_Pantheon_w0wa")==0){
  //   like.Planck18_BAO_Riess18_Pantheon_w0wa=1;  
  // } 
  // else if(strcmp(Prior4,"Planck18_BAO_w0wa")==0){
  //   like.Planck18_BAO_w0wa=1; 
  // } 
  // else if(strcmp(Prior4,"Planck18_w0")==0){
  //   like.Planck18_w0=1;          
  // } 
  else {
    printf("Error from like_fourier.c: Prior4 can only be Planck, Planck15_BAO_w0wa, Planck15_BAO_H070p6_JLA_w0wa, Planck18_BAO_Riess18_Pantheon_w0wa, Planck18_BAO_w0wa, or Planck18_w0."); //CH: no real error handling.
  }
  // Initializing IA priors
  prior.A_ia[0]=5.92; 
  prior.A_ia[1]=A_ia_Prior; 
  
  prior.beta_ia[0]=1.1; 
  prior.beta_ia[1]=beta_ia_Prior; 
  
  prior.eta_ia[0]=-0.47; 
  prior.eta_ia[1]=eta_ia_Prior; 
  
  prior.eta_ia_highz[0]=0.0; 
  prior.eta_ia_highz[1]=etaZ_ia_Prior; 
  like.IA=1;
  // Confirmation
  printf("\n");
  printf("---------------------------------------\n");
  printf("Shear Calibration Prior\n");
  printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[0][0],prior.shear_calibration_m[0][1]);

  printf("\n");
  printf("---------------------------------------\n");
  printf("Photo-z priors Weak Lensing\n");
  printf("Delta_z=%le, Sigma (Delta_z)=%le\n",prior.bias_zphot_shear[0][0],prior.bias_zphot_shear[0][1]);
  printf("Sigma_z=%le, Sigma (Sigma_z)=%le\n",prior.sigma_zphot_shear[0][0],prior.sigma_zphot_shear[0][1]); 

  printf("\n");
  printf("---------------------------------------\n");
  printf("Photo-z priors Clustering\n");
  printf("Delta_z=%le, Sigma (Delta_z)=%le\n",prior.bias_zphot_clustering[0][0],prior.bias_zphot_clustering[0][1]);
  printf("Sigma_z=%le, Sigma (Sigma_z)=%le\n",prior.sigma_zphot_clustering[0][0],prior.sigma_zphot_clustering[0][1]);

  printf("\n");
  printf("---------------------------------------\n");
  printf("IA Priors\n");
  printf("A_IA=%le, A_IA_Prior=%le\n",prior.A_ia[0],prior.A_ia[1]);
  printf("beta_ia=%le, betaIA_Prior=%le\n",prior.beta_ia[0],prior.beta_ia[1]);
  printf("eta_ia=%le, etaIA_Prior=%le\n",prior.eta_ia[0],prior.eta_ia[1]);
  printf("eta_ia_highz=%le, etaZIA_Prior=%le\n",prior.eta_ia_highz[0],prior.eta_ia_highz[1]);
}

void init_priors_KL(char *Prior1, char *Prior2, char *Prior3, char *Prior4)
{
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing priors for marginalization\n");
  printf("---------------------------------------\n");
    
  if(strcmp(Prior1,"spec_WFIRST")==0){
    set_wlphotoz_WFIRST_KL();
    set_clphotoz_WFIRST_KL();
  }
  if(strcmp(Prior1,"spec_DESI2")==0){
    set_wlphotoz_DESI2_KL();
    set_clphotoz_DESI2_KL();
  }
  if(strcmp(Prior1,"spec_DSA_allsky")==0){
    set_wlphotoz_DSA_allsky();
    set_clphotoz_DSA_allsky();
  }
  if(strcmp(Prior2,"shear_KL_WFIRST")==0){
    set_shear_priors_WFIRST_KL();
  }
  if(strcmp(Prior2,"shear_KL_DESI2")==0){
    set_shear_priors_DESI2_KL();
  }
  if(strcmp(Prior2,"shear_KL_DSA")==0){
    set_shear_priors_DSA_allsky();
  }

  if(strcmp(Prior3,"GRS")==0) like.GRS=1;
  // external probes
  sprintf(like.ext_data,"%s",Prior4);
  if(strcmp(Prior4, "Planck15_BAO_H070p6_JLA_w0wa")==0){
    printf("Set external prior to Planck+BAO+JLA (w0-wa)\n");
    like.Planck15_BAO_H070p6_JLA_w0wa=1;  
  } 
  // else if(strcmp(Prior4,"Planck18_BAO_Riess18_Pantheon_w0wa")==0){
  //   like.Planck18_BAO_Riess18_Pantheon_w0wa=1;  
  // } 
  // else if(strcmp(Prior4,"Planck18_BAO_w0wa")==0){
  //   like.Planck18_BAO_w0wa=1; 
  // } 
  // else if(strcmp(Prior4,"Planck18_w0")==0){
  //   like.Planck18_w0=1;          
  // } 
  else {
    printf("Error from like_fourier.c: Prior4 can only be Planck, Planck15_BAO_w0wa, Planck15_BAO_H070p6_JLA_w0wa, Planck18_BAO_Riess18_Pantheon_w0wa, Planck18_BAO_w0wa, or Planck18_w0."); //CH: no real error handling.
  }
}


void init_survey(char *surveyname)
{
  printf("\n");
  printf("-------------------------------\n");
  printf("Initializing Survey Parameters\n");
  printf("-------------------------------\n");

  if(strcmp(surveyname,"LSST")==0) set_survey_parameters_to_LSST();
  if(strcmp(surveyname,"LSST_Y1")==0) set_survey_parameters_to_LSST_Y1();
  if(strcmp(surveyname,"LSST_Y10")==0) set_survey_parameters_to_LSST_Y10();
  if(strcmp(surveyname,"Euclid")==0) set_survey_parameters_to_Euclid();
  if(strcmp(surveyname,"WFIRST")==0) set_survey_parameters_to_WFIRST();
  if(strcmp(surveyname,"WFIRST_WL")==0) set_survey_parameters_to_WFIRST_WL();
  if(strcmp(surveyname,"WFIRST_KL")==0) set_survey_parameters_to_WFIRST_KL();
  if(strncmp(surveyname,"DESI2_KL",8)==0) set_survey_parameters_to_DESI2_KL(surveyname);
  if(strcmp(surveyname,"DSA_allsky")==0) set_survey_parameters_to_DSA_allsky();
  if(strcmp(surveyname,"SKA_WL")==0) set_survey_parameters_to_SKA_WL();

  printf("Survey set to %s\n",survey.name);
  printf("Survey area: %f deg^2\n",survey.area);
  printf("Source Galaxy Density: %f galaxies/arcmin^2\n",survey.n_gal);
  printf("Shape noise rms: %f", survey.sigma_e);
}


void init_galaxies(char *SOURCE_ZFILE, char *LENS_ZFILE, char *lensphotoz, char *sourcephotoz, char *galsample)
{
  printf("\n");
  printf("-----------------------------------\n");
  printf("Initializing galaxy samples\n");
  printf("-----------------------------------\n");
  
  sprintf(redshift.shear_REDSHIFT_FILE,"%s",SOURCE_ZFILE);
  printf("PATH TO SOURCE_ZFILE: %s\n",redshift.shear_REDSHIFT_FILE);
  init_source_sample(sourcephotoz);
  
  sprintf(redshift.clustering_REDSHIFT_FILE,"%s",LENS_ZFILE);
  printf("\n");
  printf("PATH TO LENS_ZFILE: %s\n",redshift.clustering_REDSHIFT_FILE);
  init_lens_sample(lensphotoz,galsample);
}

void init_clusters()
{
  printf("\n");
  printf("-----------------------------------\n");
  printf("Initializing clusters\n");
  printf("-----------------------------------\n");

  set_clusters_WFIRST();
  prior.cluster_Mobs_lgM0[0]=nuisance.cluster_Mobs_lgN0;
  prior.cluster_Mobs_lgM0[1]=0.045;

  prior.cluster_Mobs_alpha[0]=nuisance.cluster_Mobs_alpha;
  prior.cluster_Mobs_alpha[1]=0.045;

  prior.cluster_Mobs_beta[0]=nuisance.cluster_Mobs_beta;
  prior.cluster_Mobs_beta[1]=0.3;

  prior.cluster_Mobs_sigma0[0]=nuisance.cluster_Mobs_sigma0;
  prior.cluster_Mobs_sigma0[1]=0.045;
  
  prior.cluster_Mobs_sigma_qm[0]=nuisance.cluster_Mobs_sigma_qm;
  prior.cluster_Mobs_sigma_qm[1]=0.03;
  
  prior.cluster_Mobs_sigma_qz[0]=nuisance.cluster_Mobs_sigma_qz; 
  prior.cluster_Mobs_sigma_qz[1]=0.1;

  like.clusterMobs=1;

}


void init_probes(char *probes)
{
  printf("\n");
  printf("------------------------------\n");
  printf("Initializing Probes\n");
  printf("------------------------------\n"); 
  printf("tomo.cluster_Nbin=%d\n",tomo.cluster_Nbin);
  printf("tomo.cgl_Npowerspectra=%d\n",tomo.cgl_Npowerspectra);
  printf("Cluster.lbin=%d\n",Cluster.lbin);
  printf("Cluster.N200_Nbin=%d\n",Cluster.N200_Nbin);
  printf("like.Ncl=%d\n",like.Ncl);
  printf("tomo.shear_Npowerspectra=%d\n",tomo.shear_Npowerspectra);
  printf("tomo.ggl_Npowerspectra=%d\n",tomo.ggl_Npowerspectra);
  printf("tomo.clustering_Npowerspectra=%d\n",tomo.clustering_Npowerspectra);

  sprintf(like.probes,"%s",probes);
  if(strcmp(probes,"clusterN")==0){
    like.Ndata=tomo.cluster_Nbin*Cluster.N200_Nbin;
    like.clusterN=1;
    printf("Cluster Number Counts computation initialized\n");
  }
  if(strcmp(probes,"clusterN_clusterWL")==0){
    like.Ndata=tomo.cluster_Nbin*Cluster.N200_Nbin+tomo.cgl_Npowerspectra*Cluster.N200_Nbin*Cluster.lbin;
    like.clusterN=1;
    like.clusterWL=1;
    printf("Cluster Number Counts computation initialized\n");
    printf("Cluster weak lensing computation initialized\n");
  }
  if(strcmp(probes,"3x2pt_clusterN")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+tomo.cluster_Nbin*Cluster.N200_Nbin;
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    like.clusterN=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("Cluster Number Counts computation initialized\n");
  }

  if(strcmp(probes,"shear_shear")==0){
    like.Ndata=like.Ncl*tomo.shear_Npowerspectra;
    like.shear_shear=1;
    printf("Shear-Shear computation initialized\n");
  }
  if(strcmp(probes,"pos_pos")==0){
    like.Ndata= like.Ncl*tomo.clustering_Npowerspectra;
    like.pos_pos=1;
    printf("Position-Position computation initialized\n");
  }
  if(strcmp(probes,"ggl_cl")==0){
    like.Ndata=like.Ncl*(tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  }
  if(strcmp(probes,"3x2pt")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  } 
  if(strcmp(probes,"3x2pt_clusterN_clusterWL")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+tomo.cluster_Nbin*Cluster.N200_Nbin+tomo.cgl_Npowerspectra*Cluster.N200_Nbin*Cluster.lbin;
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    like.clusterN=1;
    like.clusterWL=1;
    printf("%d\n",like.Ndata);
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("Cluster Number Counts computation initialized\n");
    printf("Cluster weak lensing computation initialized\n");
  }
   if (strcmp(probes,"LSSxCMB")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    like.kk = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Initializing: gg, gk, gs, kk, ks, ss\n");
   }
   if (strcmp(probes,"gg_gk_gs")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    printf("Initializing: gg, gk, gs\n");
  }
  if (strcmp(probes,"kk_ks_ss")==0) {
    like.Ndata = like.Ncl * (1+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.kk = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Initializing: kk, ks, ss\n");
  }
  printf("Total number of data points like.Ndata=%d\n",like.Ndata);
}


void init_data_inv(char *INV_FILE, char *DATA_FILE)
{
  double init;
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing data vector and covariance\n");
  printf("---------------------------------------\n");

  sprintf(like.INV_FILE,"%s",INV_FILE);
  printf("PATH TO INVCOV: %s\n",like.INV_FILE);
  sprintf(like.DATA_FILE,"%s",DATA_FILE);
  printf("PATH TO DATA: %s\n",like.DATA_FILE);
  init=data_read(0,1);
  init=invcov_read(0,1,1);
}


void init_lens_sample(char *lensphotoz, char *galsample)
{
  if(strcmp(lensphotoz,"none")==0) redshift.clustering_photoz=0;
  if(strcmp(lensphotoz,"voigt")==0) redshift.clustering_photoz=1;
  if(strcmp(lensphotoz,"voigt_out")==0) redshift.clustering_photoz=2;
  if(strcmp(lensphotoz,"gaussian")==0) redshift.clustering_photoz=3;
  if(strcmp(lensphotoz,"multihisto")==0) redshift.clustering_photoz=4;
  
  if ((redshift.clustering_photoz !=0) && (redshift.clustering_photoz !=1) && (redshift.clustering_photoz !=2) && (redshift.clustering_photoz !=3)) 
  {
    printf("init.c: init_lens_sample: redshift.clustering_photoz = %d not set properly!\nEXIT!\n",redshift.clustering_photoz);
    exit(1);
  }
  printf("Lens Sample Redshift Errors set to %s: redshift.clustering_photoz=%d\n",lensphotoz,redshift.clustering_photoz);
  
  if (strcmp(survey.name,"LSST")==0 || strcmp(survey.name,"WFIRST")==0){
    if(strcmp(galsample,"SN10")==0){
      set_galaxies_SN10();
    }
    if(strcmp(galsample,"DES_Y1")==0){
      set_galaxies_DES_Y1();
    }
  }
  //call test_kmax once to initialize look-up tables at reference cosmology
  test_kmax(1000.,1);
}


void init_source_sample(char *sourcephotoz)
{
  if(strcmp(sourcephotoz,"none")==0) redshift.shear_photoz=0;
  if(strcmp(sourcephotoz,"voigt")==0) redshift.shear_photoz=1;
  if(strcmp(sourcephotoz,"voigt_out")==0) redshift.shear_photoz=2;
  if(strcmp(sourcephotoz,"gaussian")==0) redshift.shear_photoz=3;
  if(strcmp(sourcephotoz,"multihisto")==0) {
    printf("redshift.shear_photoz=4 not supported\n"); 
    exit(1);
  }
  if ((redshift.shear_photoz !=0) && (redshift.shear_photoz !=1) && (redshift.shear_photoz !=2) && (redshift.shear_photoz !=3)) 
  {
    printf("init.c: init_source_sample: redshift.shear_photoz = %d not set properly!\nEXIT!\n",redshift.shear_photoz);
    exit(1);
  }

  printf("Source Sample Redshift Errors set to %s: redshift.shear_photoz=%d\n",sourcephotoz,redshift.shear_photoz);
    set_galaxies_source();
}


void set_galaxies_source()
{
  int k,j;
  double frac, zi;
  
  tomo.shear_Npowerspectra=(int) (tomo.shear_Nbin*(tomo.shear_Nbin+1)/2);
  
  zdistr_histo_1(0.1, NULL);
  int zbins =2000;
  double da = (redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.shear_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+zdistr_histo_1(zi, NULL);
  }
  
  tomo.shear_zmin[0] = redshift.shear_zdistrpar_zmin;
  tomo.shear_zmax[tomo.shear_Nbin-1] = redshift.shear_zdistrpar_zmax;
  printf("\n");
  printf("Source Sample - Tomographic Bin limits:\n");
  for(k=0;k<tomo.shear_Nbin-1;k++){
    frac=(k+1.)/(1.*tomo.shear_Nbin)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.shear_zmax[k] = redshift.shear_zdistrpar_zmin+j*da;
    tomo.shear_zmin[k+1] = redshift.shear_zdistrpar_zmin+j*da;
    printf("min=%le max=%le\n",tomo.shear_zmin[k],tomo.shear_zmax[k]);
  }
  printf("min=%le max=%le\n",tomo.shear_zmin[tomo.shear_Nbin-1],tomo.shear_zmax[tomo.shear_Nbin-1]);
  printf("redshift.shear_zdistrpar_zmin=%le max=%le\n",redshift.shear_zdistrpar_zmin,redshift.shear_zdistrpar_zmax);
  free_double_vector(sum,0,zbins);
}

void set_galaxies_SN10()
{
  int k,j,n,i;
  double frac, zi;
  redshift.clustering_zdistrpar_zmin = 0.25;
  redshift.clustering_zdistrpar_zmax = 4.0;
  tomo.clustering_Npowerspectra=tomo.clustering_Nbin;
  tomo.clustering_zmin[0] = redshift.clustering_zdistrpar_zmin;
  tomo.clustering_zmax[tomo.clustering_Nbin-1] = redshift.clustering_zdistrpar_zmax;
  pf_histo(zi, NULL);
  int zbins =2000;
  double da = (redshift.clustering_zdistrpar_zmax-redshift.clustering_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.clustering_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+pf_histo(zi, NULL);
  }
  
  printf("\n");
  printf("Source Sample - Tomographic Bin limits:\n");
  for(k=0;k<tomo.clustering_Nbin-1;k++){
    frac=(k+1.)/(1.*tomo.clustering_Nbin)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.clustering_zmax[k] = redshift.clustering_zdistrpar_zmin+j*da;
    tomo.clustering_zmin[k+1] = redshift.clustering_zdistrpar_zmin+j*da;
    printf("min=%le max=%le\n",tomo.clustering_zmin[k],tomo.clustering_zmax[k]);
  }
  printf("min=%le max=%le\n",tomo.clustering_zmin[tomo.clustering_Nbin-1],tomo.clustering_zmax[tomo.clustering_Nbin-1]);
  printf("redshift.clustering_zdistrpar_zmin=%le max=%le\n",redshift.clustering_zdistrpar_zmin,redshift.clustering_zdistrpar_zmax);
  free_double_vector(sum,0,zbins);
    gbias.b1_function = &b1_per_bin;
  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 1.3+0.1*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n=0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}

void set_galaxies_redmagic()
{

  int i,j,n;
  redshift.clustering_zdistrpar_zmin = 0.2;
  redshift.clustering_zdistrpar_zmax = 1.0;
  //redshift.clustering_histogram_zbins=75;
  
  survey.n_lens = 0.25; //guestimate of stage 4 lens sample with excellent photo-z
  printf("Number density of lens galaxies=%le\n",survey.n_lens);

  tomo.clustering_Nbin        = 4;
  tomo.clustering_Npowerspectra = 4;
  tomo.clustering_zmax[0]      = .4;
  tomo.clustering_zmax[1]      = .6;
  tomo.clustering_zmax[2]      = .8;
  tomo.clustering_zmax[3]      = 1.;
  
  tomo.clustering_zmin[0]      = 0.2;
  tomo.clustering_zmin[1]      = 0.4;
  tomo.clustering_zmin[2]      = 0.6;
  tomo.clustering_zmin[3]      = 0.8;
  printf("\n");
  printf("Lens Sample: LSST redmagic- Tomographic Bin limits:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
  }
  printf("\n");
  printf("Setting Galaxy Bias - Passive Evolution in z-bins:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 1.35+0.15*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
  printf("redshift.clustering_histogram_zbins=%d\n",redshift.clustering_histogram_zbins); 
  
}

void set_galaxies_DES_Y1()
{
  int i,j,n;
  redshift.clustering_zdistrpar_zmin = 0.2;
  redshift.clustering_zdistrpar_zmax = 1.0;
  //redshift.clustering_histogram_zbins=75;
  
  survey.n_lens = 0.15; //guestimate of stage 4 lens sample with excellent photo-z
  printf("Number density of lens galaxies=%le\n",survey.n_lens);

  tomo.clustering_Nbin        = 5;
  tomo.clustering_Npowerspectra = 5;
  tomo.clustering_zmax[0]      = .35;
  tomo.clustering_zmax[1]      = .5;
  tomo.clustering_zmax[2]      = .65;
  tomo.clustering_zmax[3]      = .8;
  tomo.clustering_zmax[4]      = .95;
  
  tomo.clustering_zmin[0]      = 0.2;
  tomo.clustering_zmin[1]      = 0.35;
  tomo.clustering_zmin[2]      = 0.5;
  tomo.clustering_zmin[3]      = 0.65;
  tomo.clustering_zmin[4]      = 0.8;
  printf("\n");
  printf("Lens Sample: LSST redmagic- Tomographic Bin limits:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
  }
  printf("\n");
  printf("Setting Galaxy Bias - Passive Evolution in z-bins:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 1.35+0.15*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
  printf("redshift.clustering_histogram_zbins=%d\n",redshift.clustering_histogram_zbins); 
  
}


// void set_galaxies_source()
// {
//   //take out first lens bin   
//   int i,j,n;
//   tomo.clustering_Nbin        = tomo.shear_Nbin-1;
//   tomo.clustering_Npowerspectra = tomo.clustering_Nbin;
//   for (i = 0; i < tomo.clustering_Nbin; i++){
//     tomo.clustering_zmax[i]      = tomo.shear_zmax[i+1];
//     tomo.clustering_zmin[i]      = tomo.shear_zmin[i+1];
//   }
//   //tomo.clustering_zmin[0]=0.15; //multi-probe NG covs are very likely ill-conditioned if lenses at very low redshift is included 

//   redshift.clustering_zdistrpar_zmin = tomo.clustering_zmin[0];
//   redshift.clustering_zdistrpar_zmax = redshift.shear_zdistrpar_zmax;

//   printf("\n");
//   printf("Lens Sample: Source - Tomographic Bin limits:\n");
//   for (i =0; i < tomo.clustering_Nbin ; i++){
//     printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
//   }
//   gbias.b1_function = &b1_per_bin;
//   for (i =0; i < tomo.clustering_Nbin ; i++){
//     gbias.b[i] = 1.3+0.1*i;
//     printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
//   }
//   n = 0;
//   pf_photoz(0.5,1);
//   printf("redshift.clustering_zdistrpar_zmin %le redshift.clustering_zdistrpar_zmax %le\n",redshift.clustering_zdistrpar_zmin,redshift.clustering_zdistrpar_zmax);
//   for (i = 0; i < tomo.clustering_Nbin; i++){
//     for(j = 0; j<tomo.shear_Nbin;j++){
//       n += test_zoverlap(i,j);
//       printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
//     }
//   }
//   tomo.ggl_Npowerspectra = n;
//   printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
// }


void set_clusters_WFIRST(){
  int i,j;
  //N200->M relationship from Murata et al. (2018)
  nuisance.cluster_Mobs_lgN0 = 3.207; //fiducial: 3.207, flat prior [0.5, 5.0]
  nuisance.cluster_Mobs_alpha = 0.993; //fiducial: 0.993, flat prior [0.0, 2.0]
  nuisance.cluster_Mobs_beta = 0.0; //fiducial: 0.0, flat prior [-1.5, 1.5]
  nuisance.cluster_Mobs_sigma0 = 0.456; //fiducial: 0.456, flat prior [0.0, 1.5]
  nuisance.cluster_Mobs_sigma_qm = 0.0; //fiducial: 0.0, flat prior [-1.5, 1.5]
  nuisance.cluster_Mobs_sigma_qz = 0.0; //fiducial: 0.0, flat prior [-1.5, 1.5]
  //Compliteness parameters are not marginalized, but just fixed to 1.
  nuisance.cluster_completeness[0] = 1.0;
  nuisance.cluster_completeness[1] = 1.0;
  nuisance.cluster_completeness[2] = 1.0;
  nuisance.cluster_completeness[3] = 1.0;

  //no miscentering so far
  nuisance.cluster_centering_f0 = 1.0;
  nuisance.cluster_centering_alpha = 0;
  nuisance.cluster_centering_sigma = 0;
  nuisance.cluster_centering_M_pivot = 1.e+14;
  printf("%e %e %e %e %e %e\n",nuisance.cluster_Mobs_lgN0, nuisance.cluster_Mobs_alpha, nuisance.cluster_Mobs_beta, nuisance.cluster_Mobs_sigma0,  nuisance.cluster_Mobs_sigma_qm, nuisance.cluster_Mobs_sigma_qz);
  tomo.cluster_Nbin = 4; // number of cluster redshift bins
  tomo.cluster_zmin[0] = 0.4;
  tomo.cluster_zmax[0] = 0.6;
  tomo.cluster_zmin[1] = 0.6;
  tomo.cluster_zmax[1] = 0.8;
  tomo.cluster_zmin[2] = 0.8;
  tomo.cluster_zmax[2] = 1.0;
  tomo.cluster_zmin[3] = 1.0;
  tomo.cluster_zmax[3] = 1.2;
  tomo.cgl_Npowerspectra = 0;// number of cluster-lensing tomography combinations
  for (i = 0; i < tomo.cluster_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      tomo.cgl_Npowerspectra += test_zoverlap_c(i,j);
    }
  }
  
  Cluster.N200_min = 25.; //formerly 20
  Cluster.N200_max = 220.;
  Cluster.N200_Nbin = 4;
  strcpy(Cluster.model,"Murata_etal_2018");
  //upper bin boundaries - note that bin boundaries need to be integers!
  int Nlist[4] = {40,80,120,Cluster.N200_max}; //formerly 40,80,120,...
  Cluster.N_min[0] = Cluster.N200_min;
  Cluster.N_max[0] = Nlist[0];
  for (i = 1; i < Cluster.N200_Nbin; i++){
    Cluster.N_min[i] = Nlist[i-1];
    Cluster.N_max[i] = Nlist[i];
  }
 for (i = 0; i < Cluster.N200_Nbin; i++){
    printf ("Richness bin %d: %e - %e, N(z = 0.3) = %e, N(z = 0.7) = %e\n", i,Cluster.N_min[i],Cluster.N_max[i],N_N200(0,i),N_N200(2,i));
  }
  printf("Clusters set to WFIRST\n");
  printf("Clusters cgl_Npowerspectra=%d\n",tomo.cgl_Npowerspectra);
}


void init_IA(char *model,char *lumfct)
{  
  if(strcmp(lumfct,"GAMA")==0) set_LF_GAMA();
  else if(strcmp(lumfct,"DEEP2")==0) set_LF_DEEP2();
  else {
    printf("init.c:init_IA: %s lumfct not defined\n",lumfct);
    printf("USING GAMA LF INSTEAD\n");
    set_LF_GAMA();
  }
  printf("SET LUMINOSITY FUNCTION=%s\n",lumfct);
  
  nuisance.oneplusz0_ia=1.3; 
  //z0=0.3 is arbitrary pivot redshift J11 p18
  nuisance.c1rhocrit_ia=0.0134; 
  // J11 p.8
  
  if(strcmp(model,"none")==0)  like.IA=0;
  else if(strcmp(model,"NLA_HF")==0)  like.IA=1;
  else if(strcmp(model,"lin")==0)  like.IA=2;
  else{
    printf("init.c:init_IA: %s IA model not defined\n",model);
    exit(1);
  }
  printf("SET IA MODEL=%s\n",model);
  set_ia_priors();
  log_like_f_red();
}


void set_wlphotoz_WFIRST_opti()
{
  int i;
  printf("\n");
  printf("Source sample: WFIRST opti photoz uncertainty initialized\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=0.01; 
    printf("nuisance.bias_zphot_shear[%d]=%le\n",i,nuisance.bias_zphot_shear[i]);
    printf("nuisance.sigma_zphot_shear[%d]=%le\n",i,nuisance.sigma_zphot_shear[i]);
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i];
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = 0.002;
    prior.sigma_zphot_shear[i][1]= 0.002;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_shear[i][0],prior.sigma_zphot_shear[i][1]); 
  }
  like.wlphotoz=1;
}

void set_wlphotoz_WFIRST_KL()
{
  int i;
  printf("\n");
  printf("Source sample: WFIRST KL spec-z uncertainty initialized\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=0.002; 
    printf("nuisance.bias_zphot_shear[%d]=%le\n",i,nuisance.bias_zphot_shear[i]);
    printf("nuisance.sigma_zphot_shear[%d]=%le\n",i,nuisance.sigma_zphot_shear[i]);
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i];
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = 0.0004;
    prior.sigma_zphot_shear[i][1]= 0.0004;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_shear[i][0],prior.sigma_zphot_shear[i][1]); 
  }
  like.wlphotoz=1;
}
void set_wlphotoz_DESI2_KL()
{
  int i;
  printf("\n");
  printf("Source sample: DESI-II KL spec-z uncertainty initialized\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=0.0005; 
    printf("nuisance.bias_zphot_shear[%d]=%le\n",i,nuisance.bias_zphot_shear[i]);
    printf("nuisance.sigma_zphot_shear[%d]=%le\n",i,nuisance.sigma_zphot_shear[i]);
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i];
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = 0.0001;
    prior.sigma_zphot_shear[i][1]= 0.0001;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_shear[i][0],prior.sigma_zphot_shear[i][1]); 
  }
  like.wlphotoz=1;
}



void set_wlphotoz_WFIRST_pessi()
{
  int i;
  printf("\n");
  printf("Source sample: WFIRST pessi photoz uncertainty initialized\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=0.05; 
    printf("nuisance.bias_zphot_shear[%d]=%le\n",i,nuisance.bias_zphot_shear[i]);
    printf("nuisance.sigma_zphot_shear[%d]=%le\n",i,nuisance.sigma_zphot_shear[i]);
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i];
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = 0.02;
    prior.sigma_zphot_shear[i][1]= 0.02;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_shear[i][0],prior.sigma_zphot_shear[i][1]); 
  }
  like.wlphotoz=1;
}


void set_clphotoz_WFIRST_opti()
{
  int i;
  printf("\n");
  printf("Lens sample: WFIRST opti photoz uncertainty initialized\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.01; 
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]);
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    // rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.002;
    prior.sigma_zphot_clustering[i][1]= 0.002;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}

void set_clphotoz_WFIRST_KL()
{
  int i;
  printf("\n");
  printf("Lens sample: WFIRST KL spec-z uncertainty initialized\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.002; 
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]);
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    // rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.0004;
    prior.sigma_zphot_clustering[i][1]= 0.0004;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}
void set_clphotoz_DESI2_KL()
{
  int i;
  printf("\n");
  printf("Lens sample: DESI-II KL spec-z uncertainty initialized\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.0005; 
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]);
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    // rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.0001;
    prior.sigma_zphot_clustering[i][1]= 0.0001;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}

void set_clphotoz_WFIRST_pessi()
{
  int i;
  printf("\n");
  printf("Lens sample: WFIRST pessi photoz uncertainty initialized\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.05; 
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]);
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    // rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.02;
    prior.sigma_zphot_clustering[i][1]= 0.02;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}


void set_shear_priors_WFIRST_opti()
{
  int i;
  printf("Setting Gaussian shear calibration Priors stage 4\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = 0.002;
    printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
  }
  like.shearcalib=1;
}

void set_shear_priors_WFIRST_KL()
{
  int i;
  printf("Setting Gaussian shear calibration Priors stage 4\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = 0.0004;
    printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
  }
  like.shearcalib=1;
}
void set_shear_priors_DESI2_KL()
{
  int i;
  printf("Setting Gaussian shear calibration Priors stage 4\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = 0.0004;
    printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
  }
  like.shearcalib=1;
}


void set_shear_priors_WFIRST_pessi()
{
  int i;
  printf("Setting Gaussian shear calibration Priors stage 4\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = 0.01;
    printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
  }
  like.shearcalib=1;
}


void init_HOD_rm(){
  set_HOD_redmagic_priors();
  like.Rmin_bias = 0.1;//use halo+HOD model down to 100 kpc/h
  redm.parameterization = 0; //Zehavi et al. 2011 HOD parameterization
  redm.cg = 1.0;
  redm.fc = 0.2;
  redm.hod[0] = 12.1;
  redm.hod[1] = 0.4;
  redm.hod[2] = 13.65;
  redm.hod[3] = 12.2;
  redm.hod[4] = 1.0;
}

void set_survey_parameters_to_WFIRST_WL()
{
  survey.area   = 2000.;
  survey.n_gal   = 51.;
  survey.sigma_e   = 0.37;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=28.0;
  sprintf(survey.Kcorrect_File,"../zdistris/k+e.dat");
  sprintf(survey.name,"WFIRST_WL");
}

void set_survey_parameters_to_WFIRST_KL()
{
  survey.area   = 2000.;
  survey.n_gal   = 8.;
  survey.sigma_e   = 0.05;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=28.0;
  sprintf(survey.Kcorrect_File,"../zdistris/k+e.dat");
  sprintf(survey.name,"WFIRST_KL");
}

void set_survey_parameters_to_DESI2_KL(char *surveyname)
{
  // example surveyname: DESI2_KL_00
  double source_density[6] = {0.4761, 0.1629, 0.1553, 0.0881, 0.1006, 0.0740};
  double shape_noise_rms[6] = {0.02*1.4142, 0.04*1.4142, 0.06*1.4142, 
                               0.10*1.4142, 0.20*1.4142, 0.30*1.4142};
  char _iSelect[2];
  char _iSN[2];
  strncpy(_iSelect, surveyname+9, 1);
  strncpy(_iSN, surveyname+10, 1);
  int iSelect = atoi(_iSelect);
  int iSN = atoi(_iSN);
  printf("Setting target selection %d and shape noise %d\n", iSelect, iSN);
  survey.area   = 14000.;
  survey.n_gal   = source_density[iSelect];
  survey.sigma_e   = shape_noise_rms[iSN];
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=19.5;
  sprintf(survey.Kcorrect_File,"../zdistris/k+e.dat");
  sprintf(survey.name,surveyname);
}

void set_survey_parameters_to_LSST_Y1()
{ // Table F1 and Appendix C1 in DESC SRD
  survey.area   = 12300.0;
  survey.n_gal   = 11.112;
  survey.sigma_e   = 0.37;  
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=24.5;
  sprintf(survey.name,"LSST_Y1");
}

void set_survey_parameters_to_LSST_Y10()
{ // Table F1 and Appendix C1 in DESC SRD
  survey.area   = 14300.0;
  survey.n_gal   = 27.737;
  survey.sigma_e   = 0.37;  
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=26.0;
  sprintf(survey.name,"LSST_Y10");
}

void set_wlphotoz_LSST_Y1()
{
  int i;
  printf("\n");
  printf("Source sample: LSST Y1 photoz uncertainty initialized\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=0.05; 
    printf("nuisance.bias_zphot_shear[%d]=%le\n",i,nuisance.bias_zphot_shear[i]);
    printf("nuisance.sigma_zphot_shear[%d]=%le\n",i,nuisance.sigma_zphot_shear[i]);
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i];
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = 0.002;
    prior.sigma_zphot_shear[i][1]= 0.006;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_shear[i][0],prior.sigma_zphot_shear[i][1]); 
  }
  like.wlphotoz=1;
}
void set_wlphotoz_LSST_Y10()
{
  int i;
  printf("\n");
  printf("Source sample: LSST Y10 photoz uncertainty initialized\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=0.05; 
    printf("nuisance.bias_zphot_shear[%d]=%le\n",i,nuisance.bias_zphot_shear[i]);
    printf("nuisance.sigma_zphot_shear[%d]=%le\n",i,nuisance.sigma_zphot_shear[i]);
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i];
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = 0.001;
    prior.sigma_zphot_shear[i][1]= 0.003;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_shear[i][0],prior.sigma_zphot_shear[i][1]); 
  }
  like.wlphotoz=1;
}
void set_clphotoz_LSST_Y1()
{
  int i;
  printf("\n");
  printf("Lens sample: LSST Y1 photoz uncertainty initialized\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.03; 
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]);
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    // rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.002;
    prior.sigma_zphot_clustering[i][1]= 0.006;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}
void set_clphotoz_LSST_Y10()
{
  int i;
  printf("\n");
  printf("Lens sample: LSST Y10 photoz uncertainty initialized\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.03; 
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]);
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    // rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.001;
    prior.sigma_zphot_clustering[i][1]= 0.003;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}
void set_shear_priors_LSST_Y1()
{
  int i;
  printf("\nSetting Gaussian shear calibration Priors LSST Y1\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = 0.013;
    printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
  }
  like.shearcalib=1;
}
void set_shear_priors_LSST_Y10()
{
  int i;
  printf("\nSetting Gaussian shear calibration Priors LSST Y10\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = 0.003;
    printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
  }
  like.shearcalib=1;
}

double log_L_PlanckBAOJLA_w0wa()
{
  double log_L = 0.;
  int n_param = 7;
  // https://github.com/CosmoLike/DESC_SRD/blob/master/fisher.py
  double table[7][7] = {
    {8.52889541e+03,9.34266060e-03,2.21513407e+00,1.85008170e+01,6.12185662e+00,1.17720010e+02,5.30670581e+01},
    {9.34266060e-03,3.34475078e+01,3.79135616e-01,-8.02192934e-02,5.88736315e-02,-8.69202865e+01,-7.08029577e+00},
    {2.21513407e+00,3.79135616e-01,2.48739187e+02,-7.46022981e+00,-2.69922217e+00,-1.19693897e+02,4.78085262e+00},
    {1.85008170e+01,-8.02192934e-02,-7.46022981e+00,8.42937127e+01,1.68856256e+01,-4.89063661e+01,-8.77194357e+00},
    {6.12185662e+00,5.88736315e-02,-2.69922217e+00,1.68856256e+01,9.55489400e+00,5.27704214e+00,-1.38597499e+00},
    {1.17720010e+02,-8.69202865e+01,-1.19693897e+02,-4.89063661e+01,5.27704214e+00,7.17777988e+04,3.08491105e+03},
    {5.30670581e+01,-7.08029577e+00,4.78085262e+00,-8.77194357e+00,-1.38597499e+00,3.08491105e+03,5.32751808e+02}
  };
  double param_fid[n_param], param_diff[n_param];
    
  param_diff[0] = cosmology.Omega_m-prior.Omega_m; 
  param_diff[1] = cosmology.sigma_8-prior.sigma_8; 
  param_diff[2] = cosmology.n_spec-prior.n_spec; 
  param_diff[3] = cosmology.w0-prior.w0;
  param_diff[4] = cosmology.wa-prior.wa;
  param_diff[5] = cosmology.omb-prior.omb;
  param_diff[6] = cosmology.h0-prior.h0;
  
  log_L = -0.5*do_matrix_mult_invcov(n_param,table, param_diff);

  return log_L;
}

// DSA configuration
void set_survey_parameters_to_DSA_allsky()
{
  survey.area   = 30000.0;    // DSA all-sky
  survey.n_gal  = 0.0459;     // DSA all-sky
  survey.sigma_e  = 0.05;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor = 1.0/constants.arcmin/constants.arcmin;
  survey.m_lim = 24.5;
  sprintf(survey.name,"DSA_allsky");
}
void set_wlphotoz_DSA_allsky()
{
  int i;
  printf("\n");
  printf("Source sample: DSA all-sky spec-z uncertainty initialized\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=0.0005; 
    printf("nuisance.bias_zphot_shear[%d]=%le\n",i,nuisance.bias_zphot_shear[i]);
    printf("nuisance.sigma_zphot_shear[%d]=%le\n",i,nuisance.sigma_zphot_shear[i]);
    // center of Gaussian priors: copy from DESI-II
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i];
    // rms width of Gaussian priors: randomly pick a number
    prior.bias_zphot_shear[i][1] = 0.0001;
    prior.sigma_zphot_shear[i][1]= 0.0001;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_shear[i][0],prior.sigma_zphot_shear[i][1]); 
  }
  like.wlphotoz=1;
}
void set_clphotoz_DSA_allsky()
{
  // copy from DESI-II
  int i;
  printf("\n");
  printf("Lens sample: DSA all-sky KL spec-z uncertainty initialized\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.0005; 
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]);
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    // rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = 0.0001;
    prior.sigma_zphot_clustering[i][1]= 0.0001;
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}
void set_shear_priors_DSA_allsky() 
{
  // copy from DESI-II
  int i;
  printf("\nSetting Gaussian shear calibration Priors DSA All-sky\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = 0.0004;
    printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
  }
  like.shearcalib=1;
}

// Full SKA configuration
void set_survey_parameters_to_SKA_WL()
{
  // Harrison et al. 2016, Table 1
  survey.area   = 30000.0;
  survey.n_gal  = 10;
  survey.sigma_e  = 0.3;
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor = 1.0/constants.arcmin/constants.arcmin;
  survey.m_lim = 24.5;
  sprintf(survey.name,"SKA_WL");
}
void set_wlphotoz_SKA_WL()
{
  int i;
  printf("\n");
  printf("Source sample: SKA WL spec-z uncertainty initialized\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    // Harrison et al. 2016, Table 1
    if (tomo.shar_zmax[i]>2.0){
      nuisance.bias_zphot_shear[i]=0.0;
      nuisance.sigma_zphot_shear[i]=0.3;         
      // rms width of Gaussian priors: 1/5 of the nuisance
      prior.bias_zphot_shear[i][1] = 0.06;
      prior.sigma_zphot_shear[i][1]= 0.06;
    }
    else {
      nuisance.bias_zphot_shear[i]=0.0;
      nuisance.sigma_zphot_shear[i]=0.03;         
      // rms width of Gaussian priors: 1/5 of the nuisance
      prior.bias_zphot_shear[i][1] = 0.006;
      prior.sigma_zphot_shear[i][1]= 0.006;
    }      
    printf("nuisance.bias_zphot_shear[%d]=%le\n",i,nuisance.bias_zphot_shear[i]);
    printf("nuisance.sigma_zphot_shear[%d]=%le\n",i,nuisance.sigma_zphot_shear[i]);
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i];
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_shear[i][0],prior.bias_zphot_shear[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_shear[i][0],prior.sigma_zphot_shear[i][1]); 
  }
  like.wlphotoz=1;
}
void set_clphotoz_SKA_WL()
{
  int i;
  printf("\n");
  printf("Lens sample: SKA WL spec-z uncertainty initialized\n");
  for (i=0;i<tomo.clustering_Nbin; i++){
    // Harrison et al. 2016, Table 1
    if (tomo.cluster_zmax[i]>2.0){
      nuisance.bias_zphot_clustering[i]=0.0;
      nuisance.sigma_zphot_clustering[i]=0.3; 
      // rms width of Gaussian priors: 1/5 of the nuisance
      prior.bias_zphot_clustering[i][1] = 0.06;
      prior.sigma_zphot_clustering[i][1]= 0.06;
    }
    else {
      nuisance.bias_zphot_clustering[i]=0.0;
      nuisance.sigma_zphot_clustering[i]=0.03; 
      // rms width of Gaussian priors: 1/5 of the nuisance
      prior.bias_zphot_clustering[i][1] = 0.006;
      prior.sigma_zphot_clustering[i][1]= 0.006;
    }
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]);
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    printf("Mean (of mean)=%le, Sigma (of mean)=%le\n",prior.bias_zphot_clustering[i][0],prior.bias_zphot_clustering[i][1]);
    printf("Mean (of sigma)=%le, Sigma (of sigma)=%le\n",prior.sigma_zphot_clustering[i][0],prior.sigma_zphot_clustering[i][1]); 
  }
  like.clphotoz=1;
}

void set_shear_priors_SKA_WL() 
{
  // copy from WFIRST opti
  int i;
  printf("\nSetting Gaussian shear calibration Priors SKA WL\n");
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = 0.002;
    printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[i][0],prior.shear_calibration_m[i][1]);
  }
  like.shearcalib=1;
}
