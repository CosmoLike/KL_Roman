#include "like_fourier.c"

/* Testing the KL cosmic shear functions
*/

int test_RomanPIT_WL(int i_depth, int i_ellmax, char* probe, char* bary_sce)
{
  /* Usage: ./like_fourier [i_depth] [i_ellmax] ["shear_shear"] ["dmo"/"mb2"/...]
    - i_depth: index of survey depth 0-4
    - i_ellmax: index of the ell_max 1000-4000
    - "shear_shear": flag for the likelihood to run
    - "dmo"/"mb2"/...: baryon contamination scenario to used
  */
  clock_t begin, end;
  double time_spent, loglike=0.0, init=0.0;
  int i;

  // 5 sets of survey depth results in different src density [/arcmin2]
  int N_depth = 5;
  char dndz[5][100] = {
     // "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin_norm_neff20",
     // "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin_norm_neff25", 
     // "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin_norm_neff30", 
     // "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin_norm_neff35",
     // "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin_norm_neff40",
     "zdistris/n_eff20.nz",
     "zdistris/n_eff25.nz",
     "zdistris/n_eff30.nz",
     "zdistris/n_eff35.nz",
     "zdistris/n_eff40.nz",
  };
  // 4 sets of ell max
  int N_ellmax = 4;
  double ellmax_list[4] = {1000.,2000.,3000.,4000.};
  int Nell_list[4] = {15, 15, 15, 15};
  double delta_z_src = 0.01;
  double delta_z_lens = 0.01;

  int Ntomo_source = 10;
  
  // Lens galaxies not used, set to random value
  double lens_density = 66.0;
  // Lens galaxies not used, set to random value
  int Ntomo_lens = 10;
  double Rmin_bias = 21.0; // not used 
  // 15 ell bins in Fourier space, from 20 to 3000
  int Nell = Nell_list[i_ellmax];
  double ell_min = 20.0;
  double ell_max = ellmax_list[i_ellmax];
  double ell_max_shear = ellmax_list[i_ellmax];

  char strat[20];
  char invcov_fn[500], dv_fn[500], PCs_fn[500];
  sprintf(invcov_fn, "/xdisk/timeifler/jiachuanxu/RomanPIT/invcov/Roman_WL_%d%d_ssss_invcov_Ncl15_Ntomo10", i_depth, i_ellmax);
  sprintf(dv_fn, "datav/Roman_WL_%d%d_shear_shear_Ntomo10_Ncl15_%s", i_depth, i_ellmax, bary_sce);
  sprintf(PCs_fn, "datav/Roman_WL_%d%d_shear_shear_Ntomo10_Ncl15_9sim.pca", i_depth, i_ellmax);
  sprintf(strat, "Roman_WL_%d0%d", i_depth, i_ellmax);
  /* here, do your time-consuming job */

  init_cosmo_runmode("halofit_split");
  // baryon effects initialization
  // This one is used for applying baryon effects from specific simulation
  // Available choices:
  // "dmo","mb2","illustris","eagle","HzAGN","TNG100","owls_AGN",...
  init_bary(bary_sce);

  init_binning_fourier(Nell, ell_min, ell_max, ell_max_shear, 
    Rmin_bias, Ntomo_source, Ntomo_lens);
  init_priors_IA_bary("photo_opti","shear_opti","none","none",
    false, 3.0, 1.2, 3.8, 2.0, true, 40.0, 10.0, 0.8);
  init_survey(strat);
  init_galaxies(dndz[i_depth], "zdistris/lens_LSSTY1", 
      "gaussian", "gaussian", "SN10");// the last arg is lens sample
  
  #if _WRITE_NZ_TOMO_ == 1
    // write redshift boundary of each tomo bin
    FILE *tomo_zdist;
    char tomo_zdist_fname[500];
    sprintf(tomo_zdist_fname, 
      "zdistris/tomo_zdist_src_%s_Haley_dndz", strat);
    tomo_zdist = fopen(tomo_zdist_fname, "w");
    if(tomo_zdist!=NULL){
      fprintf(tomo_zdist, "# tomo_id\tshear_zmin\tshear_zmax\n");
      for(int i=0; i<tomo.shear_Nbin; i++){
        fprintf(tomo_zdist,"%d\t%.6f\t%.6f\n",i,tomo.shear_zmin[i],tomo.shear_zmax[i]);
      }
    }
    fclose(tomo_zdist);
  #endif

  init_IA("none", "GAMA");
  init_probes(probe);
  
  /* compute fiducial data vector */
  // NOTE: different target selections have different data vectors
  #if _COMPUTE_DATAVECTOR_ == 1
  printf("like.IA = %d\n", like.IA);
  printf("like.baryons = %d\n", like.baryons);
  double fid_sigma8 = 0.831;
  double fid_omegam = 0.3156;
  compute_data_vector("",
    // cosmology+MG: Om, sigma_8, ns, w0, wa, Ob, h0, MG_sigma, MG_mu
    fid_omegam,fid_sigma8,0.9645,-1.0,0.0,0.0491685,0.6727,0.,0.,
    // galaxy bias: b[0-9]
    1.3,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,
    // source galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_src,
    // lens galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_lens,
    // additive shear calibration bias[0-9]
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    // IA: A_ia, beta_ia, eta_ia, eta_ia_highz
    5.92,1.1,-0.47,0.0,
    // lumi function: LF_alpha, LF_P, LF_Q, LF_red_alpha, LF_red_P, LF_red_Q
    0.0,0.0,0.0,0.0,0.0,0.0,
    // clus: mass_obs_norm, mass_obs_slope, mass_z_slope, mass_obs_scatter_norm
    3.207,0.993,0.0,0.456,
    // mass_obs_scatter_mass_slope, mass_obs_scatter_z_slope
    0.0, 0.0,
    // baryon scenario,
    bary_sce,
    // sigma8 split at low-z
    fid_sigma8, 1.0);
  // finite difference dv
  char details_FIM_sigma8[4][50] = {
    "_sigma8--", "_sigma8-", "_sigma8+", "_sigma8++",
  };
  char details_FIM_omegam[4][50] = {
    "_omegam--", "_omegam-", "_omegam+", "_omegam++"
  };
  double finite_diff_scale[4] = {-2., -1., 1., 2.};
  double finite_diff_h = 0.0001;
  for (int ifd=0; ifd<4; ifd++){
    printf("\n\n\n %s \n\n\n", details_FIM_omegam[ifd]);
    compute_data_vector(details_FIM_omegam[ifd],
    fid_omegam+finite_diff_scale[ifd]*finite_diff_h,fid_sigma8,0.9645,-1.0,0.0,0.0491685,0.6727,0.,0.,
    1.3,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_src,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_lens,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    5.92,1.1,-0.47,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,
    3.207,0.993,0.0,0.456,
    0.0, 0.0,
    bary_sce,
    fid_sigma8, 1.0);
    compute_data_vector(details_FIM_sigma8[ifd],
    fid_omegam,fid_sigma8+finite_diff_scale[ifd]*finite_diff_h,0.9645,-1.0,0.0,0.0491685,0.6727,0.,0.,
    1.3,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_src,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_lens,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    5.92,1.1,-0.47,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,
    3.207,0.993,0.0,0.456,
    0.0, 0.0,
    bary_sce,
    fid_sigma8+finite_diff_scale[ifd]*finite_diff_h, 1.0);
  }
  #endif
  /* compute example likelihood evaluation */
  #if _COMPUTE_LIKELIHOOD_ == 1
  init_data_inv_bary(invcov_fn, dv_fn, PCs_fn);
  begin = clock();
  for(int ii=0; ii<1; ii++)
  loglike=log_multi_like(
    // cosmology+MG: Om, S8, ns, w0, wa, Ob, h0, MG_sigma, MG_mu
    0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,
    // galaxy bias: b[0-9]
    1.3,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,
    // source galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_src,
    // lens galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_lens,
    // additive shear calibration bias[0-9]
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    // IA: A_ia, beta_ia, eta_ia, eta_ia_highz
    5.92,1.1,-0.47,0.0,
    // lumi function: LF_alpha, LF_P, LF_Q, LF_red_alpha, LF_red_P, LF_red_Q
    0.0,0.0,0.0,0.0,0.0,0.0,
    // clus: mass_obs_norm, mass_obs_slope, mass_z_slope, mass_obs_scatter_norm
    3.207,0.993,0.0,0.456,
    // mass_obs_scatter_mass_slope, mass_obs_scatter_z_slope,
    0.0,0.0,
    // GRS B1, B2, B3, B4, B5, B6, B7
    1.538026692020565,1.862707210288686,2.213131761595241,2.617023657038295,
    2.975011712138650,3.376705680190931,3.725882076395691,
    // SIGMA P1, P2, P3, P4, P5, P6, P7, Z
    290, 290, 290, 290, 290, 290, 290, 0.001,
    // Pshot, Kstar
    0.0, 0.24,
    // Q1, Q2, Q3
    0.0, 0.0, 0.0,
    // sigma8 split at low-z
    0.831, 1.0);
  printf("loglike = %le\n",loglike);
  // printf("knonlin %le\n",nonlinear_scale_computation(1.0));
  // printf("knonlin %le\n",nonlinear_scale_computation(0.5));
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC/10;
  printf("timespent %.3f\n",time_spent);
  #endif

  // print lensing kernel
  // ====================
  /*
  double da = (1. - limits.a_min)/(Ntable.N_a-1.);
  double a_list[Ntable.N_a], kernel_list[tomo.shear_Nbin][Ntable.N_a];
  for (int i=0;i<Ntable.N_a;i++) a_list[i] = limits.a_min + i*da;
  for (int i=0; i<tomo.shear_Nbin; i++){
    for (int j=0; j<Ntable.N_a; j++){
      double fK = f_K(chi(a_list[j]));
      kernel_list[i][j] = W_kappa(a_list[j],fK,i);
    }
  }
  // write redshift boundary of each tomo bin
  FILE *tomo_kernel;
  char tomo_kernel_fname[500];
  sprintf(tomo_kernel_fname, "zdistris/tomo_kernel_src_%s", strat);
  tomo_kernel = fopen(tomo_kernel_fname, "w");
  if(tomo_kernel!=NULL){
    for (int j=0; j<Ntable.N_a; j++){
      fprintf(tomo_kernel, "%le", a_list[j]);
      for(int i=0; i<tomo.shear_Nbin; i++){
        fprintf(tomo_kernel,"\t%.le",kernel_list[i][j]);
      }
      fprintf(tomo_kernel, "\n");
    }
    fclose(tomo_kernel);
  }
  else{printf("Can not open file %s!\n", tomo_kernel_fname);}
  */
  return 0;
}

int test_DESI2_KL(int i_Selection, int i_SN, char* probe, char* bary_sce)
{
  /* Usage: ./like_fourier [i_selection] [i_SN] ["shear_shear"] ["dmo"/"mb2"/...]
    - i_selection: index of the target selection scenarios, [0, 5]
    - i_SN: index of the shape noise scenarios, [0, 5]
    - "shear_shear": flag for the likelihood to run
    - "dmo"/"mb2"/...: baryon contamination scenario to used
  */
  clock_t begin, end;
  double time_spent, loglike=0.0, init=0.0;
  int i;

  // 6 sets of target selections results in different src density [/arcmin2]
  int N_scenarios_selection = 6;
  char dndz[6][100] = {
    "zdistris/zdistri_DESI2_KL_LS_DR9_sample1_v2",
    "zdistris/zdistri_DESI2_KL_LS_DR9_sample2_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Any_sample1_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Any_sample2_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Bright_sample1_v2",
    "zdistris/zdistri_DESI2_KL_BGS_Bright_sample2_v2",
  };
  double delta_z_src = 0.0005;
  double delta_z_lens = 0.0005;

  int Ntomo_source = 4;
  // 6 sets of shape noise, used to refer to covariance matrix only
  // detailed settings are stored in `set_survey_parameters_to_DESI2_KL()`
  int N_scenarios_shape_noise = 6;
  // Lens galaxies not used, set to random value
  double lens_density = 66.0;
  // Lens galaxies not used, set to random value
  int Ntomo_lens = 10;
  double Rmin_bias = 21.0; // not used 
  // 15 ell bins in Fourier space, from 20 to 3000
  int Nell = 15;
  double ell_min = 20.0;
  double ell_max = 3000.0;
  double ell_max_shear = 3000.0;
  // Now count how many scenarios
  int N_scenarios = N_scenarios_selection * N_scenarios_shape_noise;

  char strat[20];
  char invcov_fn[500], dv_fn[500], PCs_fn[500];
  sprintf(invcov_fn, "/xdisk/timeifler/jiachuanxu/DESI2KL/invcov/DESI2_KL_v2_%d%d_ssss_invcov_Ncl15_Ntomo4", i_Selection, i_SN);
  sprintf(dv_fn, "datav/DESI2_KL_%d%d_shear_shear_Ntomo4_Ncl15_dmo", i_Selection, i_SN);
  sprintf(PCs_fn, "datav/DESI2_KL_%d%d_shear_shear_Ntomo4_Ncl15_9sim.pca", i_Selection, i_SN);

  sprintf(strat, "DESI2_KL_%d%d", i_Selection, i_SN);
  /* here, do your time-consuming job */

  init_cosmo_runmode("halofit_split");
  // baryon effects initialization
  // This one is used for applying baryon effects from specific simulation
  // Available choices:
  // "dmo","mb2","illustris","eagle","HzAGN","TNG100","owls_AGN",...
  init_bary(bary_sce);

  init_binning_fourier(Nell, ell_min, ell_max, ell_max_shear, 
    Rmin_bias, Ntomo_source, Ntomo_lens);
  init_priors_IA_bary("spec_DESI2","shear_KL_DESI2","none","none",
    false, 3.0, 1.2, 3.8, 2.0, true, 40.0, 10.0, 0.8);
  init_survey(strat);
  init_galaxies(dndz[i_Selection], "zdistris/lens_LSSTY1", 
      "gaussian", "gaussian", "SN10");// the last arg is lens sample
  
  #if _WRITE_NZ_TOMO_ == 1
    // write redshift boundary of each tomo bin
    FILE *tomo_zdist;
    char tomo_zdist_fname[500];
    sprintf(tomo_zdist_fname, 
      "zdistris/tomo_zdist_src_%s_%d", strat, i_Selection);
    tomo_zdist = fopen(tomo_zdist_fname, "w");
    if(tomo_zdist!=NULL){
      fprintf(tomo_zdist, "# tomo_id\tshear_zmin\tshear_zmax\n");
      for(int i=0; i<tomo.shear_Nbin; i++){
        fprintf(tomo_zdist,"%d\t%.6f\t%.6f\n",i,tomo.shear_zmin[i],tomo.shear_zmax[i]);
      }
    }
    fclose(tomo_zdist);
  #endif

  init_IA("none", "GAMA");
  init_probes(probe);
  
  /* compute fiducial data vector */
  // NOTE: different target selections have different data vectors
  #if _COMPUTE_DATAVECTOR_ == 1
  printf("like.IA = %d\n", like.IA);
  printf("like.baryons = %d\n", like.baryons);
  compute_data_vector("",
    // cosmology+MG: Om, sigma_8, ns, w0, wa, Ob, h0, MG_sigma, MG_mu
    0.3156,0.831,0.9645,-1.0,0.0,0.0491685,0.6727,0.,0.,
    // galaxy bias: b[0-9]
    1.3,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,
    // source galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_src,
    // lens galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_lens,
    // additive shear calibration bias[0-9]
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    // IA: A_ia, beta_ia, eta_ia, eta_ia_highz
    5.92,1.1,-0.47,0.0,
    // lumi function: LF_alpha, LF_P, LF_Q, LF_red_alpha, LF_red_P, LF_red_Q
    0.0,0.0,0.0,0.0,0.0,0.0,
    // clus: mass_obs_norm, mass_obs_slope, mass_z_slope, mass_obs_scatter_norm
    3.207,0.993,0.0,0.456,
    // mass_obs_scatter_mass_slope, mass_obs_scatter_z_slope
    0.0, 0.0,
    // baryon scenario,
    bary_sce,
    // sigma8 split at low-z
    0.831, 0.15);
  #endif
  /* compute example likelihood evaluation */
  #if _COMPUTE_LIKELIHOOD_ == 1
  init_data_inv_bary(invcov_fn, dv_fn, PCs_fn);
  begin = clock();
  for(int ii=0; ii<10; ii++)
  loglike=log_multi_like(
    // cosmology+MG: Om, S8, ns, w0, wa, Ob, h0, MG_sigma, MG_mu
    0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,
    // galaxy bias: b[0-9]
    1.3,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,
    // source galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_src,
    // lens galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_lens,
    // additive shear calibration bias[0-9]
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    // IA: A_ia, beta_ia, eta_ia, eta_ia_highz
    5.92,1.1,-0.47,0.0,
    // lumi function: LF_alpha, LF_P, LF_Q, LF_red_alpha, LF_red_P, LF_red_Q
    0.0,0.0,0.0,0.0,0.0,0.0,
    // clus: mass_obs_norm, mass_obs_slope, mass_z_slope, mass_obs_scatter_norm
    3.207,0.993,0.0,0.456,
    // mass_obs_scatter_mass_slope, mass_obs_scatter_z_slope,
    0.0,0.0,
    // GRS B1, B2, B3, B4, B5, B6, B7
    1.538026692020565,1.862707210288686,2.213131761595241,2.617023657038295,
    2.975011712138650,3.376705680190931,3.725882076395691,
    // SIGMA P1, P2, P3, P4, P5, P6, P7, Z
    290, 290, 290, 290, 290, 290, 290, 0.001,
    // Pshot, Kstar
    0.0, 0.24,
    // Q1, Q2, Q3
    0.0, 0.0, 0.0,
    // sigma8 split at low-z
    0.831, 0.15);
  printf("%le\n",loglike);
  // printf("knonlin %le\n",nonlinear_scale_computation(1.0));
  // printf("knonlin %le\n",nonlinear_scale_computation(0.5));
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC/10;
  printf("timespent %.3f\n",time_spent);
  #endif

  // print lensing kernel
  // ====================
  /*
  double da = (1. - limits.a_min)/(Ntable.N_a-1.);
  double a_list[Ntable.N_a], kernel_list[tomo.shear_Nbin][Ntable.N_a];
  for (int i=0;i<Ntable.N_a;i++) a_list[i] = limits.a_min + i*da;
  for (int i=0; i<tomo.shear_Nbin; i++){
    for (int j=0; j<Ntable.N_a; j++){
      double fK = f_K(chi(a_list[j]));
      kernel_list[i][j] = W_kappa(a_list[j],fK,i);
    }
  }
  // write redshift boundary of each tomo bin
  FILE *tomo_kernel;
  char tomo_kernel_fname[500];
  sprintf(tomo_kernel_fname, "zdistris/tomo_kernel_src_%s", strat);
  tomo_kernel = fopen(tomo_kernel_fname, "w");
  if(tomo_kernel!=NULL){
    for (int j=0; j<Ntable.N_a; j++){
      fprintf(tomo_kernel, "%le", a_list[j]);
      for(int i=0; i<tomo.shear_Nbin; i++){
        fprintf(tomo_kernel,"\t%.le",kernel_list[i][j]);
      }
      fprintf(tomo_kernel, "\n");
    }
    fclose(tomo_kernel);
  }
  else{printf("Can not open file %s!\n", tomo_kernel_fname);}
  */
  return 0;
}

int test_LSST_WL(int iYear, char* probe, char* bary_sce)
{
  clock_t begin, end;
  double time_spent, loglike=0.0, init=0.0;
  int i;

  int N_scenarios_selection = 2;
  char survey_names[2][100] = {"LSST_Y1", "LSST_Y10"};
  double delta_z_src[2] = {0.05, 0.05};
  double delta_z_lens[2] = {0.03, 0.03};
  char dndz[2][100] = {"zdistris/src_LSSTY1", "zdistris/src_LSSTY10"};
  int Ntomo_source = 10;

  int N_scenarios_shape_noise = 1;
  // Lens galaxies not used, set to random value
  double lens_density = 66.0;
  // Lens galaxies not used, set to random value
  int Ntomo_lens = 10;
  double Rmin_bias = 21.0; // not used 
  // 15 ell bins in Fourier space, from 20 to 3000
  int Nell = 15;
  double ell_min = 20.0;
  double ell_max = 3000.0;
  double ell_max_shear = 3000.0;
  // Now count how many scenarios
  int N_scenarios = N_scenarios_selection * N_scenarios_shape_noise;

  char strat[20];
  sprintf(strat, survey_names[iYear]);
  char invcov_fn[500], dv_fn[500], PCs_fn[500];
  sprintf(invcov_fn, "/xdisk/timeifler/jiachuanxu/DESI2KL/invcov/%s_ssss_invcov_Ncl15_Ntomo10", strat);
  sprintf(dv_fn, "datav/%s_shear_shear_Ntomo10_Ncl15_dmo", strat);
  sprintf(PCs_fn, "datav/%s_shear_shear_Ntomo10_Ncl15_9sim.pca", strat);
  /* here, do your time-consuming job */

  init_cosmo_runmode("halofit_split");
  // baryon effects initialization
  // This one is used for applying baryon effects from specific simulation
  // Available choices:
  // "dmo","mb2","illustris","eagle","HzAGN","TNG100","owls_AGN",...
  init_bary(bary_sce);

  init_binning_fourier(Nell, ell_min, ell_max, ell_max_shear, 
    Rmin_bias, Ntomo_source, Ntomo_lens);
  char _photoz_prior[100];
  char _shearm_prior[100];
  sprintf(_photoz_prior, "photo_%s", strat);
  sprintf(_shearm_prior, "shear_%s", strat);
  init_priors_IA_bary(_photoz_prior, _shearm_prior,"none","none",
    true, 3.0, 1.2, 3.8, 2.0, true, 40.0, 10.0, 0.8);
  init_survey(strat);
  init_galaxies(dndz[iYear], 
      "zdistris/lens_LSSTY1", 
      "gaussian", "gaussian", "SN10");// the last arg is lens sample
  
  #if _WRITE_NZ_TOMO_ == 1
    // write redshift boundary of each tomo bin
    FILE *tomo_zdist;
    char tomo_zdist_fname[500];
    sprintf(tomo_zdist_fname, 
      "zdistris/tomo_zdist_src_%s_%d", strat, iYear);
    tomo_zdist = fopen(tomo_zdist_fname, "w");
    if(tomo_zdist!=NULL){
      fprintf(tomo_zdist, "# tomo_id\tshear_zmin\tshear_zmax\n");
      for(int i=0; i<tomo.shear_Nbin; i++){
        fprintf(tomo_zdist,"%d\t%.6f\t%.6f\n",i,tomo.shear_zmin[i],tomo.shear_zmax[i]);
      }
    }
    fclose(tomo_zdist);
  #endif
  printf("Setting IA...\n");
  init_IA("NLA_HF", "GAMA");
  printf("Probes: %s\n", probe);
  init_probes(probe);
  printf("\n Printing data vector \n");
  /* compute fiducial data vector */
  // NOTE: different target selections have different data vectors
  #if _COMPUTE_DATAVECTOR_ == 1
  printf("like.IA = %d\n", like.IA);
  printf("like.baryons = %d\n", like.baryons);
  compute_data_vector("",
    // cosmology+MG: Om, sigma_8, ns, w0, wa, Ob, h0, MG_sigma, MG_mu
    0.3156,0.831,0.9645,-1.0,0.0,0.0491685,0.6727,0.,0.,
    // galaxy bias: b[0-9]
    1.3,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,
    // source galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_src[iYear],
    // lens galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_lens[iYear],
    // additive shear calibration bias[0-9]
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    // IA: A_ia, beta_ia, eta_ia, eta_ia_highz
    5.92,1.1,-0.47,0.0,
    // lumi function: LF_alpha, LF_P, LF_Q, LF_red_alpha, LF_red_P, LF_red_Q
    0.0,0.0,0.0,0.0,0.0,0.0,
    // clus: mass_obs_norm, mass_obs_slope, mass_z_slope, mass_obs_scatter_norm
    3.207,0.993,0.0,0.456,
    // mass_obs_scatter_mass_slope, mass_obs_scatter_z_slope
    0.0, 0.0,
    // baryon scenario,
    bary_sce,
    // sigma8 split at low-z
    0.831, 0.15);
  #endif
  /* compute example likelihood evaluation */
  #if _COMPUTE_LIKELIHOOD_ == 1
  init_data_inv_bary(invcov_fn, dv_fn, PCs_fn);
  begin = clock();
  for(int ii=0; ii<10; ii++)
  loglike=log_multi_like(
    // cosmology+MG: Om, S8, ns, w0, wa, Ob, h0, MG_sigma, MG_mu
    0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,
    // galaxy bias: b[0-9]
    1.3,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,
    // source galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_src[iYear],
    // lens galaxy photo-z bias[0-9] + std
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,delta_z_lens[iYear],
    // additive shear calibration bias[0-9]
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    // IA: A_ia, beta_ia, eta_ia, eta_ia_highz
    5.92,1.1,-0.47,0.0,
    // lumi function: LF_alpha, LF_P, LF_Q, LF_red_alpha, LF_red_P, LF_red_Q
    0.0,0.0,0.0,0.0,0.0,0.0,
    // clus: mass_obs_norm, mass_obs_slope, mass_z_slope, mass_obs_scatter_norm
    3.207,0.993,0.0,0.456,
    // mass_obs_scatter_mass_slope, mass_obs_scatter_z_slope,
    0.0,0.0,
    // GRS B1, B2, B3, B4, B5, B6, B7
    1.538026692020565,1.862707210288686,2.213131761595241,2.617023657038295,
    2.975011712138650,3.376705680190931,3.725882076395691,
    // SIGMA P1, P2, P3, P4, P5, P6, P7, Z
    290, 290, 290, 290, 290, 290, 290, 0.001,
    // Pshot, Kstar
    0.0, 0.24,
    // Q1, Q2, Q3
    0.0, 0.0, 0.0,
    // sigma8 split at low-z
    0.831, 0.15);
  printf("%le\n",loglike);
  // printf("knonlin %le\n",nonlinear_scale_computation(1.0));
  // printf("knonlin %le\n",nonlinear_scale_computation(0.5));
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC/10;
  printf("timespent %.3f\n",time_spent);
  #endif

  // print lensing kernel
  // ====================
  /*
  double da = (1. - limits.a_min)/(Ntable.N_a-1.);
  double a_list[Ntable.N_a], kernel_list[tomo.shear_Nbin][Ntable.N_a];
  for (int i=0;i<Ntable.N_a;i++) a_list[i] = limits.a_min + i*da;
  for (int i=0; i<tomo.shear_Nbin; i++){
    for (int j=0; j<Ntable.N_a; j++){
      double fK = f_K(chi(a_list[j]));
      kernel_list[i][j] = W_kappa(a_list[j],fK,i);
    }
  }
  // write redshift boundary of each tomo bin
  FILE *tomo_kernel;
  char tomo_kernel_fname[500];
  sprintf(tomo_kernel_fname, "zdistris/tomo_kernel_src_%s", strat);
  tomo_kernel = fopen(tomo_kernel_fname, "w");
  if(tomo_kernel!=NULL){
    for (int j=0; j<Ntable.N_a; j++){
      fprintf(tomo_kernel, "%le", a_list[j]);
      for(int i=0; i<tomo.shear_Nbin; i++){
        fprintf(tomo_kernel,"\t%.le",kernel_list[i][j]);
      }
      fprintf(tomo_kernel, "\n");
    }
    fclose(tomo_kernel);
  }
  else{printf("Can not open file %s!\n", tomo_kernel_fname);}
  */
  return 0;
}

int main(int argc, char** argv)
{	// ./test_cosmic_shear DESI2/LSST/RomanPIT id_sample (id_ellmax) id_baryon 
	char probe[500]="shear_shear";
	char bary_scenarios[12][500]={"dmo", "mb2", "illustris", "eagle", "HzAGN", 
		"TNG100", "cowls_AGN", "cowls_AGN_T8p5", "cowls_AGN_T8p7", 
		"BAHAMAS", "BAHAMAS_T7p6", "BAHAMAS_T8p0"};
	if (strcmp(argv[1], "DESI2")==0){
    int i = atoi(argv[2]);
    int k = atoi(argv[3]);
		assert(i<6);assert(k<12);
		// test DESI2-KL
		test_DESI2_KL(i, 0, probe, bary_scenarios[k]);
	}
	else if (strcmp(argv[1], "LSST")==0){
		int i = atoi(argv[2]);
    int k = atoi(argv[3]);
    assert(i<2);assert(k<12);
		// test LSST cosmic shear
		test_LSST_WL(i, probe, bary_scenarios[k]);
	}
  else if (strcmp(argv[1], "RomanPIT")==0){
    printf("Calculate dv for RomanPIT!\n");
    int i = atoi(argv[2]); // n_eff (dndz)
    int j = atoi(argv[3]); // ell_max
    int k = atoi(argv[4]); // baryons
    // id_sample here means ell_max; n_eff doesn't matter
    assert(i<5);assert(j<4);assert(k<12);
    // test Roman PIT cosmic shear
    test_RomanPIT_WL(i, j, probe, bary_scenarios[k]);
  }
	return 0;
}
