#!/bin/bash

COV_DIR=/home/u15/yhhuang/cosmology/CosmoLike/3Dx2D/cov
cd ${COV_DIR}

echo Combine Roman_KL_llll_cov_Ncl20_Ntomo8
cat Roman_KL_llll_cov_Ncl20_Ntomo8_* > Roman_KL_llll_cov_Ncl20_Ntomo8
rm Roman_KL_llll_cov_Ncl20_Ntomo8_*

echo Combine Roman_KL_llls_cov_Ncl20_Ntomo8
cat Roman_KL_llls_cov_Ncl20_Ntomo8_* > Roman_KL_llls_cov_Ncl20_Ntomo8
rm Roman_KL_llls_cov_Ncl20_Ntomo8_*

echo Combine Roman_KL_llss_cov_Ncl20_Ntomo8
cat Roman_KL_llss_cov_Ncl20_Ntomo8_* > Roman_KL_llss_cov_Ncl20_Ntomo8
rm Roman_KL_llss_cov_Ncl20_Ntomo8_*

echo Combine Roman_KL_lsss_cov_Ncl20_Ntomo8
cat Roman_KL_lsss_cov_Ncl20_Ntomo8_* > Roman_KL_lsss_cov_Ncl20_Ntomo8
rm Roman_KL_lsss_cov_Ncl20_Ntomo8_*

echo Combine Roman_KL_lsls_cov_Ncl20_Ntomo8
cat Roman_KL_lsls_cov_Ncl20_Ntomo8_* > Roman_KL_lsls_cov_Ncl20_Ntomo8
rm Roman_KL_lsls_cov_Ncl20_Ntomo8_*

echo Combine Roman_KL_ssss_cov_Ncl20_Ntomo8
cat Roman_KL_ssss_cov_Ncl20_Ntomo8_* > Roman_KL_ssss_cov_Ncl20_Ntomo8
rm Roman_KL_ssss_cov_Ncl20_Ntomo8_*

echo Done
date