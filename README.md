# SSFPCA: Survival-Supervised Functional Principal Component Analysis

Overview
Traditional FPCA extracts features that capture the maximum variance in longitudinal trajectories, but these features may not be optimal for predicting survival outcomes. 
SSFPCA incorporates survival information during the dimension reduction process to extract functional features that are most relevant to survival outcomes, thereby enhancing prediction accuracy. 
This repository provides a complete implementation of SSFPCA along with an application to the Mayo Clinic Primary Biliary Cirrhosis Data, demonstrating how longitudinal serum bilirubin measurements can be used to improve survival prediction.

File Descriptions
1. SSFPCA_analysis.R
Full analysis pipeline for survival prediction based on SSFPCA, main analysis script that demonstrates:
Loading and preparing the PBC2 dataset
Cross-validation for supervision intensity parameter (theta) tuning
Training SSFPCA model with optimal theta and calculate CE scores
Fit Cox model and evaluating model performance on test data

3. functional_utils.R
Core utility functions for SSFPCA implementation—designed to extract survival-relevant functional components:
GetCEScores(): Calculate Conditional Expectation (CE) scores for functional longitudinal data (core step to quantify subject-specific functional component contributions)
GetMuPhiSig(): Interpolate population-level mean/eigenfunctions to subject-specific observation times (align functional estimates with individual measurement timelines)
GetIndCEScores(): Compute individual CE scores (wraps C++ implementation from fdapace for efficiency)
sup_basis_fun_use(): Generate survival-supervised basis functions by integrating FPCA results with survival pseudo-observations (the core of SSFPCA—links functional features to survival outcomes)
dat.to.mats(): Convert long-format longitudinal data to subject-specific lists (values + time points) for FPCA compatibility
cond.prob.pec.cox(): Calculate conditional survival probability using a fitted Cox proportional hazards model (predicts survival to t_pred given survival to t_star)

Data Description
The analysis uses th pbc2 dataset from the JM R package. The used dataset includes:
Longitudinal data: Repeated measurements of serum bilirubin, a key biomarker for liver function
Survival data: Time to death or liver transplantation with censoring information
Baseline covariates: drug (D-penicillamine vs. placebo), age, and sex

The data has been preprocessed and split into training and testing sets:
surv_train.RData: Survival training data (212 patients)
surv_test.RData: Survival testing data (100 patients)
long_train.RData: Longitudinal training data (serum bilirubin measurements)
long_test.RData: Longitudinal testing data (serum bilirubin measurements)



