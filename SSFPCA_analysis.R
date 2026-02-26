# ==============================================================================
# SSFPCA (Supervised Sparse Functional Principal Component Analysis)
# Survival Prediction with Functional Data - Simulation Analysis
# 
# This script implements:
# 1. 5-fold cross-validation to tune the theta parameter for SSFPCA
# 2. Survival prediction using SSFPCA vs standard FPCA
# 3. Evaluation metrics: AUC (Time-Dependent ROC) and Brier Score (BS)
# 
# Key parameters:
# - t_star: Time point for conditional prediction (baseline time)
# - t_delta: Prediction horizon (time after t_star)
# - t_pred: Target prediction time (t_star + t_delta)
# - theta: Tuning parameter balancing FPC and survival information
# ==============================================================================

# --------------------------
# 1. Load Data Files
# --------------------------
# Load survival and longitudinal data (training/testing sets)
load("data/surv_train.RData")  # Training survival data (id, time, event)
load("data/surv_test.RData")   # Testing survival data
load("data/long_train.RData")  # Training longitudinal functional data (id, obstime, Z1)
load("data/long_test.RData")   # Testing longitudinal functional data


# --------------------------
# 2. Load Custom Functions
# --------------------------
source("R/functional_utils.R")  # Load SSFPCA core functions


# --------------------------
# 3. Define Core Parameters
# --------------------------
L = 14                          # Maximum follow-up time
t_delta = 5                    # Prediction horizon (time after t_star)
t_star = 5                     # Baseline time for conditional prediction
t_pred = t_star + t_delta       # Target prediction time (survival to t_pred given survival to t_star)


# --------------------------
# 4. Data Preprocessing
# --------------------------
# Convert longitudinal data to subject-specific lists (Ly: response values, Lt: time points)
raw_obs = dat.to.mats(long_train)

# Filter test data to subjects surviving past t_star (conditional prediction population)
m = long_test[long_test$time > t_star,]          # Test subjects with follow-up beyond t_star
mid = unique(m$id)                               # Unique IDs of eligible test subjects

# Restrict test longitudinal data to observations BEFORE t_star (prediction window)
testd = long_test[long_test$id %in% mid,]
testd = testd[testd$obstime <= t_star,]          # Only use observations up to baseline time t_star

# Filter test survival data to match eligible subjects
surv_test = surv_test[surv_test$id %in% mid,]

# Convert filtered test data to subject-specific lists
test_r = dat.to.mats(testd)

# Re-define training data (redundant assignment, kept for consistency with original code)
tdat = long_train  
tdat.id = surv_train
raw_obs = dat.to.mats(tdat)

# Re-filter test data (redundant step, kept for consistency)
m = long_test[long_test$time > t_star,]
mid = unique(m$id)
testd = long_test[long_test$id %in% mid,]
testd = testd[testd$obstime <= t_star,]
surv_test = surv_test[surv_test$id %in% mid,]
test_r = dat.to.mats(testd)

# Re-define prediction time (redundant assignment)
t_pred = t_star + t_delta

# --------------------------
# 5. Hyperparameter Tuning (theta)
# --------------------------
# 5-fold cross-validation to optimize theta parameter
nfold = 5                                         # Number of cross-validation folds
set.seed(1234)                                    # Set seed for reproducibility
id.bank = surv_train$id                           # All training subject IDs
ntest.id = NULL                                   # List to store test IDs for each fold

# Split training IDs into 5 folds (stratified random sampling)
for(i in 1:(nfold-1)){
  # Sample 1/nfold of IDs for test fold (without replacement)
  ntest.id[[i]] = sample(id.bank, size = (1/nfold)*length(surv_train$id), replace = F)
  # Remove sampled IDs from remaining pool
  id.bank = id.bank[-which(id.bank %in% unlist(ntest.id))]
}
# Assign remaining IDs to the last fold
ntest.id[[nfold]] = id.bank 

# Define theta tuning grid (logarithmic spacing + linear steps)
t_tune = c(0.0001, 0.001, 0.01, seq(0.1, 0.9, by=0.1))
npc = 1                                           # Number of supervised principal components to retain
nbasis = 5                                        # Number of basis functions for FPCA

# Initialize results table to store CV metrics
pred_tab = NULL

# --------------------------
# 6. Cross-Validation Loop
# --------------------------
for(d in 1:nfold){
  # Initialize metrics for current fold
  AUC_tune = BS_tune = NULL
  
  # Test each theta value in tuning grid
  for(z in 1:length(t_tune)){
    # --------------------------
    # 6.1 Split Data for Current Fold
    # --------------------------
    # Training data: exclude current fold IDs
    long_train1 = long_train[!(long_train$id %in% ntest.id[[d]]),]
    surv_train1 = surv_train[!(surv_train$id %in% ntest.id[[d]]),]
    
    # Validation data: include only current fold IDs
    long_test1 = long_train[long_train$id %in% ntest.id[[d]],]
    surv_test1 = surv_train[surv_train$id %in% ntest.id[[d]],]
    
    # Re-define training data for current fold
    tdat1 = long_train1 
    tdat.id1 = surv_train1
    raw_obs1 = dat.to.mats(tdat1)
    
    # Filter validation data to subjects surviving past t_star
    m1 = long_test1[long_test1$time >= t_star,]
    mid1 = unique(m1$id)
    
    # Restrict validation longitudinal data to observations BEFORE t_star
    testd1 = long_test1[long_test1$id %in% mid1,]
    testd1 = testd1[testd1$obstime <= (t_star),]  # Predict using only pre-t_star observations
    surv_test1 = surv_test1[surv_test1$id %in% mid1,]
    test_r1 = dat.to.mats(testd1)
    
    # --------------------------
    # 6.2 Run Standard FPCA
    # --------------------------
    D = nbasis
    # Perform Sparse FPCA on training data
    fpc = FPCA(Ly = raw_obs1[[1]],                # List of subject-specific response vectors
               Lt = raw_obs1[[2]],                # List of subject-specific time points
               list(maxK = D, dataType = "Sparse"))# FPCA options: max basis functions, sparse data type
    
    # Calculate time step for numerical integration
    tdf = fpc$workGrid[2] - fpc$workGrid[1]
    basis = fpc$phi                               # Estimated eigenfunctions (basis functions)
    noft = nrow(basis)                            # Number of time grid points
    
    # --------------------------
    # 6.3 Generate Supervised Basis Functions (SSFPCA)
    # --------------------------
    sup_fun = sup_basis_fun_use(K = ncol(fpc$phi),        # Number of basis functions
                                grid = fpc$workGrid,      # Time grid for basis evaluation
                                basis = fpc$phi[,c(1:ncol(fpc$phi))],  # Original FPCA basis
                                tdf = diff(fpc$workGrid)[1],           # Time step size
                                theta = t_tune[z],                     # Current theta value
                                fpc = fpc,                             # FPCA results object
                                surv_train = surv_train1,              # Survival data for supervision
                                npc = npc,                             # Number of supervised PCs
                                end_time = L)                          # Maximum follow-up time
    
    # Extract supervised basis functions and eigenvalues
    o_tfs = sup_fun[[1]]                          # Supervised basis functions matrix
    eigen_val = as.matrix(sup_fun[[2]], nrow = 1) # Eigenvalues for supervised components
    
    # --------------------------
    # 6.4 Construct Covariance Matrix for SSFPCA
    # --------------------------
    # Initialize covariance matrix (same dimension as FPCA work grid)
    covZZ = matrix(0, nrow = length(fpc$workGrid), ncol = length(fpc$workGrid))
    
    # Calculate covariance matrix using supervised basis and eigenvalues
    for(rr in 1:length(fpc$workGrid)){
      for(cc in 1:length(fpc$workGrid)){
        rr_mat = matrix(o_tfs[rr,]*as.numeric(eigen_val), nrow = 1)
        cc_mat = matrix(o_tfs[cc,], ncol = 1)
        covZZ[rr, cc] = as.numeric(rr_mat %*% cc_mat)
      }
    }
    
    # --------------------------
    # 6.5 Calculate CE Scores (Training Data)
    # --------------------------
    # Compute Conditional Expectation scores using SSFPCA estimates
    scoresObj <- GetCEScores(y = raw_obs1[[1]],            # Training response values
                             t = raw_obs1[[2]],            # Training time points
                             mu = fpc$mu,                  # FPCA mean function
                             obsGrid = fpc$workGrid,       # FPCA time grid
                             fittedCov = covZZ,            # SSFPCA covariance matrix
                             lambda = eigen_val,           # SSFPCA eigenvalues
                             phi=o_tfs,                    # SSFPCA eigenfunctions
                             sigma2=fpc$sigma2)            # Measurement error variance
    
    # Extract CE scores into matrix format (subjects x components)
    o_gam = NULL
    for(ii in 1:length(scoresObj[1,])){
      o_gam = rbind(o_gam,as.vector(unlist(scoresObj[1,ii])))
    }
    
    # Assign column names to score matrix
    score_names = c()
    for(q in 1:(ncol(o_gam))){
      tname = paste("score", as.character(q), sep = "")
      score_names = c(score_names, tname)
    }
    
    # Merge scores with survival data for Cox model fitting
    tdat.id1 = cbind(tdat.id1, o_gam)
    tdat.id1 <- na.omit(tdat.id1)                     # Remove subjects with missing scores
    colnames(tdat.id1)[(ncol(tdat.id1) - ncol(o_gam) + 1) : ncol(tdat.id1)] = score_names
    
    # --------------------------
    # 6.6 Fit Cox Proportional Hazards Model
    # --------------------------
    # Formula: Survival ~ clinical covariates + SSFPCA scores
    fmla1 = as.formula(paste("Surv(time,event) ~ drug + age +sex +", paste(score_names, collapse= "+")))
    
    # Fit Cox model with SSFPCA scores
    fitted_obj_SSFPCA = coxph(fmla1, data = tdat.id1, x = T, y = T)
    
    # --------------------------
    # 6.7 Calculate CE Scores (Validation Data)
    # --------------------------
    # Compute CE scores for validation data using training SSFPCA model
    scoresObj_test <- GetCEScores(y = test_r1[[1]],        # Validation response values
                                  t = test_r1[[2]],        # Validation time points
                                  mu = fpc$mu,             # Training mean function
                                  obsGrid = fpc$workGrid,  # Training time grid
                                  fittedCov = covZZ,       # Training SSFPCA covariance
                                  lambda = eigen_val,      # Training SSFPCA eigenvalues
                                  phi=o_tfs,               # Training SSFPCA eigenfunctions
                                  sigma2=fpc$sigma2)       # Training measurement error
    
    # Extract validation scores into matrix format
    test_gam = NULL
    for(ii in 1:length(scoresObj_test[1,])){
      test_gam = rbind(test_gam,as.vector(unlist(scoresObj_test[1,ii])))
    }
    
    # --------------------------
    # 6.8 Prepare Validation Data for Prediction
    # --------------------------
    # Re-assign score column names
    score_names = c()
    for(q in 1:(ncol(test_gam))){
      tname = paste("score", as.character(q), sep = "")
      score_names = c(score_names, tname)
    }
    
    # Merge scores with validation survival data
    surv_test1 = cbind(surv_test1, test_gam)
    surv_test1 <- na.omit(surv_test1)
    colnames(surv_test1)[(ncol(surv_test1) - ncol(o_gam) + 1) : ncol(surv_test1)] = score_names
    
    # --------------------------
    # 6.9 Calculate Prediction Metrics
    # --------------------------
    # 6.9.1 Conditional Survival Probability (P(survive to t_pred | survive to t_star))
    tmp.surv.data = surv_test1
    DP.prob = cond.prob.pec.cox(fitted_obj_SSFPCA, tmp.surv.data, t_star, t_pred)
    
    # 6.9.2 Kaplan-Meier estimator for baseline survival
    surv.train =  tdat.id1
    km = survfit(Surv(time, event)~1, data=surv.train)  # KM curve for training data
    survest = stepfun(km$time, c(1, km$surv))           # Step function for KM survival estimates
    
    # 6.9.3 Brier Score (BS) calculation
    surv.test = surv_test1
    N_vali = nrow(surv.test)                            # Number of validation subjects
    tp = t_pred                                         # Prediction time
    D = rep(0, N_vali)                                  # Event indicator (1 = event by t_pred)
    D[surv.test$time<=tp&surv.test$event==1] = 1        # Mark subjects with event before/at t_pred
    
    pi = 1-DP.prob                                      # Predicted event probability (1 - survival probability)
    
    # Weight calculation for weighted Brier score
    km_pts = survest(surv.test$time)/survest(t_star)    # Conditional KM survival
    W2 <- D/km_pts                                      # Weight for events
    W1 <- as.numeric(surv.test$time>tp)/(survest(tp)/survest(t_star))  # Weight for censored subjects
    W <- W1 + W2                                        # Combined weight
    
    # Calculate individual Brier score points
    BS_pts <- W * (D - pi)^2
    BS.SSFPCA = sum(na.omit(BS_pts)) / N_vali           # Average Brier score for current theta
    
    # Store Brier score for current theta
    BS_tune = c(BS_tune, BS.SSFPCA)
    
    # 6.9.4 Time-Dependent AUC calculation
    roc = tdROC::tdROC(
      X = 1 - cond.prob.pec.cox(model = fitted_obj_SSFPCA,  # Predictor: 1 - survival probability
                                newdata = surv_test1, 
                                Tstart = t_star, 
                                Tpred = t_pred),
      Y= surv_test1$time,                                   # Survival time
      delta = surv_test1$event,                             # Event indicator
      tau = t_pred,                                         # Prediction time
      span = 0.05,                                          # Smoothing parameter
      nboot = 0,                                            # No bootstrap for CI
      alpha = 0.05,                                         # Significance level
      n.grid = 1000,                                        # Grid points for ROC
      cut.off = 0.5)                                        # Threshold for classification
    
    # Extract AUC value
    auc.SSFPCA = roc$AUC$value
    AUC_tune = c(AUC_tune, auc.SSFPCA)                     # Store AUC for current theta
  }
  
  # Store fold results (fold ID, theta, AUC, Brier Score)
  pred_tab = rbind(pred_tab, cbind(rep(d, length(t_tune)), t_tune, AUC_tune, BS_tune))
}

# --------------------------
# 7. Post-Processing CV Results
# --------------------------
# Convert results to data frame and rename columns
pred_tab = as.data.frame(pred_tab)
colnames(pred_tab)[1] = 'fold'

# Plot tuning results (Note: Fixed column name typos from original code)
plot(pred_tab$t_tune, pred_tab$BS_tune, main = "Brier Score vs Theta", 
     xlab = "Theta", ylab = "Brier Score", type = "b")
plot(pred_tab$t_tune, pred_tab$AUC_tune, main = "AUC vs Theta", 
     xlab = "Theta", ylab = "AUC", type = "b")

# --------------------------
# 8. Select Optimal Theta
# --------------------------
pick_me = NULL
# For each fold, select theta with minimum Brier Score
for(d in 1:nfold){
  sub = pred_tab[pred_tab$fold == d,]
  # Skip if all Brier scores are identical
  if(!all(sub$BS_tune == sub$BS_tune[1])){
    # Select theta with minimum Brier Score
    pick_me = c(pick_me, t_tune[which(sub$BS_tune == min(sub$BS_tune))])
  }
}

# Count frequency of selected thetas across folds
pick_me = as.data.frame(table(pick_me))
colnames(pick_me) = c('tune', 'freq')
pick_me$tune = as.numeric(levels(pick_me$tune))

# Select the most frequent optimal theta (majority vote)
tune_picked_SSFPCA = pick_me[which(pick_me$freq == max(pick_me$freq)),]$tune[1]

# --------------------------
# 9. Final Model with Optimal Theta
# --------------------------
# 9.1 Run Standard FPCA on full training data
D = nbasis
fpc = FPCA(Ly = raw_obs[[1]], 
           Lt = raw_obs[[2]], 
           list(maxK = D, dataType = "Sparse"))
tdf = fpc$workGrid[2] - fpc$workGrid[1]
basis = fpc$phi
noft = nrow(basis)

# 9.2 Generate SSFPCA with optimal theta
sup_fun_SSFPCA = sup_basis_fun_use(K = ncol(fpc$phi),
                                   grid = fpc$workGrid,
                                   basis = fpc$phi[,c(1:ncol(fpc$phi))],
                                   tdf = diff(fpc$workGrid)[1],
                                   theta = tune_picked_SSFPCA,  # Optimal theta from CV
                                   fpc = fpc,
                                   surv_train = tdat.id,
                                   npc = npc,
                                   end_time = L)

# Extract SSFPCA results
o_tfs_SSFPCA = sup_fun_SSFPCA[[1]]
eigen_val_SSFPCA = as.matrix(sup_fun_SSFPCA[[2]], nrow = 1)

# Initialize final metrics storage
auc.FPCA = auc.SSFPCA = NULL
BS.FPCA = BS.SSFPCA = NULL

# --------------------------
# 10. Final Prediction with SSFPCA
# --------------------------
# 10.1 Reconstruct covariance matrix with optimal theta
covZZ = matrix(0, nrow = length(fpc$workGrid), ncol = length(fpc$workGrid))
for(rr in 1:length(fpc$workGrid)){
  for(cc in 1:length(fpc$workGrid)){
    rr_mat = matrix(o_tfs_SSFPCA[rr,]*as.numeric(eigen_val_SSFPCA), nrow = 1)
    cc_mat = matrix(o_tfs_SSFPCA[cc,], ncol = 1)
    covZZ[rr, cc] = as.numeric(rr_mat %*% cc_mat)
  }
}

# 10.2 Calculate CE scores for full training data
scoresObj <- GetCEScores(y = raw_obs[[1]], 
                         t = raw_obs[[2]], 
                         mu = fpc$mu, 
                         obsGrid = fpc$workGrid, 
                         fittedCov = covZZ, 
                         lambda = matrix(eigen_val_SSFPCA, ncol = 1), 
                         phi= matrix(o_tfs_SSFPCA, ncol = npc),
                         sigma2=fpc$sigma2)

# Extract training scores
o_gam = NULL
for(ii in 1:length(scoresObj[1,])){
  o_gam = rbind(o_gam,as.vector(unlist(scoresObj[1,ii])))
}

# Assign score column names
score_names = c()
for(q in 1:npc){
  tname = paste("score", as.character(q), sep = "")
  score_names = c(score_names, tname)
}

# Merge scores with full training survival data
tdat.id_SSFPCA = cbind(surv_train, o_gam[,1:npc])
tdat.id_SSFPCA <- na.omit(tdat.id_SSFPCA)
colnames(tdat.id_SSFPCA)[(ncol(tdat.id_SSFPCA) - npc + 1) : ncol(tdat.id_SSFPCA)] = score_names

# 10.3 Fit final Cox model with optimal SSFPCA
fmla = as.formula(paste("Surv(time,event) ~ drug + age +sex +", paste(score_names, collapse= "+")))
fitted_obj_SSFPCA = coxph(fmla, data = tdat.id_SSFPCA, x = T, y = T)

# 10.4 Calculate CE scores for test data
scoresObj_test <- GetCEScores(y = test_r[[1]], 
                              t = test_r[[2]], 
                              mu = fpc$mu, 
                              obsGrid = fpc$workGrid, 
                              fittedCov = covZZ, 
                              lambda = matrix(eigen_val_SSFPCA, ncol = 1), 
                              phi= matrix(o_tfs_SSFPCA, ncol = npc),
                              sigma2=fpc$sigma2)

# Extract test scores (only keep top npc components)
test_gam = NULL
for(ii in 1:length(scoresObj_test[1,])){
  test_gam = rbind(test_gam,as.vector(unlist(scoresObj_test[1,ii])))
}
test_gam = test_gam[,1:npc]

# Assign score column names for test data
score_names = c()
for(q in 1:npc){
  tname = paste("score", as.character(q), sep = "")
  score_names = c(score_names, tname)
}

# Merge scores with test survival data
surv_test_SSFPCA = cbind(surv_test, test_gam)
surv_test_SSFPCA <- na.omit(surv_test_SSFPCA)
colnames(surv_test_SSFPCA)[(ncol(surv_test_SSFPCA) - npc + 1) : ncol(surv_test_SSFPCA)] = score_names

# 10.5 Final AUC calculation
roc = tdROC::tdROC(
  X = 1 - cond.prob.pec.cox(model = fitted_obj_SSFPCA, 
                            newdata = surv_test_SSFPCA, 
                            Tstart = t_star, 
                            Tpred = t_pred),
  Y= surv_test_SSFPCA$time,
  delta = surv_test_SSFPCA$event,
  tau = t_pred, 
  span = 0.05,
  nboot = 0, 
  alpha = 0.05,
  n.grid = 1000, 
  cut.off = 0.5)

auc.SSFPCA = roc$AUC$value

# 10.6 Final Brier Score calculation
tmp.surv.data = surv_test_SSFPCA
DP.prob = cond.prob.pec.cox(fitted_obj_SSFPCA, tmp.surv.data, t_star, t_pred)

# KM estimator for baseline survival
surv.train =  tdat.id_SSFPCA
km = survfit(Surv(time, event)~1, data=surv.train )
survest = stepfun(km$time, c(1, km$surv))

# Brier Score components
surv.test = surv_test_SSFPCA
N_vali = nrow(surv.test)
tp = t_pred
D = rep(0, N_vali)
D[surv.test$time<=tp&surv.test$event==1] = 1
pi = 1-DP.prob

# Weight calculation
km_pts = survest(surv.test$time)/survest(t_star)
W2 <- D/km_pts
W1 <- as.numeric(surv.test$time>tp)/(survest(tp)/survest(t_star))
W <- W1 + W2

# Final Brier Score
BS_pts <- W * (D - pi)^2
BS.SSFPCA = sum(na.omit(BS_pts)) / N_vali