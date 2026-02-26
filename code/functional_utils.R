# Load required libraries
require(dplyr)          # For data manipulation
require(tidyverse)      # For comprehensive data analysis
require(fdapace)        # For functional data analysis (FDA) with PACE method
require(pseudo)         # For pseudo-observations in survival analysis


# ============================================================
#' Function: GetCEScores
#' Calculate Conditional Expectation (CE) Scores for Functional Data
#' 
#' This function computes CE scores for functional data by integrating subject-specific 
#' observations with mean functions, eigenfunctions, and covariance estimates.
#' 
#' @param y List of numeric vectors: Observed functional responses for each subject
#' @param t List of numeric vectors: Time points corresponding to observations in y
#' @param optns List: Options controlling the computation (e.g., verbose mode)
#' @param mu Numeric vector: Estimated mean function evaluated on observation grid
#' @param obsGrid Numeric vector: Regular grid where mean/covariance functions are evaluated
#' @param fittedCov Matrix: Estimated covariance matrix on obsGrid
#' @param lambda Numeric vector: Eigenvalues from functional principal component analysis
#' @param phi Matrix: Eigenfunctions (columns) evaluated on obsGrid
#' @param sigma2 Numeric (optional): Measurement error variance (default = 0)
#' 
#' @return List of lists: Each element contains xiEst (estimated scores), xiVar (score variance), 
#'         and fittedY (fitted values) for a subject
# ============================================================

GetCEScores <- function(y, t, optns, mu, obsGrid, fittedCov, lambda, phi, sigma2) {
  
  # Validate dimension consistency between eigenvalues and eigenfunctions
  if (length(lambda) != ncol(phi))
    stop('No of eigenvalues is not the same as the no of eigenfunctions.')
  
  # Set default measurement error variance if not provided
  if (is.null(sigma2))
    sigma2 <- 0
  
  # Total covariance (systematic + measurement error)
  Sigma_Y <- fittedCov + diag(sigma2, nrow(phi))
  
  # Get subject-specific mean, eigenfunctions, and covariance matrices
  MuPhiSig <- GetMuPhiSig(t, obsGrid, mu, phi, Sigma_Y)
  
  # Apply individual CE score calculation to each subject
  ret <- mapply(function(yVec, muphisig){
    rst = GetIndCEScores(yVec, muphisig$muVec, lambda, muphisig$phiMat, muphisig$Sigma_Yi, verbose=optns$verbose)
    return(rst)
  }, 
  y, MuPhiSig)
  
  return(ret)
}


# ============================================================
#' Function: GetMuPhiSig
#' Prepare Subject-Specific Mean, Eigenfunctions and Covariance
#' 
#' Interpolates population-level mean function, eigenfunctions, and covariance matrix 
#' to match the actual observation time points of each subject.
#' 
#' @param t List of numeric vectors: Observation time points for each subject
#' @param obsGrid Numeric vector: Regular grid of population-level estimates
#' @param mu Numeric vector: Population mean function on obsGrid
#' @param phi Matrix: Population eigenfunctions (columns) on obsGrid
#' @param Sigma_Y Matrix: Population covariance matrix on obsGrid
#' 
#' @return List of lists: Each element contains muVec (interpolated mean), 
#'         phiMat (interpolated eigenfunctions), and Sigma_Yi (interpolated covariance)
# ============================================================

GetMuPhiSig <- function(t, obsGrid, mu, phi, Sigma_Y) {
  
  # Interpolate population-level functions to subject-specific time points
  ret <- lapply(t, function(tvec) {
    if(length(tvec)!=0){
      # Interpolate mean function
      muVec <- approx(obsGrid, mu, tvec)$y
      
      # Interpolate each eigenfunction and arrange as matrix
      phiMat <- matrix(apply(phi, 2, function(phivec){
        return(approx(obsGrid, phivec, tvec)$y)
      }), nrow = length(tvec))
      
      # Create grid for 2D interpolation of covariance matrix
      grid1 <- sapply(tvec, function(x){return(rep(x, length(tvec)))})
      grid2 <- rep(tvec, length(tvec))
      
      # Interpolate covariance matrix
      Sigma_Yi <- matrix(pracma::interp2(obsGrid, obsGrid, Sigma_Y, 
                                         as.numeric(as.vector(grid1)), 
                                         grid2), 
                         length(tvec), length(tvec))
      
      return(list(muVec=muVec, phiMat=phiMat, Sigma_Yi=Sigma_Yi))
    } else {
      # Return empty matrices for subjects with no observations
      return(list(muVec=numeric(0), phiMat=numeric(0), Sigma_Yi=numeric(0)))
    }
  })
  
  return(ret)
}


# ============================================================
#' Function: GetIndCEScores
#' Calculate Individual Conditional Expectation (CE) Scores
#' 
#' Computes CE scores for a single subject, handling special cases with new observation indices
#' and empty observations. Wraps C++ implementation for core computation.
#' 
#' @param yVec Numeric vector: Observed values for a single subject
#' @param muVec Numeric vector: Interpolated mean function values at observation times
#' @param lamVec Numeric vector: Eigenvalues from FPCA
#' @param phiMat Matrix: Interpolated eigenfunctions at observation times
#' @param Sigma_Yi Matrix: Interpolated covariance matrix for the subject
#' @param newyInd Numeric vector (optional): Indices of new observations to handle separately
#' @param verbose Logical: Whether to print warning messages (default = FALSE)
#' 
#' @return List: Contains xiEst (estimated FPC scores), xiVar (variance of scores), 
#'         and fittedY (fitted functional values)
# ============================================================

GetIndCEScores <- function(yVec, muVec, lamVec, phiMat, Sigma_Yi, newyInd=NULL, verbose=FALSE) {
  
  # Handle empty observation case
  if (length(yVec) == 0) {
    if (verbose) {
      warning('Empty observation found, possibly due to truncation')
    }
    return(list(xiEst=matrix(NA, length(lamVec)), 
                xiVar=matrix(NA, length(lamVec), length(lamVec)), 
                fittedY=matrix(NA, 0, 0)))
  }
  
  # Handle special case with new observation indices
  if (!is.null(newyInd)) {    
    if (length(yVec) != 1){ 
      # Extract new observations and update matrices
      newPhi <- phiMat[newyInd, , drop=FALSE]
      newMu <- muVec[newyInd]
      yVec <- yVec[-newyInd]
      muVec <- muVec[-newyInd]
      phiMat <- phiMat[-newyInd, , drop=FALSE]
      Sigma_Yi <- Sigma_Yi[-newyInd, -newyInd, drop=FALSE]  
      
      # Call specialized C++ function for new indices
      return ( GetIndCEScoresCPPnewInd( yVec, muVec, lamVec, phiMat, Sigma_Yi, newPhi, newMu) )
    } else {   
      # Handle edge case with single observation
      Lam <- diag(x=lamVec, nrow = length(lamVec))
      LamPhi <- Lam %*% t(phiMat)
      LamPhiSig <- LamPhi %*% solve(Sigma_Yi)
      xiEst <- LamPhiSig %*% matrix(yVec - muVec, ncol=1)
      xiVar <- Lam - LamPhi %*% t(LamPhiSig) 
      return( list(xiEst=xiEst, xiVar = xiVar, fittedY=NA) )
    }
  } 
  
  # Call core C++ implementation for standard case
  return( GetIndCEScoresCPP( yVec, muVec, lamVec, phiMat, Sigma_Yi) )
}


# ============================================================
#' Function: GetIndCEScoresCPP
#' C++ Wrapper for Individual CE Score Calculation
#' 
#' Low-level function that calls the fdapace package's C++ implementation 
#' for computing individual CE scores (internal use only).
#' 
#' @param yVec Numeric vector: Subject's observed values
#' @param muVec Numeric vector: Subject's mean function values
#' @param lamVec Numeric vector: Eigenvalues
#' @param phiMat Matrix: Subject's eigenfunction matrix
#' @param SigmaYi Matrix: Subject's covariance matrix
#' 
#' @return List: CE score estimates and variance
#' @keywords internal
#' Note: This function is internal to the fdapace package
# ============================================================

GetIndCEScoresCPP <- function(yVec, muVec, lamVec, phiMat, SigmaYi) {
  .Call('_fdapace_GetIndCEScoresCPP', PACKAGE = 'fdapace', yVec, muVec, lamVec, phiMat, SigmaYi)
}


# ============================================================
#' Function: sup_basis_fun_use
#' Generate Survival Supervised Basis Functions  in functional data analysis
#' 
#' Constructs supervised basis functions by integrating functional principal components (FPCs)
#' with survival pseudo-observations, optimizing for prediction performance in survival analysis.
#' 
#' @param K Integer: Number of basis functions
#' @param grid Numeric vector: Time grid for basis function evaluation
#' @param basis Matrix: Original basis functions (columns) evaluated on grid
#' @param tdf Numeric vector: Time step weights for inner product calculation
#' @param theta Numeric: Weight parameter balancing FPC and pseudo-observation contributions (0 ≤ theta ≤ 1)
#' @param fpc List: FPCA results containing xiEst (FPC scores)
#' @param surv_train Data frame: Training survival data with columns (time, event)
#' @param npc Integer: Number of supervised principal components to retain
#' @param plot Logical: Whether to plot cumulative variance explained (default = FALSE)
#' @param end_time Numeric (optional): Maximum follow-up time for pseudo-observation calculation
#' @param covariates Data frame (optional): Additional covariates (not used in current implementation)
#' 
#' @return List: Contains sup_basis (supervised basis functions), 
#'         eigenvalues, and fd_list (basis function coefficients)
# ============================================================

sup_basis_fun_use = function(K, grid, basis, tdf, theta, fpc, surv_train, npc, plot=FALSE, end_time = NULL, covariates = NULL){
  
  # Internal function: Calculate inner product of two basis functions (using trapezoidal rule)
  inprod = function(x, y, tdf){
    sum(x * y * c(0, rep(tdf, (length(x) - 1))))
  }
  
  # Calculate inner product matrix W of basis functions
  W = matrix(0, nrow = K, ncol = K)
  for(i in 1:K){
    for(j in 1:K){
      W[i,j] = inprod(x = basis[,i], y = basis[,j], tdf = tdf)
    }
  }
  
  # Extract principal component scores
  S = t(fpc$xiEst)
  
  # ----------------- Process survival data -----------------
  surv_dat = surv_train
  
  # Calculate pseudo-observations (for survival analysis)
  surv_dat$pseudo <- pseudomean(time = surv_dat$time, 
                                event = surv_dat$event, 
                                tmax = end_time)
  surv_dat$pseudo_use <- log(surv_dat$pseudo)
  
  # Center the response variable
  Y_bar = mean(surv_dat$pseudo_use)
  surv_dat$Y_de_mean = surv_dat$pseudo_use - Y_bar
  train_y = surv_dat$Y_de_mean
  
  # Calculate matrix M (used for supervised basis function calculation)
  maty = matrix(rep(train_y, each = nrow(S)), nrow = nrow(S))
  M = rowSums(maty * W %*% S)
  MM = as.matrix(M) %*% t(as.matrix(M))
  
  # Internal function: Calculate inverse square root of matrix
  sqrM = function (X) {
    EX <- eigen(X)
    VX <- EX$values
    QX <- EX$vectors
    YX <- QX %*% diag(1/sqrt(VX)) %*% t(QX)
    return(YX)
  }
  
  # Calculate U matrix (supervised information matrix)
  U = theta/length(train_y) * W %*% S %*% t(S) %*% W + 
    (1 - theta) * MM/(length(train_y)^2)
  
  # Assume G = W (regularization parameter λ = 0)
  G = W
  halfG_inv = sqrM(G)
  
  # Eigen decomposition
  tM = t(halfG_inv) %*% U %*% halfG_inv
  eigen_res = eigen(tM)
  
  # Optional: Plot variance explained proportion
  if(plot == TRUE){
    print(tM)
    varprp = round(eigen_res$values/sum(eigen_res$values), 4)
    print(varprp)
  }
  
  # Extract first npc supervised basis functions
  fd_list = lapply(1:npc, function(ipc){
    coef_pc = halfG_inv %*% as.matrix(eigen_res$vectors[, ipc])
  })
  
  # Optional: Plot Q statistic
  if(plot == TRUE){
    Q = c()
    Q[1] = (t(fd_list[[1]]) %*% U %*% fd_list[[1]]) / 
      (t(fd_list[[1]]) %*% G %*% fd_list[[1]])
    for(q in 2:length(fd_list)){
      Q[q] = (t(fd_list[[q]]) %*% U %*% fd_list[[q]]) / 
        (t(fd_list[[q]]) %*% G %*% fd_list[[q]])
    }
    
    plot(x = 1:length(fd_list), y = cumsum(Q), cex = 0,
         xlab = 'Number of supervised basis functions',
         ylab = expression('Q('~phi~')'))
    lines(x = 1:length(fd_list), y = cumsum(Q), lty = 2)
  }
  
  # Convert coefficients to basis functions
  sup_basis = NULL
  for(k in 1:npc){
    kth = matrix(fd_list[[k]], nrow = 1) %*% t(basis)
    sup_basis = cbind(sup_basis, t(kth))
  }
  
  return(list(sup_basis = sup_basis, 
              eigenvalues = eigen_res$values[1:npc], 
              coefficients = fd_list))
}


# ============================================================
#' Function: dat.to.mats
#' Convert Long Format Functional Data to Matrix Lists
#' 
#' Transforms longitudinal functional data (long format with subject IDs) 
#' into a list of subject-specific observation vectors (values and time points).
#' 
#' @param data Data frame: Long format data with columns (id, Z1, obstime)
#'        - id: Subject identifier
#'        - Z1: Functional response values
#'        - obstime: Observation time points
#' 
#' @return List: Contains Ly (list of response vectors) and Lt (list of time vectors)
# ============================================================

dat.to.mats = function(data){
  data = as.data.frame(data)
  cens = data
  ids = cens$id
  uids = unique(ids)       # Unique subject IDs
  n_row = length(unique(ids))
  
  dy = list()              # Store subject-wise response values
  dt = list()              # Store subject-wise time points
  ns = c()
  sums = 1
  
  # Extract observations for each subject
  for(i in uids){
    td = cens[cens$id == i,]  # Extract data for current subject
    ty = td$Z1                # Functional response values
    tt = td$obstime           # Corresponding time points
    dy[[sums]] = ty
    dt[[sums]] = tt
    sums = sums + 1
  }
  
  return(list("Ly" = dy, "Lt" = dt))
}


# ============================================================
#' Function: cond.prob.pec.cox
#' Calculate Conditional Survival Probability Using Cox Model
#' 
#' Computes the conditional probability of surviving from Tstart to Tpred
#' using a fitted Cox proportional hazards model (Pseudo-value Estimation for Censored data).
#' 
#' @param model Cox model object: Fitted Cox PH model from survival package
#' @param newdata Data frame: New data for prediction (same covariates as model)
#' @param Tstart Numeric: Start time for conditional probability
#' @param Tpred Numeric: Prediction time (must be ≥ Tstart)
#' 
#' @return Numeric vector: Conditional survival probabilities (risk.Tpred / risk.Tstart)
#' @export
# ============================================================

cond.prob.pec.cox = function(model, newdata, Tstart, Tpred){
  # Survival probability at start time
  risk.Tstart = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tstart)$surv)
  
  # Survival probability at prediction time
  risk.Tpred = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tpred)$surv)
  
  # Conditional probability = P(survive to Tpred | survive to Tstart)
  return(risk.Tpred/risk.Tstart)
}
