#' Fit Inverse Probability Weighting using Covariate Balancing Methods: IPT, Graham et al. (2012), or CBPS, Imai and Ratkovic (2014)
#'
#' Fitting function called by \link[teffects2]{teffect} when the method is
#' set to \code{'IPW'} and \code{'cov.balancing'} is set to \code{'IPT'} or \code{'CBPS'}.
#' 
#' Inverse Probability Weighting fits a model for the exposure model to obtain observation weights using either:
#' \code{'IPT'}, Inverse Probability Tilting according to Graham et al. (2012) which solves a method-of-moments problem to derive estimates
#' for the exposure (propensity score model).
#' \code{'CBPS'}, Covariate balancing according to Imai and Ratkovic (2014) which effectively weighs the general maximum-likelihood estimation
#' moment conditions by their estimated inverse conditional variance to estimate the exposure (propensity score model).
#' Predictions from this model are used to weight the observations and derive
#' potential outcome means.
#'
#' @param X numeric vector of (0,1) specifiying the treatment variable.
#' @param Y numeric vector sepcifying the outcome variable
#' @param Zx design matrix for the exposure (propensity score) model.
#' @param weights an optional numeric vector of ‘observation weights’ to be used
#'  in the fitting process.
#' @param cov.balancing character string indicating which covariate balancing
#' method should be applied. Possible values are \code{"CBPS"} or \code{"IPT"}.
#' @param cov.balancing.norm character string indicating whether to normalize the estimator. Results in NIPW.
#' @param treatment.effect the treatment effect estimand of interest.
#' Can be either \code{"ATE"} or \code{"ATT"}. For the ATT, IPT and CBPS use the same
#' moment conditions (the weighted covariate means of the controls equal those of the
#' treated), so both give the same estimate, and the weights are automatically normalized.
#' @return List of fit parameters, which is used to derive an object of class
#' \code{teffects} when called by \link[teffects2]{teffect}.
#' @examples
#' # generate some data
#' N <- 50
#' X <- rnorm(N) # confounder
#' A <- rbinom(N, 1, plogis(X)) # treatment variable
#' Y <- X + 0.5 * A # continuous outcome
#' Z <- rbinom(N, 1, plogis(Y)) # binary outcome
#' df <- data.frame(X = X, A = A, Y = Y, Z = Z)
#'
#' teffect(A ~ X, Y ~ 1, data = df, method = "IPW", cov.balancing = "CBPS")
#' teffect(A ~ X, Z ~ 1, data = df, method = "IPW", cov.balancing = "IPT", cov.balancing.norm = TRUE)
#' teffect(A ~ X, Y ~ 1, data = df, method = "IPW", cov.balancing = "IPT", treatment.effect = "ATT")
#' 
#' @importFrom stats quasibinomial family plogis
#' @import geex
#' @export
fit_ipw_cb <- function(X, Y, Zx, weights = rep(1, N), cov.balancing, cov.balancing.norm,
                       treatment.effect = "ATE") {
  N <- NROW(X)
  if (is.null(weights)) {
    weights <- rep.int(1, N)
  }
  wt <- as.vector(weights)
  df <- data.frame(X = X, Y = Y, wt = wt)
  df$Zx <- I(Zx)
  
  #Standardization of non-binary covariates, propensity scores and treatment effect estimates remain the same, only the coeficients of the exposure model are naturally rescaled.
  for (col in colnames(df$Zx)){
    unique_values <- unique(df$Zx[, col])
    if (length(unique_values) > 2){
      df$Zx[, col] <- scale(df$Zx[, col])
    }
  }
  Zx <- df$Zx
  
  initial.mod <- stats::glm.fit(Zx, X, family = quasibinomial(), weights = wt)
  initial.coef <- initial.mod$coefficients

  if (treatment.effect == "ATT") {
    return(fit_ipw_cb_att(df, X, Y, Zx, wt, initial.coef, cov.balancing.norm))
  }
  
  #Propensity Score Estimation with CBPS
  if (cov.balancing == "CBPS"){
  cbps_moment <- function(data){
    X <- data$X
    Zx <- data$Zx
    function(theta){
      Xhat <- Zx %*% theta
      px <- plogis(Xhat)
      c(((X - px)/(px*(1-px)))%*%Zx)  
    }
  }
  cbps_results <- m_estimate(
    estFUN = cbps_moment,
    data   = df,
    weights = wt,
    root_control = setup_root_control(start = initial.coef)
  )
  bips <- coef(cbps_results)
  ps <- plogis(Zx %*% bips)
  ps_c <- ps
  }
  
  #Propensity Score Estimation with IPT
  if (cov.balancing == "IPT"){
  ipt_moment_1 <- function(data){
    X <- data$X
    Zx <- data$Zx
    function(theta){
      Xhat <- Zx %*% theta
      px <- plogis(Xhat)
      c((X/(px))%*%Zx-Zx)
      
    }
  }
  ipt_results_1 <- m_estimate(
    estFUN = ipt_moment_1,
    data   = df,
    weights = wt,
    root_control = setup_root_control(start = initial.coef)
  )
  bips_1 <- coef(ipt_results_1)
  ps <- plogis(Zx %*% bips_1)
  
  ipt_moment_0 <- function(data){
    X <- data$X
    Zx <- data$Zx
    function(theta){
      Xhat <- Zx %*% theta
      px <- plogis(Xhat)
      c(((1-X)/(1-px))%*%Zx-Zx)  
    }
  }
  ipt_results_0 <- m_estimate(
    estFUN = ipt_moment_0,
    data   = df,
    weights = wt,
    root_control = setup_root_control(start = initial.coef)
  )
  bips_0 <- coef(ipt_results_0)
  ps_c <- plogis(Zx %*% bips_0)
  }
 
  X_ps <- X / ps
  cX_ps <- (1 - X) / (1 - ps_c)
  
  
  ## Unnormalized
  Po1.mean <- sum(Y * X_ps * wt) / sum(wt)
  Po0.mean <- sum(Y * cX_ps * wt) / sum(wt)
  
  ## Normalized -- NIPW
  if (cov.balancing.norm == TRUE){
  N.wt1 <- sum(wt * X_ps)
  N.wt0 <- sum(wt * cX_ps)
  Po1.mean <- sum(Y * X_ps * wt) / N.wt1
  Po0.mean <- sum(Y * cX_ps * wt) / N.wt0
  }
  TE <- Po1.mean - Po0.mean
  
  # Moment Conditions for ATE with CBPS
  if (cov.balancing == "CBPS"){
    eq_cbps <- expression(((X - px)/(px*(1-px)))%*%Zx)
    eq_mu1 <- expression((X/px)*Y-mu1)
    eq_mu0 <- expression(((1-X)/(1-px))*Y-mu0)
    #NIPW
    if (cov.balancing.norm == TRUE){
    eq_mu1 <- expression((X/px)*(Y-mu1))
    eq_mu0 <- expression(((1-X)/(1-px))*(Y-mu0))
    }
    eq_ate <- expression(ate - mu1 + mu0)
    
  ipw_ate_moments <- function(data){
    X <- data$X
    Zx <- data$Zx
    Y <- data$Y
    function(theta){
      Xhat <- Zx %*% theta[1:ncol(Zx)]
      mu1 <- theta[ncol(Zx)+1]
      mu0 <- theta[ncol(Zx)+2]
      ate <- theta[ncol(Zx)+3]
      px <- plogis(Xhat)
      c(eval(eq_cbps), eval(eq_mu1), eval(eq_mu0), eval(eq_ate))
    }
  }
  
  ipw_ate_estimation <- m_estimate(
    estFUN = ipw_ate_moments,
    data   = df,
    weights = wt,
    root_control = setup_root_control(start = c(bips,Po1.mean,Po0.mean,TE)))
  }
  
  # Moment conditions with IPT
  if (cov.balancing == "IPT"){
    eq_ipt_1 <- expression((X/(px_1))%*%Zx-Zx)
    eq_ipt_0 <- expression(((1-X)/(1-px_0))%*%Zx-Zx)
    eq_mu1 <- expression((X/px_1)*Y-mu1)
    eq_mu0 <- expression(((1-X)/(1-px_0))*Y-mu0)
    # NIPW
    if (cov.balancing.norm == TRUE){
    eq_mu1 <- expression((X/px_1)*(Y-mu1))
    eq_mu0 <- expression(((1-X)/(1-px_0))*(Y-mu0))
    }
    eq_ate <- expression(ate - mu1 + mu0)
  
  ipw_ate_moments <- function(data){
    X <- data$X
    Zx <- data$Zx
    Y <- data$Y
    function(theta){
      Xhat_1 <- Zx %*% theta[1:ncol(Zx)]
      Xhat_0 <- Zx %*% theta[(ncol(Zx)+1):(2*ncol(Zx))]
      mu1 <- theta[2*ncol(Zx)+1]
      mu0 <- theta[2*ncol(Zx)+2]
      ate <- theta[2*ncol(Zx)+3]
      px_1 <- plogis(Xhat_1)
      px_0 <- plogis(Xhat_0)
      c(eval(eq_ipt_1),eval(eq_ipt_0), eval(eq_mu1), eval(eq_mu0), eval(eq_ate))
    }
  }
  ipw_ate_estimation <- m_estimate(
    estFUN = ipw_ate_moments,
    data   = df,
    weights = wt,
    root_control = setup_root_control(start = c(bips_1,bips_0,Po1.mean,Po0.mean,TE)))  
  }
  
  
  
  
  
  coef_ate <- coef(ipw_ate_estimation)
  ate <- coef_ate[length(coef_ate)]
  vcov_ate <- vcov(ipw_ate_estimation)
  var_ate <- vcov_ate[nrow(vcov_ate),ncol(vcov_ate)]
  mu1 <- coef_ate[length(coef_ate)-2]
  mu0 <- coef_ate[length(coef_ate)-1]
  var_mu1 <- vcov_ate[nrow(vcov_ate)-2,ncol(vcov_ate)-2]
  var_mu0 <- vcov_ate[nrow(vcov_ate)-1,ncol(vcov_ate)-1]
  
  a <- list()
  a$coefs <- c(ATE = ate)
  a$std.err <- c(ATE = sqrt(var_ate))
  a$Wald <- c(ATE = ate^2 / var_ate)
  a$PO.means <- c(Po1 = mu1, Po0 = mu0)
  a$PO.std.err <- c(Po1 = sqrt(var_mu1), Po0 = sqrt(var_mu0))
  a
}


# IPW estimator of the ATT with covariate balancing weights (internal).
# The propensity score solves the ATT balancing conditions
#   sum_i [ (1 - W_i) p_i / (1 - p_i) - W_i ] X_i = 0,
# which are the IPT conditions for the ATT and also the CBPS conditions for the ATT
# (Imai and Ratkovic, 2014). The controls are weighted by p/(1-p).
fit_ipw_cb_att <- function(df, X, Y, Zx, wt, initial.coef, cov.balancing.norm) {
  bips <- fit_att_ps_cb(df, wt, initial.coef)
  ps <- as.vector(plogis(Zx %*% bips))
  om <- (1 - X) * ps / (1 - ps)

  Po1.mean <- sum(Y * X * wt) / sum(X * wt)
  Po0.mean <- sum(Y * om * wt) / sum(X * wt)
  if (cov.balancing.norm == TRUE){
    Po0.mean <- sum(Y * om * wt) / sum(om * wt)
  }
  TE <- Po1.mean - Po0.mean

  eq_ps <- expression(((1 - X) * px / (1 - px) - X) %*% Zx)
  eq_mu1 <- expression(X * (Y - mu1))
  eq_mu0 <- expression((1 - X) * px / (1 - px) * Y - X * mu0)
  if (cov.balancing.norm == TRUE){
    eq_mu0 <- expression((1 - X) * px / (1 - px) * (Y - mu0))
  }
  eq_att <- expression(att - mu1 + mu0)

  ipw_att_moments <- function(data){
    X <- data$X
    Zx <- data$Zx
    Y <- data$Y
    function(theta){
      px <- plogis(Zx %*% theta[1:ncol(Zx)])
      mu1 <- theta[ncol(Zx) + 1]
      mu0 <- theta[ncol(Zx) + 2]
      att <- theta[ncol(Zx) + 3]
      c(eval(eq_ps), eval(eq_mu1), eval(eq_mu0), eval(eq_att))
    }
  }
  est <- m_estimate(
    estFUN = ipw_att_moments,
    data   = df,
    weights = wt,
    root_control = setup_root_control(start = c(bips, Po1.mean, Po0.mean, TE)))

  cf <- coef(est)
  V <- vcov(est)
  k <- length(cf)
  att <- cf[k]
  var_att <- V[k, k]

  a <- list()
  a$coefs <- c(ATT = unname(att))
  a$std.err <- c(ATT = sqrt(var_att))
  a$Wald <- c(ATT = unname(att^2 / var_att))
  a$PO.means <- c(Po1 = unname(cf[k - 2]), Po0 = unname(cf[k - 1]))
  a$PO.std.err <- c(Po1 = sqrt(V[k - 2, k - 2]), Po0 = sqrt(V[k - 1, k - 1]))
  a
}


# Propensity score for the ATT with covariate balancing (internal).
# Solves sum_i [ (1 - W_i) p_i / (1 - p_i) - W_i ] X_i = 0 for the logit index,
# i.e. the odds-weighted covariate means of the controls equal those of the treated.
# These are the IPT conditions for the ATT and also the CBPS conditions for the ATT.
fit_att_ps_cb <- function(df, wt, start) {
  att_ps_moment <- function(data){
    X <- data$X
    Zx <- data$Zx
    function(theta){
      px <- plogis(Zx %*% theta)
      c(((1 - X) * px / (1 - px) - X) %*% Zx)
    }
  }
  ps_results <- m_estimate(
    estFUN = att_ps_moment,
    data   = df,
    weights = wt,
    root_control = setup_root_control(start = start)
  )
  coef(ps_results)
}
