#' Fit Augmented Inverse Probability Weighting using Covariate Balancing Methods: IPT, Graham et al. (2012), or CBPS, Imai and Ratkovic (2014)
#'
#' Fitting function called by \link[teffects2]{teffect} when the method is set
#' to \code{'AIPW'} and \code{'cov.balancing'} is set to \code{'IPT'} or \code{'CBPS'}.
#' 
#' Augmented Inverse Probability Weighting combines the IPW method and the RA
#' method.
#' An exposure model to obtain observation weights is fitted using either:
#' \code{'IPT'}, Inverse Probability Tilting according to Graham et al. (2012) which solves a method-of-moments problem to derive estimates
#' for the exposure (propensity score model).
#' \code{'CBPS'}, Covariate balancing according to Imai and Ratkovic (2014) which effectively weighs the general maximum-likelihood estimation
#' moment conditions by their estimated inverse conditional variance to estimate the exposure (propensity score model).
#' Also a (\link[stats]{glm}) outcome model is fitted to each treatment
#' subgroup, and used to predict potential outcomes.
#'
#' @param X numeric vector of (0,1) specifying the treatment variable.
#' @param Y numeric vector specifying the outcome variable
#' @param Zx design matrix for the exposure (propensity score) model.
#' @param Zy design matrix for the outcome model to be fitted to each treatment
#' subgroup.
#' @param Ofam family function for the outcome model.
#' See \link[stats]{family} for details of family functions.
#' @param treatment.effect the treatment effect estimand of interest.
#' Can be either \code{"ATE"}, \code{"ATT"}, or \code{"All"} (ATE and ATT).
#' The ATC is not implemented with covariate balancing. For the ATT, IPT and CBPS use
#' the same moment conditions and give the same estimate.
#' @param weights an optional numeric vector of ‘observation weights’ to be used
#' in the fitting process.
#' @param cov.balancing character string indicating which covariate balancing
#' method should be applied. Possible values are \code{"CBPS"} or \code{"IPT"}.
#' @param cov.balancing.norm character string indicating whether to normalize the estimator. Results in NAIPW.
#'
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
#' # Fit AIPW using a covariate balancing method to estimate propensity scores
#' teffect(A ~ X, Y ~ X, data = df, method = "AIPW", cov.balancing = "CBPS")
#' teffect(A ~ X, Y ~ X, data = df, method = "AIPW", cov.balancing = "IPT", treatment.effect = "ATT")
#' teffect(A ~ X, Z ~ X, data = df, method = "AIPW", outcome.family = "binomial", cov.balancing = "IPT", cov.balancing.norm = TRUE)
#' 
#' @importFrom stats gaussian quasibinomial family
#' @import geex
#' @export
fit_aipw_cb <- function(X, Y, Zx, Zy, Ofam = gaussian(), treatment.effect = "ATE",
                        weights = rep(1, N), cov.balancing, cov.balancing.norm) {
  
  #Ofam <- binomial(link = "probit")
  link <- Ofam$linkinv # invlink function for outcome model
  mu.eta <- Ofam$mu.eta # derivative of invlink function for outcome model
  var.glm <- Ofam$variance
  
  
  N <- NROW(X)
  if (is.null(weights)) {
    weights <- rep.int(1, N)
  }
  wt <- as.vector(weights)
  df <- data.frame(X = X, Y = Y, weights = wt)
  df$Zx <- I(Zx)
  df$Zy <- I(Zy)
  
  # Standardization of non-binary covariates, propensity scores and treatment effect estimates remain the same, only the coeficients of the exposure model are naturally rescaled.
  for (col in colnames(df$Zx)){
    unique_values <- unique(df$Zx[, col])
    if (length(unique_values) > 2){
      df$Zx[, col] <- scale(df$Zx[, col])
    }
  }
  Zx <- df$Zx
  for (col in colnames(df$Zy)){
    unique_values <- unique(df$Zy[, col])
    if (length(unique_values) > 2){
      df$Zy[, col] <- scale(df$Zy[, col])
    }
  }
  Zy <- df$Zy
  
  initial.X.mod <- stats::glm.fit(Zx, X, family = quasibinomial(), weights = wt)
  initial.Y.mod1 <- stats::glm.fit(Zy, Y, family = Ofam, weights = wt * X)
  initial.Y.mod0 <- stats::glm.fit(Zy, Y, family = Ofam, weights = wt * (1 - X))
  
  initial.X.coef <- initial.X.mod$coefficients
  
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
      root_control = setup_root_control(start = initial.X.coef)
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
      root_control = setup_root_control(start = initial.X.coef)
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
      root_control = setup_root_control(start = initial.X.coef)
    )
    bips_0 <- coef(ipt_results_0)
    ps_c <- plogis(Zx %*% bips_0)
  }
  
  eta1 <- initial.Y.mod1$linear.predictor
  bips_y_1 <- initial.Y.mod1$coefficients
  eta0 <- initial.Y.mod0$linear.predictor
  bips_y_0 <- initial.Y.mod0$coefficients
  
  initial.mu1 <- link(eta1)
  initial.mu0 <- link(eta0)
  
  X_ps <- X / ps
  cX_ps <- (1 - X) / (1 - ps_c)

  # Get potential outcomes
  Po1 <- X_ps * Y + (1 - X_ps) * initial.mu1
  Po0 <- cX_ps * Y + (1 - cX_ps) * initial.mu0
  
  if (!treatment.effect %in% c("ATE", "ATT", "ATC", "All")) {
    stop("'treatment.effect' not recognized")
  }
  TEs <- list(ATE = "ATE", ATT = "ATT", ATC = "ATC")
  
  if (treatment.effect != "All") TEs <- TEs[treatment.effect]
  if ("ATC" %in% names(TEs)) {
    if (treatment.effect == "ATC") stop("ATC not implemented with covariate balancing.")
    warning("ATC not implemented with covariate balancing. Estimating ATE and ATT.")
    TEs <- TEs[c("ATE", "ATT")]
  }
  
  run <- lapply(TEs, function(teffect) {

  ## ATE 
  if (teffect == "ATE"){
    Po1.mean <- sum(Po1 * wt) / sum(wt)
    Po0.mean <- sum(Po0 * wt) / sum(wt)
    # NIPW
    if (cov.balancing.norm == TRUE){
    N.wt1 <- sum(wt * X_ps)
    N.wt0 <- sum(wt * cX_ps)
    Po1.mean <- sum(Po1 * wt) / N.wt1
    Po0.mean <- sum(Po0 * wt) / N.wt0}
    
    TE <- Po1.mean - Po0.mean
    # Moment Conditions for ATE with CBPS
    if (cov.balancing == "CBPS"){
      eq_cbps <- expression(((X - px)/(px*(1-px)))%*%Zx)
      eq_y1 <- expression((X*((Y-initial.est.mu1)/var.glm1)*mu.eta1)%*%Zy)
      eq_y0 <- expression(((1-X)*((Y-initial.est.mu0)/var.glm0)*mu.eta0)%*%Zy)
      eq_mu1 <- expression((X/px)*Y + (1-(X/px))*initial.est.mu1 - mu1)
      eq_mu0 <- expression(((1-X)/(1-px))*Y + (1-(1-X)/(1-px))*initial.est.mu0 - mu0)
      #NIPW
      if (cov.balancing.norm == TRUE){
      eq_mu1 <- expression(((X/px)*Y + (1-(X/px))*initial.est.mu1)/(N.wt1/sum(wt)) - mu1)
      eq_mu0 <- expression((((1-X)/(1-px))*Y + (1-(1-X)/(1-px))*initial.est.mu0)/(N.wt0/sum(wt)) - mu0)
      }
      eq_ate <- expression(ate - mu1 + mu0)
      
      aipw_ate_moments <- function(data){
        X <- data$X
        Zx <- data$Zx
        Y <- data$Y
        Zy <- data$Zy
        function(theta){
          Xhat <- Zx %*% theta[1:ncol(Zx)]
          Yhat1 <- Zy %*% theta[(ncol(Zx)+1):(ncol(Zx)+ncol(Zy))]
          Yhat0 <- Zy %*% theta[(ncol(Zx)+ncol(Zy)+1):(ncol(Zx)+2*ncol(Zy))]
          mu1 <- theta[(ncol(Zx)+2*ncol(Zy)+1)]
          mu0 <- theta[(ncol(Zx)+2*ncol(Zy)+2)]
          ate <- theta[(ncol(Zx)+2*ncol(Zy)+3)]
          px <- plogis(Xhat)
          initial.est.mu1 <- link(Yhat1)
          initial.est.mu0 <- link(Yhat0)
          mu.eta1 <- mu.eta(Yhat1)
          mu.eta0 <- mu.eta(Yhat0)
          var.glm1 <- var.glm(initial.est.mu1)
          var.glm0 <- var.glm(initial.est.mu0)
          c(eval(eq_cbps), eval(eq_y1), eval(eq_y0), eval(eq_mu1), eval(eq_mu0), eval(eq_ate))
        }
      }
      
      aipw_ate_estimation <- m_estimate(
        estFUN = aipw_ate_moments,
        data   = df,
        weights = wt,
        root_control = setup_root_control(start = c(bips, bips_y_1,bips_y_0,Po1.mean,Po0.mean,TE)))
    }
    
    if (cov.balancing == "IPT"){
      eq_ipt_1 <- expression((X/(px_1))%*%Zx-Zx)
      eq_ipt_0 <- expression(((1-X)/(1-px_0))%*%Zx-Zx)
      eq_y1 <- expression((X*((Y-initial.est.mu1)/var.glm1)*mu.eta1)%*%Zy)
      eq_y0 <- expression(((1-X)*((Y-initial.est.mu0)/var.glm0)*mu.eta0)%*%Zy)
      eq_mu1 <- expression((X/px_1)*Y + (1-(X/px_1))*initial.est.mu1 - mu1)
      eq_mu0 <- expression(((1-X)/(1-px_0))*Y + (1-(1-X)/(1-px_0))*initial.est.mu0 - mu0)
      #NIPW
      if (cov.balancing.norm == TRUE){
        eq_mu1 <- expression(((X/px_1)*Y + (1-(X/px_1))*initial.est.mu1)/(N.wt1/sum(wt)) - mu1)
        eq_mu0 <- expression((((1-X)/(1-px_0))*Y + (1-(1-X)/(1-px_0))*initial.est.mu0)/(N.wt0/sum(wt)) - mu0)
      }
      eq_ate <- expression(ate - mu1 + mu0)
      
      aipw_ate_moments <- function(data){
        X <- data$X
        Zx <- data$Zx
        Y <- data$Y
        Zy <- data$Zy
        function(theta){
          Xhat_1 <- Zx %*% theta[1:ncol(Zx)]
          Xhat_0 <- Zx %*% theta[(ncol(Zx)+1):(2*ncol(Zx))]
          Yhat1 <- Zy %*% theta[(2*ncol(Zx)+1):(2*ncol(Zx)+ncol(Zy))]
          Yhat0 <- Zy %*% theta[(2*ncol(Zx)+ncol(Zy)+1):(2*ncol(Zx)+2*ncol(Zy))]
          mu1 <- theta[(2*ncol(Zx)+2*ncol(Zy)+1)]
          mu0 <- theta[(2*ncol(Zx)+2*ncol(Zy)+2)]
          ate <- theta[(2*ncol(Zx)+2*ncol(Zy)+3)]
          px_1 <- plogis(Xhat_1)
          px_0 <- plogis(Xhat_0)
          initial.est.mu1 <- link(Yhat1)
          initial.est.mu0 <- link(Yhat0)
          mu.eta1 <- mu.eta(Yhat1)
          mu.eta0 <- mu.eta(Yhat0)
          var.glm1 <- var.glm(initial.est.mu1)
          var.glm0 <- var.glm(initial.est.mu0)
          c(eval(eq_ipt_1),eval(eq_ipt_0), eval(eq_y1), eval(eq_y0), eval(eq_mu1), eval(eq_mu0), eval(eq_ate))
        }
      }
      
      aipw_ate_estimation <- m_estimate(
        estFUN = aipw_ate_moments,
        data   = df,
        weights = wt,
        root_control = setup_root_control(start = c(bips_1,bips_0, bips_y_1,bips_y_0,Po1.mean,Po0.mean,TE)))
    }
    
    coef_ate <- coef(aipw_ate_estimation)
    ate <- coef_ate[length(coef_ate)]
    vcov_ate <- vcov(aipw_ate_estimation)
    var_ate <- vcov_ate[nrow(vcov_ate),ncol(vcov_ate)]
    mu1 <- coef_ate[length(coef_ate)-2]
    mu0 <- coef_ate[length(coef_ate)-1]
    var_mu1 <- vcov_ate[nrow(vcov_ate)-2,ncol(vcov_ate)-2]
    var_mu0 <- vcov_ate[nrow(vcov_ate)-1,ncol(vcov_ate)-1]
    }
    
    ## ATT: the ATT propensity score is the same for IPT and CBPS.
    ## mu0 = mean over treated of m0(X) + (odds-weighted mean residual of the controls),
    ## with m0 an (unweighted) GLM fit on the controls.
    if (teffect == "ATT"){
      bips_att <- fit_att_ps_cb(df, wt, initial.X.coef)
      ps_att <- as.vector(plogis(Zx %*% bips_att))
      om <- (1 - X) * ps_att / (1 - ps_att)
      m0 <- as.vector(link(eta0))
      N1 <- sum(wt * X)
      Po1.mean <- sum(wt * X * Y) / N1
      r0 <- sum(wt * om * (Y - m0)) / N1
      if (cov.balancing.norm == TRUE) r0 <- sum(wt * om * (Y - m0)) / sum(wt * om)
      Po0.mean <- sum(wt * X * m0) / N1 + r0
      TE <- Po1.mean - Po0.mean

      eq_ps <- expression(((1 - X) * px / (1 - px) - X) %*% Zx)
      eq_y0 <- expression(((1 - X) * ((Y - initial.est.mu0) / var.glm0) * mu.eta0) %*% Zy)
      eq_mu1 <- expression(X * (Y - mu1))
      eq_r0 <- expression((1 - X) * px / (1 - px) * (Y - initial.est.mu0) - X * r0)
      if (cov.balancing.norm == TRUE){
        eq_r0 <- expression((1 - X) * px / (1 - px) * (Y - initial.est.mu0 - r0))
      }
      eq_mu0 <- expression(X * (initial.est.mu0 + r0 - mu0))
      eq_att <- expression(att - mu1 + mu0)

      aipw_att_moments <- function(data){
        X <- data$X
        Zx <- data$Zx
        Y <- data$Y
        Zy <- data$Zy
        function(theta){
          kx <- ncol(Zx); ky <- ncol(Zy)
          px <- plogis(Zx %*% theta[1:kx])
          Yhat0 <- Zy %*% theta[(kx + 1):(kx + ky)]
          mu1 <- theta[kx + ky + 1]
          r0 <- theta[kx + ky + 2]
          mu0 <- theta[kx + ky + 3]
          att <- theta[kx + ky + 4]
          initial.est.mu0 <- link(Yhat0)
          mu.eta0 <- mu.eta(Yhat0)
          var.glm0 <- var.glm(initial.est.mu0)
          c(eval(eq_ps), eval(eq_y0), eval(eq_mu1), eval(eq_r0), eval(eq_mu0), eval(eq_att))
        }
      }
      aipw_att_estimation <- m_estimate(
        estFUN = aipw_att_moments,
        data   = df,
        weights = wt,
        root_control = setup_root_control(start = c(bips_att, bips_y_0, Po1.mean, r0, Po0.mean, TE)))

      coef_ate <- coef(aipw_att_estimation)
      vcov_ate <- vcov(aipw_att_estimation)
      k <- length(coef_ate)
      ate <- coef_ate[k]
      var_ate <- vcov_ate[k, k]
      mu1 <- coef_ate[k - 3]
      mu0 <- coef_ate[k - 1]
      var_mu1 <- vcov_ate[k - 3, k - 3]
      var_mu0 <- vcov_ate[k - 1, k - 1]
      }
    
    a <- list()
    a$coefs <- c(ATE = ate)
    a$std.err <- c(ATE = sqrt(var_ate))
    a$Wald <- c(ATE = ate^2 / var_ate)
    a$PO.means <- c(Po1 = mu1, Po0 = mu0)
    a$PO.std.err <- c(Po1 = sqrt(var_mu1), Po0 = sqrt(var_mu0))
    a
  })
  
  run <- simplify2array(run)
  
  a <- list(
    coef = as.double(run["coefs", ]),
    std.err = as.double(run["std.err", ]),
    Wald = as.double(run["Wald", ]),
    Po.means = run["PO.means", ],
    Po.std.err = run["PO.std.err", ]
  )
  names(a$coef) <- names(a$std.err) <- names(a$Wald) <- names(TEs)
  a
}
