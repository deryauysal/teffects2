# Install once (the R package is in the R/ folder of the repository):
# devtools::install_github("deryauysal/teffects2", subdir = "R")
library(teffects2)

set.seed(1)


N <- 1000
X <- rnorm(N, mean = 0, sd = 1) # confounder
X_1 <- rnorm(N, mean = 1, sd = 5) # spurious covariate
A <- rbinom(N, 1, plogis(X)) # treatment variable
Y <- X + 0.5 * A # continuous outcome
Z <- rbinom(N, 1, plogis(Y)) # binary outcome
df <- data.frame(X = X, A = A, Y = Y, Z = Z, X_1 = X_1)

weights <- c(rep(1, 500), rep(1, 500))




# As described in the paper we expect the AIPW, IPW, and IPWRA estimates to be the same when using IPT weights.
model_IPT_IPW <- teffect(A ~ X + X_1, Z ~ X + X_1, data = df, weights = weights, outcome.family = "gaussian",
                    method = "IPW", treatment.effect = "ATE", cov.balancing = "IPT",
                    cov.balancing.norm = TRUE)
print(model_IPT_IPW)


model_IPT_AIPW <- teffect(A ~ X + X_1, Z ~ X + X_1, data = df, weights = weights, outcome.family = "gaussian",
                        method = "AIPW", treatment.effect = "ATE", cov.balancing = "IPT",
                        cov.balancing.norm = TRUE)
print(model_IPT_AIPW)

model_IPT_IPWRA <- teffect(A ~ X + X_1, Z ~ X + X_1, data = df, weights = weights, outcome.family = "gaussian",
                         method = "IPWRA", treatment.effect = "ATE", cov.balancing = "IPT",
                         cov.balancing.norm = TRUE)
print(model_IPT_IPWRA)


# Further we can see that the AIPW and IPW estimates using IPT weights are automatically normalized.
# For IPW using IPT
model_IPT_IPW_nrm <- teffect(A ~ X + X_1, Z ~ X + X_1, data = df, weights = weights, outcome.family = "gaussian",
                        method = "IPW", treatment.effect = "ATE", cov.balancing = "IPT",
                        cov.balancing.norm = TRUE)
print(model_IPT_IPW_nrm)

model_IPT_IPW_not_nrm <- teffect(A ~ X + X_1, Z ~ X + X_1, data = df, weights = weights, outcome.family = "gaussian",
                            method = "IPW", treatment.effect = "ATE", cov.balancing = "IPT",
                            cov.balancing.norm = FALSE)
print(model_IPT_IPW_not_nrm)


# For AIPW using IPT
model_IPT_AIPW_nrm <- teffect(A ~ X + X_1, Z ~ X + X_1, data = df, weights = weights, outcome.family = "gaussian",
                            method = "AIPW", treatment.effect = "ATE", cov.balancing = "IPT",
                            cov.balancing.norm = TRUE)
print(model_IPT_AIPW_nrm)

model_IPT_AIPW_not_nrm <- teffect(A ~ X + X_1, Z ~ X + X_1, data = df, weights = weights, outcome.family = "gaussian",
                                method = "AIPW", treatment.effect = "ATE", cov.balancing = "IPT",
                                cov.balancing.norm = FALSE)
print(model_IPT_AIPW_not_nrm)

# Normalization does affect estimates when using CBPS
model_CBPS_AIPW_nrm <- teffect(A ~ X + X_1, Z ~ X + X_1, data = df, weights = weights, outcome.family = "gaussian",
                             method = "AIPW", treatment.effect = "ATE", cov.balancing = "CBPS",
                             cov.balancing.norm = TRUE)
print(model_CBPS_AIPW_nrm)

model_CBPS_AIPW_not_nrm <- teffect(A ~ X + X_1, Z ~ X + X_1, data = df, weights = weights, outcome.family = "gaussian",
                                 method = "AIPW", treatment.effect = "ATE", cov.balancing = "CBPS",
                                 cov.balancing.norm = FALSE)
print(model_CBPS_AIPW_not_nrm)

