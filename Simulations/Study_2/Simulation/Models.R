# =================================================================
# COMPLETE POPULATION MODELS FOR SIMULATION STUDY
# Three reliability levels: 0.8, 0.6, 0.4
# All maintain R² ≈ 0.30 for endogenous variables
# =================================================================

# RELIABILITY 0.8 POPULATION MODEL

population_model_rel08 <- "
# ---------- Measurement Model ----------
# All loadings fixed to 1 for reliability = 0.8 with measurement error = 0.25
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15
eta6 =~ 1*x16 + 1*x17 + 1*x18

# ---------- Structural Model ----------
# Each equation has its own set of equal coefficients within effect types
# eta4: main effects = 0.21, interaction = 0.13, quadratic = 0.09
eta4 ~ 0.1*1 + 0.21*eta1 + 0.21*eta2 + 0.21*eta3 + 0.13*eta1:eta2 + 0.09*eta1:eta1 

# eta5: main effects = 0.18, interaction = 0.11, quadratic = 0.09  
eta5 ~ 0.1*1 + 0.18*eta4 + 0.18*eta1 + 0.18*eta2 + 0.18*eta3 + 0.11*eta2:eta4 + 0.09*eta2:eta2

# eta6: main effects = 0.15, interaction = 0.10, quadratic = 0.09
eta6 ~ 0.1*1 + 0.15*eta5 + 0.15*eta1 + 0.15*eta2 + 0.15*eta3 + 0.10*eta3:eta5 + 0.09*eta3:eta3

# ---------- Variances and Covariances ----------
# Exogenous variables
eta1 ~~ 0.3*eta2
eta2 ~~ 0.3*eta3
eta1 ~~ 0.2*eta3
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3

# Residual variances adjusted to achieve R² ≈ 0.30
# With different coefficients per equation, we can better control R²
eta4 ~~ 0.68*eta4
eta5 ~~ 0.70*eta5  
eta6 ~~ 0.70*eta6

# ---------- Measurement Errors ----------
# All measurement errors = 0.25 for reliability = 0.8
x1 ~~ 0.25*x1; x2 ~~ 0.25*x2; x3 ~~ 0.25*x3
x4 ~~ 0.25*x4; x5 ~~ 0.25*x5; x6 ~~ 0.25*x6
x7 ~~ 0.25*x7; x8 ~~ 0.25*x8; x9 ~~ 0.25*x9
x10 ~~ 0.25*x10; x11 ~~ 0.25*x11; x12 ~~ 0.25*x12
x13 ~~ 0.25*x13; x14 ~~ 0.25*x14; x15 ~~ 0.25*x15
x16 ~~ 0.25*x16; x17 ~~ 0.25*x17; x18 ~~ 0.25*x18
"

# RELIABILITY 0.6 POPULATION MODEL

population_model_rel06 <- "
# ---------- Measurement Model ----------
# All loadings fixed to 1 for reliability = 0.6 with measurement error = 0.67
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15
eta6 =~ 1*x16 + 1*x17 + 1*x18

# ---------- Structural Model ----------
# Each equation has its own set of equal coefficients within effect types
# eta4: main effects = 0.21, interaction = 0.13, quadratic = 0.09
eta4 ~ 0.1*1 + 0.21*eta1 + 0.21*eta2 + 0.21*eta3 + 0.13*eta1:eta2 + 0.09*eta1:eta1 

# eta5: main effects = 0.18, interaction = 0.11, quadratic = 0.09  
eta5 ~ 0.1*1 + 0.18*eta4 + 0.18*eta1 + 0.18*eta2 + 0.18*eta3 + 0.11*eta2:eta4 + 0.09*eta2:eta2

# eta6: main effects = 0.15, interaction = 0.10, quadratic = 0.09
eta6 ~ 0.1*1 + 0.15*eta5 + 0.15*eta1 + 0.15*eta2 + 0.15*eta3 + 0.10*eta3:eta5 + 0.09*eta3:eta3

# ---------- Variances and Covariances ----------
# Exogenous variables
eta1 ~~ 0.3*eta2
eta2 ~~ 0.3*eta3
eta1 ~~ 0.2*eta3
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3

# Residual variances adjusted to achieve R² ≈ 0.30
# With different coefficients per equation, we can better control R²
eta4 ~~ 0.68*eta4
eta5 ~~ 0.70*eta5  
eta6 ~~ 0.70*eta6

# ---------- Measurement Errors ----------
# All measurement errors = 0.67 for reliability = 0.6
# Reliability = var(eta) / (var(eta) + var(error)) = 1/(1+0.67) ≈ 0.60
x1 ~~ 0.67*x1; x2 ~~ 0.67*x2; x3 ~~ 0.67*x3
x4 ~~ 0.67*x4; x5 ~~ 0.67*x5; x6 ~~ 0.67*x6
x7 ~~ 0.67*x7; x8 ~~ 0.67*x8; x9 ~~ 0.67*x9
x10 ~~ 0.67*x10; x11 ~~ 0.67*x11; x12 ~~ 0.67*x12
x13 ~~ 0.67*x13; x14 ~~ 0.67*x14; x15 ~~ 0.67*x15
x16 ~~ 0.67*x16; x17 ~~ 0.67*x17; x18 ~~ 0.67*x18
"

# RELIABILITY 0.4 POPULATION MODEL

population_model_rel04 <- "
# ---------- Measurement Model ----------
# All loadings fixed to 1 for reliability = 0.4 with measurement error = 1.5
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15
eta6 =~ 1*x16 + 1*x17 + 1*x18

# ---------- Structural Model ----------
# Each equation has its own set of equal coefficients within effect types
# eta4: main effects = 0.21, interaction = 0.13, quadratic = 0.09
eta4 ~ 0.1*1 + 0.21*eta1 + 0.21*eta2 + 0.21*eta3 + 0.13*eta1:eta2 + 0.09*eta1:eta1 

# eta5: main effects = 0.18, interaction = 0.11, quadratic = 0.09  
eta5 ~ 0.1*1 + 0.18*eta4 + 0.18*eta1 + 0.18*eta2 + 0.18*eta3 + 0.11*eta2:eta4 + 0.09*eta2:eta2

# eta6: main effects = 0.15, interaction = 0.10, quadratic = 0.09
eta6 ~ 0.1*1 + 0.15*eta5 + 0.15*eta1 + 0.15*eta2 + 0.15*eta3 + 0.10*eta3:eta5 + 0.09*eta3:eta3

# ---------- Variances and Covariances ----------
# Exogenous variables
eta1 ~~ 0.3*eta2
eta2 ~~ 0.3*eta3
eta1 ~~ 0.2*eta3
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3

# Residual variances adjusted to achieve R² ≈ 0.30
# With different coefficients per equation, we can better control R²
eta4 ~~ 0.68*eta4
eta5 ~~ 0.70*eta5  
eta6 ~~ 0.70*eta6

# ---------- Measurement Errors ----------
# All measurement errors = 1.5 for reliability = 0.4
# Reliability = var(eta) / (var(eta) + var(error)) = 1/(1+1.5) = 0.40
x1 ~~ 1.5*x1; x2 ~~ 1.5*x2; x3 ~~ 1.5*x3
x4 ~~ 1.5*x4; x5 ~~ 1.5*x5; x6 ~~ 1.5*x6
x7 ~~ 1.5*x7; x8 ~~ 1.5*x8; x9 ~~ 1.5*x9
x10 ~~ 1.5*x10; x11 ~~ 1.5*x11; x12 ~~ 1.5*x12
x13 ~~ 1.5*x13; x14 ~~ 1.5*x14; x15 ~~ 1.5*x15
x16 ~~ 1.5*x16; x17 ~~ 1.5*x17; x18 ~~ 1.5*x18
"

# CHECKING

# ANALYSIS MODEL (SAME FOR ALL)

analysis_model <- "
# Measurement model
eta1 =~ x1 + x2 + x3
eta2 =~ x4 + x5 + x6
eta3 =~ x7 + x8 + x9
eta4 =~ x10 + x11 + x12
eta5 =~ x13 + x14 + x15
eta6 =~ x16 + x17 + x18

# Structural model
eta4 ~ eta1 + eta2 + eta3 + eta1:eta2 + eta1:eta1 + eta2:eta2
eta5 ~ eta4 + eta1 + eta2 + eta3 + eta2:eta4 + eta3:eta3
eta6 ~ eta5 + eta1 + eta2 + eta3 + eta1:eta5 + eta2:eta3
"

# generate data and check reliabilities/R² values

set.seed(123)

# 0.8
data_08 <- GenerateData(
  model = population_model_rel08,
  N = 1000,
  skewness = c(0, 0, 0),
  excesskurtosis = c(0, 0, 0),
  distr.exo = "normal",
  distr.zeta = "normal",
  distr.epsilon = "normal",
  add.eta = TRUE,
  return.info = TRUE
)

cat("\nReliabilities (target = 0.80):\n")
rel_08 <- attr(data_08, "observed_reliabilities")
for(i in names(rel_08)) {
  cat(sprintf("  %s: %.3f\n", i, mean(rel_08[[i]])))
}

cat("\nR² values (target = 0.30):\n")
r2_08 <- attr(data_08, "observed_R2")
for(i in names(r2_08)) {
  cat(sprintf("  %s: %.3f\n", i, r2_08[[i]]))
}

# 0.6

data_06 <- GenerateData(
  model = population_model_rel06,
  N = 1000,
  skewness = c(0, 0, 0),
  excesskurtosis = c(0, 0, 0),
  distr.exo = "normal",
  distr.zeta = "normal",
  distr.epsilon = "normal",
  add.eta = FALSE,
  return.info = TRUE
)

cat("\nReliabilities (target = 0.60):\n")
rel_06 <- attr(data_06, "observed_reliabilities")
for(i in names(rel_06)) {
  cat(sprintf("  %s: %.3f\n", i, mean(rel_06[[i]])))
}

cat("\nR² values (target = 0.30):\n")
r2_06 <- attr(data_06, "observed_R2")
for(i in names(r2_06)) {
  cat(sprintf("  %s: %.3f\n", i, r2_06[[i]]))
}

# 0.4

data_04 <- GenerateData(
  model = population_model_rel04,
  N = 1000,
  skewness = c(0, 0, 0),
  excesskurtosis = c(0, 0, 0),
  distr.exo = "normal",
  distr.zeta = "normal",
  distr.epsilon = "normal",
  add.eta = FALSE,
  return.info = TRUE
)

cat("\nReliabilities (target = 0.40):\n")
rel_04 <- attr(data_04, "observed_reliabilities")
for(i in names(rel_04)) {
  cat(sprintf("  %s: %.3f\n", i, mean(rel_04[[i]])))
}

cat("\nR² values (target = 0.30):\n")
r2_04 <- attr(data_04, "observed_R2")
for(i in names(r2_04)) {
  cat(sprintf("  %s: %.3f\n", i, r2_04[[i]]))
}
