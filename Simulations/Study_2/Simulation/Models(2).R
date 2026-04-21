################################################################################
# File:         Models(2)_recalibrated.R
# Description:  Recalibrated population models for Study 2 (5-factor model).
#               Adjusts theta for endogenous indicators and psi55 for R^2 ~ 0.30.
#               Calibrated using GenerateData() with N = 500,000, averaged over
#               5 runs.
# Dependencies: none (pure data definitions)
# Used by:      Simulation(2).R

all_models <- list()

# full models (nonlinear terms included)

all_models[["normal_rel04"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15

# structural model (intercepts calibrated for E[eta4] = E[eta5] = 0)
eta4 ~ -0.215*1 + 0.20*eta1 + 0.20*eta2 + 0.20*eta3
     + 0.11*eta1:eta2 + 0.11*eta1:eta3
     + 0.08*eta1:eta1 + 0.08*eta2:eta2

eta5 ~ -0.161*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
     + 0.08*eta1:eta4 + 0.08*eta2:eta4
     + 0.06*eta1:eta1 + 0.06*eta3:eta3

# variances and covariances (exogenous)
eta1 ~~ 0.3*eta2
eta2 ~~ 0.3*eta3
eta1 ~~ 0.2*eta3
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3

# structural residual variances (calibrated for R^2 approx. 0.30)
eta4 ~~ 0.530*eta4
eta5 ~~ 0.558*eta5

# indicator residual variances (rel = 0.4)
x1 ~~ 1.500*x1;  x2 ~~ 1.500*x2;  x3 ~~ 1.500*x3
x4 ~~ 1.500*x4;  x5 ~~ 1.500*x5;  x6 ~~ 1.500*x6
x7 ~~ 1.500*x7;  x8 ~~ 1.500*x8;  x9 ~~ 1.500*x9
x10 ~~ 1.135*x10; x11 ~~ 1.135*x11; x12 ~~ 1.135*x12
x13 ~~ 1.197*x13; x14 ~~ 1.197*x14; x15 ~~ 1.197*x15
"

all_models[["normal_rel06"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15

# structural model (intercepts calibrated for E[eta4] = E[eta5] = 0)
eta4 ~ -0.215*1 + 0.20*eta1 + 0.20*eta2 + 0.20*eta3
     + 0.11*eta1:eta2 + 0.11*eta1:eta3
     + 0.08*eta1:eta1 + 0.08*eta2:eta2

eta5 ~ -0.161*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
     + 0.08*eta1:eta4 + 0.08*eta2:eta4
     + 0.06*eta1:eta1 + 0.06*eta3:eta3

# variances and covariances (exogenous)
eta1 ~~ 0.3*eta2
eta2 ~~ 0.3*eta3
eta1 ~~ 0.2*eta3
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3

# structural residual variances (calibrated for R^2 approx. 0.30)
eta4 ~~ 0.530*eta4
eta5 ~~ 0.558*eta5

# indicator residual variances (rel = 0.6)
x1 ~~ 0.667*x1;  x2 ~~ 0.667*x2;  x3 ~~ 0.667*x3
x4 ~~ 0.667*x4;  x5 ~~ 0.667*x5;  x6 ~~ 0.667*x6
x7 ~~ 0.667*x7;  x8 ~~ 0.667*x8;  x9 ~~ 0.667*x9
x10 ~~ 0.505*x10; x11 ~~ 0.505*x11; x12 ~~ 0.505*x12
x13 ~~ 0.532*x13; x14 ~~ 0.532*x14; x15 ~~ 0.532*x15
"

all_models[["normal_rel08"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15

# structural model (intercepts calibrated for E[eta4] = E[eta5] = 0)
eta4 ~ -0.215*1 + 0.20*eta1 + 0.20*eta2 + 0.20*eta3
     + 0.11*eta1:eta2 + 0.11*eta1:eta3
     + 0.08*eta1:eta1 + 0.08*eta2:eta2

eta5 ~ -0.161*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
     + 0.08*eta1:eta4 + 0.08*eta2:eta4
     + 0.06*eta1:eta1 + 0.06*eta3:eta3

# variances and covariances (exogenous)
eta1 ~~ 0.3*eta2
eta2 ~~ 0.3*eta3
eta1 ~~ 0.2*eta3
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3

# structural residual variances (calibrated for R^2 approx. 0.30)
eta4 ~~ 0.530*eta4
eta5 ~~ 0.558*eta5

# indicator residual variances (rel = 0.8)
x1 ~~ 0.250*x1;  x2 ~~ 0.250*x2;  x3 ~~ 0.250*x3
x4 ~~ 0.250*x4;  x5 ~~ 0.250*x5;  x6 ~~ 0.250*x6
x7 ~~ 0.250*x7;  x8 ~~ 0.250*x8;  x9 ~~ 0.250*x9
x10 ~~ 0.189*x10; x11 ~~ 0.189*x11; x12 ~~ 0.189*x12
x13 ~~ 0.199*x13; x14 ~~ 0.199*x14; x15 ~~ 0.199*x15
"

# null models (no interaction and quadratic effects)

all_models[["null_model_rel04"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15

# structural model (same intercepts, nonlinear terms = 0)
eta4 ~ -0.215*1 + 0.20*eta1 + 0.20*eta2 + 0.20*eta3
     + 0*eta1:eta2 + 0*eta1:eta3
     + 0*eta1:eta1 + 0*eta2:eta2

eta5 ~ -0.161*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
     + 0*eta1:eta4 + 0*eta2:eta4
     + 0*eta1:eta1 + 0*eta3:eta3

# variances and covariances (exogenous)
eta1 ~~ 0.3*eta2
eta2 ~~ 0.3*eta3
eta1 ~~ 0.2*eta3
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3

# structural residual variances (calibrated for R^2 approx. 0.30)
eta4 ~~ 0.530*eta4
eta5 ~~ 0.558*eta5

# indicator residual variances (rel = 0.4)
x1 ~~ 1.500*x1;  x2 ~~ 1.500*x2;  x3 ~~ 1.500*x3
x4 ~~ 1.500*x4;  x5 ~~ 1.500*x5;  x6 ~~ 1.500*x6
x7 ~~ 1.500*x7;  x8 ~~ 1.500*x8;  x9 ~~ 1.500*x9
x10 ~~ 1.135*x10; x11 ~~ 1.135*x11; x12 ~~ 1.135*x12
x13 ~~ 1.197*x13; x14 ~~ 1.197*x14; x15 ~~ 1.197*x15
"

all_models[["null_model_rel06"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15

# structural model (same intercepts, nonlinear terms = 0)
eta4 ~ -0.215*1 + 0.20*eta1 + 0.20*eta2 + 0.20*eta3
     + 0*eta1:eta2 + 0*eta1:eta3
     + 0*eta1:eta1 + 0*eta2:eta2

eta5 ~ -0.161*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
     + 0*eta1:eta4 + 0*eta2:eta4
     + 0*eta1:eta1 + 0*eta3:eta3

# variances and covariances (exogenous)
eta1 ~~ 0.3*eta2
eta2 ~~ 0.3*eta3
eta1 ~~ 0.2*eta3
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3

# structural residual variances (calibrated for R^2 approx. 0.30)
eta4 ~~ 0.530*eta4
eta5 ~~ 0.558*eta5

# indicator residual variances (rel = 0.6)
x1 ~~ 0.667*x1;  x2 ~~ 0.667*x2;  x3 ~~ 0.667*x3
x4 ~~ 0.667*x4;  x5 ~~ 0.667*x5;  x6 ~~ 0.667*x6
x7 ~~ 0.667*x7;  x8 ~~ 0.667*x8;  x9 ~~ 0.667*x9
x10 ~~ 0.505*x10; x11 ~~ 0.505*x11; x12 ~~ 0.505*x12
x13 ~~ 0.532*x13; x14 ~~ 0.532*x14; x15 ~~ 0.532*x15
"

all_models[["null_model_rel08"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15

# structural model (same intercepts, nonlinear terms = 0)
eta4 ~ -0.215*1 + 0.20*eta1 + 0.20*eta2 + 0.20*eta3
     + 0*eta1:eta2 + 0*eta1:eta3
     + 0*eta1:eta1 + 0*eta2:eta2

eta5 ~ -0.161*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
     + 0*eta1:eta4 + 0*eta2:eta4
     + 0*eta1:eta1 + 0*eta3:eta3

# variances and covariances (exogenous)
eta1 ~~ 0.3*eta2
eta2 ~~ 0.3*eta3
eta1 ~~ 0.2*eta3
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3

# structural residual variances (calibrated for R^2 approx. 0.30)
eta4 ~~ 0.530*eta4
eta5 ~~ 0.558*eta5

# indicator residual variances (rel = 0.8)
x1 ~~ 0.250*x1;  x2 ~~ 0.250*x2;  x3 ~~ 0.250*x3
x4 ~~ 0.250*x4;  x5 ~~ 0.250*x5;  x6 ~~ 0.250*x6
x7 ~~ 0.250*x7;  x8 ~~ 0.250*x8;  x9 ~~ 0.250*x9
x10 ~~ 0.189*x10; x11 ~~ 0.189*x11; x12 ~~ 0.189*x12
x13 ~~ 0.199*x13; x14 ~~ 0.199*x14; x15 ~~ 0.199*x15
"
