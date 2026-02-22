################################################################################
# File:         Models(2).R
# Description:  Population models for Study 2 (5-factor model). Defines lavaan
#               syntax for full and null models at three reliability levels
#               (0.4, 0.6, 0.8). Intercepts calibrated so E[eta4] = E[eta5] = 0.
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

eta5 ~ -0.160*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
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
eta4 ~~ 0.533*eta4
eta5 ~~ 0.642*eta5

# indicator residual variances (rel = 0.4)
x1 ~~ 1.50*x1;  x2 ~~ 1.50*x2;  x3 ~~ 1.50*x3
x4 ~~ 1.50*x4;  x5 ~~ 1.50*x5;  x6 ~~ 1.50*x6
x7 ~~ 1.50*x7;  x8 ~~ 1.50*x8;  x9 ~~ 1.50*x9
x10 ~~ 1.50*x10; x11 ~~ 1.50*x11; x12 ~~ 1.50*x12
x13 ~~ 1.50*x13; x14 ~~ 1.50*x14; x15 ~~ 1.50*x15
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

eta5 ~ -0.160*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
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
eta4 ~~ 0.533*eta4
eta5 ~~ 0.642*eta5

# indicator residual variances (rel = 0.6)
x1 ~~ 0.67*x1;  x2 ~~ 0.67*x2;  x3 ~~ 0.67*x3
x4 ~~ 0.67*x4;  x5 ~~ 0.67*x5;  x6 ~~ 0.67*x6
x7 ~~ 0.67*x7;  x8 ~~ 0.67*x8;  x9 ~~ 0.67*x9
x10 ~~ 0.67*x10; x11 ~~ 0.67*x11; x12 ~~ 0.67*x12
x13 ~~ 0.67*x13; x14 ~~ 0.67*x14; x15 ~~ 0.67*x15
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

eta5 ~ -0.160*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
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
eta4 ~~ 0.533*eta4
eta5 ~~ 0.642*eta5

# indicator residual variances (rel = 0.8)
x1 ~~ 0.25*x1;  x2 ~~ 0.25*x2;  x3 ~~ 0.25*x3
x4 ~~ 0.25*x4;  x5 ~~ 0.25*x5;  x6 ~~ 0.25*x6
x7 ~~ 0.25*x7;  x8 ~~ 0.25*x8;  x9 ~~ 0.25*x9
x10 ~~ 0.25*x10; x11 ~~ 0.25*x11; x12 ~~ 0.25*x12
x13 ~~ 0.25*x13; x14 ~~ 0.25*x14; x15 ~~ 0.25*x15
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

eta5 ~ -0.160*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
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
eta4 ~~ 0.533*eta4
eta5 ~~ 0.642*eta5

# indicator residual variances (rel = 0.4)
x1 ~~ 1.50*x1;  x2 ~~ 1.50*x2;  x3 ~~ 1.50*x3
x4 ~~ 1.50*x4;  x5 ~~ 1.50*x5;  x6 ~~ 1.50*x6
x7 ~~ 1.50*x7;  x8 ~~ 1.50*x8;  x9 ~~ 1.50*x9
x10 ~~ 1.50*x10; x11 ~~ 1.50*x11; x12 ~~ 1.50*x12
x13 ~~ 1.50*x13; x14 ~~ 1.50*x14; x15 ~~ 1.50*x15
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

eta5 ~ -0.160*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
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
eta4 ~~ 0.533*eta4
eta5 ~~ 0.642*eta5

# indicator residual variances (rel = 0.6)
x1 ~~ 0.67*x1;  x2 ~~ 0.67*x2;  x3 ~~ 0.67*x3
x4 ~~ 0.67*x4;  x5 ~~ 0.67*x5;  x6 ~~ 0.67*x6
x7 ~~ 0.67*x7;  x8 ~~ 0.67*x8;  x9 ~~ 0.67*x9
x10 ~~ 0.67*x10; x11 ~~ 0.67*x11; x12 ~~ 0.67*x12
x13 ~~ 0.67*x13; x14 ~~ 0.67*x14; x15 ~~ 0.67*x15
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

eta5 ~ -0.160*1 + 0.16*eta4 + 0.16*eta1 + 0.16*eta2 + 0.16*eta3
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
eta4 ~~ 0.533*eta4
eta5 ~~ 0.642*eta5

# indicator residual variances (rel = 0.8)
x1 ~~ 0.25*x1;  x2 ~~ 0.25*x2;  x3 ~~ 0.25*x3
x4 ~~ 0.25*x4;  x5 ~~ 0.25*x5;  x6 ~~ 0.25*x6
x7 ~~ 0.25*x7;  x8 ~~ 0.25*x8;  x9 ~~ 0.25*x9
x10 ~~ 0.25*x10; x11 ~~ 0.25*x11; x12 ~~ 0.25*x12
x13 ~~ 0.25*x13; x14 ~~ 0.25*x14; x15 ~~ 0.25*x15
"
