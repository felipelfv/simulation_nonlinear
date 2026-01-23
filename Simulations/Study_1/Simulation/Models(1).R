# population models: based on normal distribution for all variables
# that is latent exogenous, epsilon and zeta
all_models <- list()

# full models with the population parameters 
all_models[["normal_rel04"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9

# structural model
eta3 ~ -0.255*1
     + 0.316*eta1
     + 0.316*eta2
     + 0.139*eta1:eta2
     + 0.101*eta1:eta1
     + 0.101*eta2:eta2

# latent exogenous variables covariances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2

# structural residual variance (calibrated for R^2 approx. 0.30)
eta3 ~~ 0.756*eta3

# indicator residual variances (rel = 0.4)
x1 ~~ 1.5*x1
x2 ~~ 1.5*x2
x3 ~~ 1.5*x3
x4 ~~ 1.5*x4
x5 ~~ 1.5*x5
x6 ~~ 1.5*x6
x7 ~~ 1.893*x7
x8 ~~ 1.893*x8
x9 ~~ 1.893*x9
"

all_models[["normal_rel06"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9

# structural model
eta3 ~ -0.255*1
     + 0.316*eta1
     + 0.316*eta2
     + 0.139*eta1:eta2
     + 0.101*eta1:eta1
     + 0.101*eta2:eta2

# latent exogenous variables covariances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2

# structural residual variance (calibrated for R^2 approx. 0.30)
eta3 ~~ 0.756*eta3

# indicator residual variances (rel = 0.6)
x1 ~~ 0.667*x1
x2 ~~ 0.667*x2
x3 ~~ 0.667*x3
x4 ~~ 0.667*x4
x5 ~~ 0.667*x5
x6 ~~ 0.667*x6
x7 ~~ 0.845*x7
x8 ~~ 0.845*x8
x9 ~~ 0.845*x9
"

all_models[["normal_rel08"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9

# structural model
eta3 ~ -0.255*1
     + 0.316*eta1
     + 0.316*eta2
     + 0.139*eta1:eta2
     + 0.101*eta1:eta1
     + 0.101*eta2:eta2

# latent exogenous variables covariances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2

# structural residual variance (calibrated for R^2 approx. 0.30)
eta3 ~~ 0.756*eta3

# indicator residual variances (rel = 0.8)
x1 ~~ 0.25*x1
x2 ~~ 0.25*x2
x3 ~~ 0.25*x3
x4 ~~ 0.25*x4
x5 ~~ 0.25*x5
x6 ~~ 0.25*x6
x7 ~~ 0.316*x7
x8 ~~ 0.316*x8
x9 ~~ 0.316*x9
"

# null models (no interaction and quadratic effects)
all_models[["null_normal_rel04"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9

# structural model
eta3 ~ -0.255*1
     + 0.316*eta1
     + 0.316*eta2
     + 0*eta1:eta2
     + 0*eta1:eta1
     + 0*eta2:eta2

# latent exogenous variables covariances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2

# structural residual variance (same R^2 baseline)
eta3 ~~ 0.756*eta3

# indicator residual variances (rel = 0.4)
x1 ~~ 1.5*x1
x2 ~~ 1.5*x2
x3 ~~ 1.5*x3
x4 ~~ 1.5*x4
x5 ~~ 1.5*x5
x6 ~~ 1.5*x6
x7 ~~ 1.893*x7
x8 ~~ 1.893*x8
x9 ~~ 1.893*x9
"

all_models[["null_normal_rel06"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9

# structural model 
eta3 ~ -0.255*1
     + 0.316*eta1
     + 0.316*eta2
     + 0*eta1:eta2
     + 0*eta1:eta1
     + 0*eta2:eta2

# latent exogenous variables covariances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2

# structural residual variance (same R^2 baseline)
eta3 ~~ 0.756*eta3

# indicator residual variances (rel = 0.6)
x1 ~~ 0.667*x1
x2 ~~ 0.667*x2
x3 ~~ 0.667*x3
x4 ~~ 0.667*x4
x5 ~~ 0.667*x5
x6 ~~ 0.667*x6
x7 ~~ 0.845*x7
x8 ~~ 0.845*x8
x9 ~~ 0.845*x9
"

all_models[["null_normal_rel08"]] <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9

# structural model 
eta3 ~ -0.255*1
     + 0.316*eta1
     + 0.316*eta2
     + 0*eta1:eta2
     + 0*eta1:eta1
     + 0*eta2:eta2

# latent exogenous variables covariances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2

# structural residual variance (same R^2 baseline)
eta3 ~~ 0.756*eta3

# indicator residual variances (rel = 0.8)
x1 ~~ 0.25*x1
x2 ~~ 0.25*x2
x3 ~~ 0.25*x3
x4 ~~ 0.25*x4
x5 ~~ 0.25*x5
x6 ~~ 0.25*x6
x7 ~~ 0.316*x7
x8 ~~ 0.316*x8
x9 ~~ 0.316*x9
"
