##==============================================================================
## General Information 
##
## Script name: Models(1).R
##
## Purpose of script: The models defined here (lavaan-syntax) are used to 
## generate data (specified in GenerateData() through parameter `model`) as well 
## as to fit/estimate through the functions defined in script Methods(1).R. For 
## the former, we have the population.X.model(s), for the latter we have the 
## fit.X.model(s).
##
##==============================================================================

##==============================================================================
## 1. Population Models 
##==============================================================================

##### 1.1 First Population Model (M_{lin}) #####

# No interaction terms and only used to examine the Type I error rate
# \eta = 0.316*\xi_{1} + 0.316*\xi_{2} + \zeta

population.linear.model <- "
# Measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9

# Structural model
eta3 ~ 0.316*eta1 + 0.316*eta2

# (Co)variances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
"

##### 1.2 Second Population Model (M_{int}) #####

# Note: this was not used in the simulation study in the end
# One interaction included
# \eta = -0.255 + 0.316*\xi_{1} + 0.316*\xi_{2} + 0.139*\xi_{1}\xi_{2} + \zeta

population.interaction.model <- "
# Measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9

# Structural model
eta3 ~ -0.255*1 + 0.316*eta1 + 0.316*eta2 + 0.139*eta1:eta2

# (Co)variances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
"

##### 1.3 Third Population Model (M_{full}) #####

# Two quadratic effects (on top of the interaction)
# \eta = -0.255 + 0.316*\xi_{1} + 0.316*\xi_{2} + 0.139*\xi_{1}\xi_{2} + 
# 0.101*\xi_{1}^2 + 0.101*\xi_{2}^2 + \zeta

population.full.model <- "
# Measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9

# Structural model
eta3 ~ -0.255*1 + 0.316*eta1 + 0.316*eta2 + 0.139*eta1:eta2 + 0.101*eta1:eta1 + 0.101*eta2:eta2

# (Co)variances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
"

##==============================================================================
## 2. Fit Models 
##==============================================================================


##### 2.1 First Fitted Model (A_{int}) #####

# Note: this was not used in the simulation study in the end

fit.interaction.model <- "
# Measurement model
eta1 =~ x1 + x2 + x3
eta2 =~ x4 + x5 + x6
eta3 =~ x7 + x8 + x9

# Structural model
eta3 ~ eta1 + eta2 + eta1:eta2

# (Co)variances
eta1 ~~ eta2
eta1 ~~ eta1
eta2 ~~ eta2
"

##### 2.2 Second Fitted Model (A_{full}) #####

fit.full.model <- "
# Measurement model
eta1 =~ x1 + x2 + x3
eta2 =~ x4 + x5 + x6
eta3 =~ x7 + x8 + x9

# Structural model
eta3 ~ eta1 + eta2 + eta1:eta2 + eta1:eta1 + eta2:eta2

# (Co)variances
eta1 ~~ eta2
eta1 ~~ eta1
eta2 ~~ eta2
"
