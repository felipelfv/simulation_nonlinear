#### 1. General Information ####

#### 2. Population Models - Under the Alternative ####

# Here we start with the population models generated under the alternative
# hypothesis 

##### 2.1 First Population Model (M_{lin}) #####
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

##### 2.2 Second Population Model (M_{int}) #####
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

##### 2.3 Third Population Model (M_{full}) #####
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

#### 3. Population Models - Under the Null ####

# Here we add the population models generated under the null hypothesis 

##### 3.1. Null Model for Main Effects #####
population.linear.model.null <- "
# Measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
# Structural model
eta3 ~ 0*eta1 + 0*eta2
# (Co)variances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
"

##### 3.2 Null Model for Interaction Effect #####
population.interaction.model.null <- "
# Measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
# Structural model
eta3 ~ -0.255*1 + 0.316*eta1 + 0.316*eta2 + 0*eta1:eta2
# (Co)variances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
"

##### 3.3. Null Model for all nonlinear effects #####
population.full.model.null <- "
# Measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
# Structural model
eta3 ~ -0.255*1 + 0.316*eta1 + 0.316*eta2 + 0*eta1:eta2 + 0*eta1:eta1 + 0*eta2:eta2
# (Co)variances
eta1 ~~ 0.375*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
"

#### 4. Fitted Models ####

##### 4.1 First Fitted Model (A_{int}) #####
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

##### 4.2 Second Fitted Model (A_{full}) #####
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

#### 5. Study 2 ####

##### 5.1 Population Large Model #####
population.large.model <- "
# measurement model 
eta1  =~ 1*x1 + 1*x2 + 1*x3 + 1*x4
eta2  =~ 1*x5 + 1*x6 + 1*x7
eta3  =~ 1*x8 + 1*x9 + 1*x10
eta4  =~ 1*x11 + 1*x12 + 1*x13 + 1*x14
eta5  =~ 1*x15 + 1*x16 + 1*x17

# structural part (regressions)
eta4 ~ 0.2*eta1 + 0.6*eta2 + 0.3*eta3 + 0.2*norm
             + 0.1*eta2:eta3 + 0.4*eta1:eta3 
             + 0.3*norm:eta3 + 1*eta5:eta2
eta5 ~ 2*eta1:eta2
eta3 ~ 0.7*gender + 0.3*age + 0.4*norm
norm ~ 0.2*gender + 0.4*age

# (co)variances
eta1 ~~ 1*eta1
eta1 ~~ 0.4*eta2
eta2 ~~ 1*eta2
gender ~~ 0.1*eta1
gender ~~ 0.2*eta2
gender ~~ 0.2*age
gender ~~ 1*gender
age ~~ 1*age
age ~~ 0.2*eta1
age ~~ 0.3*eta2
"

#R2 = list(
  #"norm" = 0.4,
  #"eta3" = 0.5,
  #"eta4" = 0.6,
  #"eta5" = 0.5
#)

population.large.model <- "
# measurement model 
eta1  =~ 1*x1 + 1*x2 + 1*x3 + 1*x4
eta2  =~ 1*x5 + 1*x6 + 1*x7
eta3  =~ 1*x8 + 1*x9 + 1*x10
eta4  =~ 1*x11 + 1*x12 + 1*x13 + 1*x14
eta5  =~ 1*x15 + 1*x16 + 1*x17
eta6  =~ 1*x18 + 1*x19 + 1*x20
  
# structural part (regressions)
eta4 ~ 0.2*eta1 + 0.6*eta2 + 0.3*eta3 + 
       0.1*eta2:eta3 + 0.4*eta1:eta3 + 1*eta6:eta2 +
       0.3*eta5:eta6 + 0.2*eta5:eta5 + 0.1*eta1:eta2

# (co)variances
eta1 ~~ 1*eta1
eta1 ~~ 0.4*eta2
eta1 ~~ 0.2*eta3
eta1 ~~ 0.1*eta5
eta1 ~~ 0.1*eta6
eta2 ~~ 0.1*eta3
eta2 ~~ 0.3*eta5
eta2 ~~ 0.4*eta6
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3
eta3 ~~ 0.2*eta5
eta3 ~~ 0.2*eta6
eta5 ~~ 1*eta5
eta5 ~~ 0.2*eta6
eta6 ~~ 1*eta6
"

fit.large.model <- "
# measurement model 
eta1  =~ x1 + x2 + x3 + x4
eta2  =~ x5 + x6 + x7
eta3  =~ x8 + x9 + x10
eta4  =~ x11 + x12 + x13 + x14
eta5  =~ x15 + x16 + x17
eta6  =~ x18 + x19 + x20
  
# structural part (regressions)
eta4 ~ eta1 + eta2 + eta3 + 
       eta2:eta3 + eta1:eta3 + eta6:eta2 +
       eta5:eta6 + eta5:eta5 + eta1:eta2

# (co)variances
eta1 ~~ eta1
eta1 ~~ eta2
eta1 ~~ eta3
eta1 ~~ eta5
eta1 ~~ eta6
eta2 ~~ eta3
eta2 ~~ eta5
eta2 ~~ eta6
eta2 ~~ eta2
eta3 ~~ eta3
eta3 ~~ eta5
eta3 ~~ eta6
eta5 ~~ eta5
eta5 ~~ eta6
eta6 ~~ eta6
"

##### 5.2 Fit Large Model #####
# See: https://modsem.org/articles/interaction_two_etas.html

