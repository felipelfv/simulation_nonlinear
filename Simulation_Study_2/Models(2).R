#### 1. General Information ####

#################
#### Study 2 ####
#################

# For using modsem:
# See: https://modsem.org/articles/interaction_two_etas.html

##### 2.1 Population Large Model #####

###### 2.1.1 Model ######
if (FALSE) {
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
}

###### 2.1.2 Model ######
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

##### 3.2 Fit Large Model #####

###### 3.2.1 Model ######
# For the second population model (!)
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
# note: this part is redundant, but included for clarity
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


