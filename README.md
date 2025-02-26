## README 

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://github.com/felipelfv/Simulation_First_PhD_Article">A simulation study comparing structural-after-measurement versus traditional approaches to estimate nonlinear effects in structural equation modeling</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://github.com/felipelfv">Felipe Fontana Vieira</a> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>

Repository for the first article of my PhD: A simulation study comparing 
structural-after-measurement versus traditional approaches 
to estimate nonlinear effects in structural equation modeling

# Functions

## GenerateData.R

This script contains the function that will parse the model information (in lavaan syntax) 
and generate data based on some user-specified arguments.

```{r}
GenerateData <- function(model,
                         N = 1000L,
                         skewness = NULL,
                         excesskurtosis = NULL,
                         exo.mean = NULL,
                         distr.exo = "rIG",
                         distr.zeta = "normal",
                         distr.epsilon = "normal",
                         center.exogenous.latent = TRUE,
                         center.exogenous.manifest = TRUE,
                         center.lv.dependent = FALSE,
                         center.lv.prod = FALSE,
                         center.indicators = FALSE,
                         target.var = NULL,
                         R2 = NULL,
                         rel = 0.64,
                         seed = NULL,
                         add.eta = FALSE)
```

## Methods.R

This file contains the functions dedicated for the estimation based on the different approaches. 

```{r}
method_uca <- function(Data = NULL, model.fit = NULL)

method_sam <- function(Data = NULL, estimator = "ML",
                       joint = FALSE, add.attr = FALSE, 
                       model.fit = NULL)
                       
method_analytic <- function(Data = NULL, model.fit = NULL, 
                            standardized = FALSE, method = "lms")
```

## Models.R

This file contains all the models specified in lavaan syntax that are passed to GenerateData() in GenerateData.R

## Design.R

This file contains the parameters that we vary in the simulation study.

## Simulation.R

This file contains the actual script for running the simulation. Importantly, this script is dependent on the abovementioned. 