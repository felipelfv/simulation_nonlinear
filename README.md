## README 

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://github.com/felipelfv/simulation_nonlinear.git">Evaluating Local Structural-After-Measurement (LSAM) and Traditional Approaches for the Estimation of Complex Nonlinear Effects Among Latent Variables</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://github.com/felipelfv">Felipe Fontana Vieira</a>, Kjell Solem Slupphaug, and Yves Rosseel is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>

Repository for the article: Evaluating Local Structural-After-Measurement (LSAM) and Traditional Approaches for the Estimation of Complex Nonlinear Effects Among Latent Variables

## Root Directory Structure

```{r}
Simulation_First_PhD_Article/
├── Debug/                          # Debugging scripts
├── Manuscript/                     # All manuscript and supplemental materials
├── Simulations/                    # Simulation code and analyses
├── renv/                           # R environment management
├── .gitattributes                  # .RData files handled by Git LFS
├── .gitignore                      # Git ignore rules
├── README.html                     # HTML version of README
├── README.md                       # Repository documentation
├── SimulationStudy_Nonlinear.Rproj # RStudio project file
└── renv.lock                       # R dependency lock file
```

### 1. Manuscript 

Contains all manuscript materials converted to Quarto format.

```{r}
Manuscript/
├── _extensions/wjschne/apaquarto/  # APA Quarto extension
├── references.bib                  # Bibliography file
├── appendix_a.qmd                  # Appendix A 
├── appendix_b.qmd                  # Appendix B 
├── appendix_c.qmd                  # Appendix C 
├── design.qmd                      # Simulation design section 
├── discussion.qmd                  # Discussion section 
├── introduction.qmd                # Introduction section
├── manuscript.qmd                  # Main manuscript file
├── modelsframework.qmd             # Models and framework section
├── results.qmd                     # Results section
└── supplemental.Rmd                # Supplemental materials (R markdown)
```

### 2. Simulations

Contains all simulation studies and supporting scripts. 

```{r}
Simulations/
├── Study_1/                             # First simulation study
│   ├── Data/                            # Study 1 raw simulation output (Git LFS)
│   │   └── Data_Study_1_final.RData     # All iterations data from Simulation(1).R
│   ├── Results/                         # Study 1 processed results (Git LFS)
│   │   ├── Results(1).R                 # Calculates performance metrics from Data/
│   │   ├── results_study_1.rds          # Output: CalculatePerformance() with warnings
│   │   └── results_study_1_excluding_warnings.rds  # Output: CalculatePerformance() excluding warnings
│   └── Simulation/                      # Study 1 simulation code
│       ├── Models(1).RData              # Model specifications
│       └── Simulation(1).R              # Runs simulation → saves to Data/
│
├── Study_2/                             # Second simulation study
│   ├── Data/                            # Study 2 raw simulation output (Git LFS)
│   │   └── Data_Study_2_final.RData     # All iterations data from Simulation(2).R
│   ├── Results/                         # Study 2 processed results (Git LFS)
│   │   ├── Results(2).R                 # Calculates performance metrics from Data/
│   │   ├── results_study_2.rds          # Output: CalculatePerformance() with warnings
│   │   └── results_study_2_excluding_warnings.rds  # Output: CalculatePerformance() excluding warnings
│   └── Simulation/                      # Study 2 simulation code
│       ├── Models(2).RData              # Model specifications
│       └── Simulation(2).R              # Runs simulation → saves to Data/
│
├── GenerateData.R                       # Data generation function used for both studies
├── Methods.R                            # Common methods/approaches for both studies
└── Plots.R                              # Plotting functions for both results
```

Workflow:

1. Simulation(X).R → Data_Study_X_final.RData (raw iterations)

2. Results(X).R → results_study_x.rds (processed metrics)

3. Manuscript loads results_study_x.rds for tables/figures

### 3. Debug

This folder contains two scripts which I have used to check the data generation mechanism as well as retrieve a dataset generated under the simulation conditions.
Hence, the interested reader may find these functions helpful. However, they are not relevant for the manuscript and results reported. 

1. `Data_Retrieval(debug).R` contains the R code and information needed to reproduce the datasets generated for any particular iteration across conditions for both simulation studies.
2. `GenerateData(debug).R` contains the R code and information needed to generate data given specific parameters. This is useful to check if the data is generated correctly (i.e., according to our expectation).

### 4. renv

Manages R package dependencies and ensures consistent package versions across environments. Configured via renv.lock file

## Cloning this repository 

Note that when you clone this GitHub repository with LFS files (i.e., Data_Study_X_final.RData files) without the Git LFS installed, you will get pointer files instead of the actual .RData files.

You will need to install Git LFS. See https://git-lfs.github.com/ for more information. Once the repository is cloned and LFS is installed, you will need to run 'git lfs pull'.

## Reproducing the simulation 

(...Ignoring dockerfiles and/or Nix/rix...)

To reproduce the .RData files included in the folders ../Simulations/Study_1/Data and ../Simulations/Study_2/Data you just need to run (given the correctly specified folder structure) the Simulation(1).R and Simulation(2).R, respectively.

Those latter mentioned files are dependent on three files: GenerateData.R, Methods.R and Models(X).R, where X = 1,2 for Simulation 1 and Simulation 2.

### Packages needed for this

| Package        | Version     | Citation                       |
|----------------|-------------|--------------------------------|
| lavaan         | 0.6.21.2400 | Rosseel (2012)                 |
| modsem         | 1.0.12      | Slupphaug et al. (2024)        |
| doParallel     | 1.0.17      | Corporation and Weston (2022)  |
| doRNG          | 1.8.6.2     | Gaujoux (2025)                 |
| covsim         | 1.1.0       | Gronneberg et al. (2022)       |
| rvinecopulib   | 0.7.3.1.0   | Nagler & Vatter (2025)         |

## Reproducing the results (i.e., convergence/outliers information, performance metrics and visualization)

(...Ignoring dockerfiles and/or Nix/rix...)

### Convergence, outliers statistics, and performance metrics

For the results from Simulation 1:

```{r}
results_study_1 <- CalculatePerformance(
  all_results,
  parameters_of_interest = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"),
  remove_outliers = TRUE,
  outlier_threshold = 3,
  alpha = 0.05,
  min_reps = 10,
  exclude_warnings = FALSE,  # TRUE to exclude iterations with warnings
  return_convergence_details = FALSE
)
```

For the results from Simulation 2:

```{r}

results_study_2 <- CalculatePerformance_Study2(
  all_results,
  parameters_of_interest = list(
    eta4 = c("eta1", "eta2", "eta3", "eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2"),
    eta5 = c("eta4", "eta1", "eta2", "eta3", "eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
  ),
  remove_outliers = TRUE,
  outlier_threshold = 3,
  alpha = 0.05,
  min_reps = 10,
  exclude_warnings = FALSE,
  return_convergence_details = FALSE
)
```

### Figures

In the results.qmd file under labels 'plots-first-sim-study' and 'plots-second-sim-study'. 

For figures associated with simulation 1:

```{r}
#| label: plots-first-sim-study

source("../Simulations/Plots.R")
plot_study_1_data <- prepare_study1_data(results_study_1)

# RELATIVE BIAS 
study_1_relative_bias <- plot_study_1_data$bias_relative
dat_relative_bias_400_full_study_1 <- study_1_relative_bias %>%
   filter(Model == "Full", SampleSize == "400") %>%
   droplevels()

p_relative_bias_400_full_study_1 <- plot_bias(
   dat_relative_bias_400_full_study_1, 
   shapes = SHAPES_4, 
   ltys = LTYS_4,
   bias_type = "relative",  
   facet_formula = Distribution + Parameter ~ .,  
   y_breaks = seq(-1, 1.5, by = 0.5),
   y_limits = c(-1, 1.5)
)

# Relative RMSE 
study_1_relative_rmse <- plot_study_1_data$rmse_relative

dat_relative_rmse_400_full_study_1 <- study_1_relative_rmse %>%
  filter(Model == "Full", SampleSize == "400") %>%
  droplevels()

p_relative_rmse_400_full_study_1 <- plot_rmse(
  dat_relative_rmse_400_full_study_1, 
  shapes = SHAPES_4, 
  ltys = LTYS_4,
  facet_formula = Distribution + Parameter ~ ., 
  #y_breaks = seq(0, 0.30, by = 0.05),           
  #y_limits = c(0, 0.30)                        
)

# SE/SD RATIO (Right-skewed & Uniform, Full, N=400 & 1000)
study_1_sesd <- plot_study_1_data$sesd

dat_sesd_nonnorm_unif_study_1 <- study_1_sesd %>%
  filter(Model == "Full", 
         SampleSize %in% c("400", "1000"), 
         Distribution %in% c("Right-skewed", "Uniform")) %>%
  droplevels()

p_slope_nonnorm_unif_study_1 <- plot_sesd(
  dat_sesd_nonnorm_unif_study_1,
  shapes = SHAPES_4,
  ltys = LTYS_4,
  facet_formula = Distribution + Parameter ~ SampleSize
)

# COVERAGE (Full, N=400 & 1000, all distributions)
study_1_coverage <- plot_study_1_data$coverage

dat_cov_nonnorm_unif_study_1 <- study_1_coverage %>%
  filter(Model == "Full", 
         SampleSize %in% c("400", "1000"),
         Distribution %in% c("Right-skewed", "Uniform")) %>% 
  droplevels()

p_cov_nonnorm_unif_study_1 <- plot_coverage(
  dat_cov_nonnorm_unif_study_1, 
  shapes = SHAPES_4,
  ltys = LTYS_4,
  facet_formula = Distribution + Parameter ~ SampleSize
)

# TYPE I ERROR (Linear, N=1000 only) 
study_1_type1 <- plot_study_1_data$type1

dat_t1_1000_study_1 <- study_1_type1 %>%
  filter(grepl("N=1000", Condition)) %>%
  droplevels()

p_t1_1000_study_1 <- plot_type1(
  dat_t1_1000_study_1,  
  greys = GREYS_4,
  facet_formula = Distribution + Parameter ~ Condition
)
```

For figures associated with simulation 2:

```{r}
#| label: plots-second-sim-study

source("../Simulations/Plots.R")
study_2_normal_data <- prepare_study2_data(results_study_2, "normal")
study_2_nonnormal_data <- prepare_study2_data(results_study_2, "nonnormal")
study_2_uniform_data <- prepare_study2_data(results_study_2, "uniform")

# RELATIVE BIAS (full model, N=400, all Distributions) 
dat_relative_bias_400_full_all_study_2 <- bind_rows(
  study_2_normal_data$bias_relative %>%
    filter(Model == "Full", SampleSize == "400"),
  study_2_nonnormal_data$bias_relative %>%
    filter(Model == "Full", SampleSize == "400"),
  study_2_uniform_data$bias_relative %>%
    filter(Model == "Full", SampleSize == "400")
) %>%
  droplevels()

p_relative_bias_400_full_all_study_2 <- plot_bias(
  dat_relative_bias_400_full_all_study_2,
  shapes = SHAPES_3,
  ltys = LTYS_3,
  facet_formula = Parameter ~ Distribution,
  y_breaks = seq(-1, 1.5, by = 0.5),
  y_limits = c(-1, 1.5),
  bias_type = "relative"
)

# RELATIVE RMSE (full model, N=400, all Distributions) 
dat_relative_rmse_400_full_all_study_2 <- bind_rows(
  study_2_normal_data$rmse_relative %>%
    filter(Model == "Full", SampleSize == "400"),
  study_2_nonnormal_data$rmse_relative %>%
    filter(Model == "Full", SampleSize == "400"),
  study_2_uniform_data$rmse_relative %>%
    filter(Model == "Full", SampleSize == "400")
) %>%
  droplevels()

p_relative_rmse_400_full_all_study_2 <- plot_rmse(
    dat_relative_rmse_400_full_all_study_2,
    shapes = SHAPES_3,
    ltys = LTYS_3,
    facet_formula = Parameter ~ Distribution,
    y_breaks = seq(0, 5, by = 1),  
    y_limits = c(0, 5.2),              
)


# SE/SD RATIO (Full model, N=400, All Distributions) 
dat_sesd_400_full_all_study_2 <- bind_rows(
  study_2_normal_data$sesd %>%
    filter(Model == "Full", SampleSize == "400"),
  study_2_nonnormal_data$sesd %>%
    filter(Model == "Full", SampleSize == "400"),
  study_2_uniform_data$sesd %>%
    filter(Model == "Full", SampleSize == "400")
) %>%
  droplevels()

p_sesd_400_full_all_study_2 <- plot_sesd(
  dat_sesd_400_full_all_study_2,
  shapes = SHAPES_3,
  ltys = LTYS_3,
  facet_formula = Parameter ~ Distribution,
  title = ""
)

#  COVERAGE (Full model, N=400, All Distributions)
dat_coverage_400_full_all_study_2 <- bind_rows(
  study_2_normal_data$coverage %>%
    filter(Model == "Full", SampleSize == "400"),
  study_2_nonnormal_data$coverage %>%
    filter(Model == "Full", SampleSize == "400"),
  study_2_uniform_data$coverage %>%
    filter(Model == "Full", SampleSize == "400")
) %>%
  droplevels()

p_coverage_400_full_all_study_2 <- plot_coverage(
  dat_coverage_400_full_all_study_2,
  shapes = SHAPES_3,
  ltys = LTYS_3,
  facet_formula = Parameter ~ Distribution,
  title = ""
)

# TYPE I ERROR (Linear, N=1000, Nonnormal & Uniform)
dat_t1_1000_all_study_2 <- bind_rows(
  study_2_nonnormal_data$type1 %>%
    filter(grepl("N=1000", Condition)),
  study_2_uniform_data$type1 %>%
    filter(grepl("N=1000", Condition))
) %>%
  droplevels()

p_t1_1000_all_study_2 <- plot_type1(
  dat_t1_1000_all_study_2,
  greys = GREYS_3,
  facet_formula = Parameter ~ Distribution + Condition,
  title = ""
)
```

### Packages needed for this

| Package        | Version     | Citation                         |
|----------------|-------------|----------------------------------|
| simhelpers     | 0.3.1       | Joshi & Pustejovsky (2025)       |
| dplyr          | 1.1.4       | Wickham et al. (2023)            |
| ggplot2        | 3.5.2       | Wickham (2025)                   |
| tibble         | 3.3.0       | Müller & Wickham (2025)          |
