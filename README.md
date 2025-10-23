## README 

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://github.com/felipelfv/simulation_nonlinear.git">Local Structural-After-Measurement (LSAM) and Traditional Approaches for Nonlinear Effects among Latent Variables: A Simulation Study</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://github.com/felipelfv">Felipe Fontana Vieira</a> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>

Repository for the article: Local Structural-After-Measurement (LSAM) and Traditional Approaches for Nonlinear Effects among Latent Variables: A Simulation Study

## Root Directory Structure

```{r}
Simulation_First_PhD_Article/
├── Debug/                          # Debugging scripts
├── Manuscript/                     # All manuscript and supplementary materials
├── Simulations/                    # Simulation code and studies
├── renv/                           # R environment management
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
├── appendix_files/figure-html/     # (Some) Appendix figure files
├── design_files/libs/              # Design support files
├── figure/                         # Main manuscript figures
├── Phd - Initial Readings.bib      # Bibliography file
├── appendix_a.qmd                  # Appendix A (Quarto markdown)
├── appendix_b.qmd                  # Appendix B (Quarto markdown)
├── design.qmd                      # Design section
├── discussion.qmd                  # Discussion section
├── introduction.qmd                # Introduction section
├── manuscript.qmd                  # Main manuscript file
├── modelsframework.qmd             # Models and framework section
├── results.qmd                     # Results section
└── supplemental.Rmd                # Supplemental materials (R markdown)
```

### 2. Simulations

Contains all simulation studies and supporting scripts. Note that the usage of "results" may sound confusing. I realized this after all the simulations were done.
Inside the folder "Data", results refer to the simulation runs, meaning estimates obtained for each iteration throughout conditions. Inside the folder "Results", results in this case refer to the performance metrics (i.e., bias, rmse, coverage, se/sd, etc).

```{r}
Simulations/
├── Study_1/                             # First simulation study
│   ├── Data/                            # Study 1 raw simulation output (Git LFS)
│   │   └── Results_Study_1_final.RData  # All iterations from Simulation(1).R
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
│   │   └── Results_Study_2_final.RData  # All iterations from Simulation(2).R
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

1. Simulation(X).R → Results_Study_X_final.RData (raw iterations)

2. Results(X).R → results_study_x.rds (processed metrics)

3. Manuscript loads results_study_x.rds for tables/figures

### 3. Debug

This folder contains two scripts which I have used to check the data generation mechanism as well as retrieve a dataset generated under the simulation conditions.
Hence, the interested reader may find these functions helpful. However, they are not relevant for the manuscript and results reported. 

1. `Data_Retrieval(debug).R` contains the R code and information needed to reproduce the datasets generated for any particular iteration across conditions for both simulation studies.
2. `GenerateData(debug).R` contains the R code and information needed to generate data given specific parameters. This is useful to check if the data is generated correctly (i.e., according to our expectation).

### 4. renv

Manages R package dependencies and ensures consistent package versions across environments. Configured via renv.lock file

## Reproducing the simulation 

To reproduce the .RData files included in the folders ../Simulations/Study_1/Data and ../Simulations/Study_2/Data you just need to run (given the correctly specified folder structure) the Simulation(1).R and Simulation(2).R, respectively.

Those latter mentioned files are dependent on three files: GenerateData.R, Methods.R and Models(X).R, where X = 1,2 for simulation 1 and simulation 2.

### Packages needed for this

| Package        | Version     | Citation                                                                                      |
|----------------|-------------|-----------------------------------------------------------------------------------------------|
| lavaan         |             |                                                                                               |
| modsem         |             |                                                                                               |
| doParallel     |             |                                                                                               |
| doRNG          |             |                                                                                               |
| covsim         |             |                                                                                               |
| rvinecopulib   |             |                                                                                               |

## Reproducing the results (i.e., performance metrics and visualization)



