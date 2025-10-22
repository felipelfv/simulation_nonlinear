## README 

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://github.com/felipelfv/Simulation_First_PhD_Article">A simulation study comparing structural-after-measurement versus traditional approaches to estimate nonlinear effects in structural equation modeling</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://github.com/felipelfv">Felipe Fontana Vieira</a> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>

Repository for the first article of my PhD: A simulation study comparing 
structural-after-measurement versus traditional approaches 
to estimate nonlinear effects in structural equation modeling

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

Contains all simulation studies and supporting scripts.

```{r}
Simulations/
├── Study_1/                       # First simulation study
│   ├── Results/                   # Study 1 results
│   │   └── Results(1).R           # Results analysis script
│   └── Simulation/                # Study 1 simulation code
│       ├── Models(1).RData        # Saved model objects
│       └── Simulation(1).R        # Main simulation script
│
├── Study_2/                       # Second simulation study
│   ├── Results/                   # Study 2 results
│   │   └── Results(2).R           # Results analysis script
│   └── Simulation/                # Study 2 simulation code
│       ├── Models(2).RData        # Saved model objects
│       └── Simulation(2).R        # Main simulation script
│
├── GenerateData.R                 # Data generation function used for both studies
├── Methods.R                      # Common methods/approaches for both studies
└── Plots.R                        # Plotting functions for both results
```

### 3. Debug

This folder contains two scripts which I have used to check the data generation mechanism as well as retrieve a dataset generated under the simulation.
Hence, the interested reader may find these functions helpful. However, they are not relevant for the manuscript and results reported. 

1. `Data_Retrieval(debug).R` contains the R code and information needed to reproduce the datasets generated for any particular iterations across conditions for both simulation studies.
2. `GenerateData(debug).R` contains the R code and information needed to generate data given specific parameters. This is useful to check if the data is generated correctly (i.e., according to our expectation).

### 4. renv

Manages R package dependencies and ensures consistent package versions across environments. Configured via renv.lock file

## Reproducing the simulation 



## Reproducing the results (i.e., performance metrics and visualization)