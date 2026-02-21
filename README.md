## README

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://github.com/felipelfv/simulation_nonlinear.git">Evaluating Local Structural-After-Measurement (LSAM) and Traditional Approaches for the Estimation of Complex Nonlinear Effects Among Latent Variables</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://github.com/felipelfv">Felipe Fontana Vieira</a>, Kjell Solem Slupphaug, and Yves Rosseel is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>

Repository for the article: Evaluating Local Structural-After-Measurement (LSAM) and Traditional Approaches for the Estimation of Complex Nonlinear Effects Among Latent Variables.

## Repository Structure

```
simulation_nonlinear/
├── Manuscript/                     # Manuscript and supplemental materials (.qmd)
├── Simulations/                    # Simulation code, data, and results
├── Debug/                          # Debugging/data-retrieval helper scripts
├── renv/                           # R environment management
├── SimulationStudy_Nonlinear.Rproj # RStudio project file
└── renv.lock                       # R dependency lock file
```

### 1. Simulations

Contains all simulation code and output for both studies. Each study follows the same structure.

```
Simulations/
├── GenerateData.R                       # Data generation functions (both studies)
├── Methods.R                            # Estimation methods: LMS, QML, UPI, LSAM
├── Plots.R                              # Plotting functions for results
│
├── Study_1/
│   ├── Simulation/
│   │   ├── Models(1).R                  # lavaan model specifications
│   │   └── Simulation(1).R              # Runs 1000 replications per condition
│   ├── Data/                            # Per-condition .RData files (Git LFS)
│   └── Results/
│       ├── Results(1).R                 # Computes performance metrics from Data/
│       ├── results_study_1_basic.rds
│       ├── results_study_1_basic_robust.rds
│       ├── results_study_1_strict.rds
│       └── results_study_1_strict_robust.rds
│
└── Study_2/
    ├── Simulation/
    │   ├── Models(2).R                  # lavaan model specifications
    │   └── Simulation(2).R              # Runs 1000 replications per condition
    ├── Data/                            # Per-condition .RData files (Git LFS)
    └── Results/
        ├── Results(2).R                 # Computes performance metrics from Data/
        ├── results_study_2_basic.rds
        ├── results_study_2_basic_robust.rds
        ├── results_study_2_strict.rds
        └── results_study_2_strict_robust.rds
```

**Results variants:** Each study produces four result files based on the combination of convergence checking and standard error type:
- `basic` / `strict` — basic convergence checks (NA/NaN/Inf only) vs. strict checks (positive SEs, positive variances, positive definite factor covariance matrix)
- `robust` suffix — results computed with robust standard errors

### 2. Manuscript

Contains all manuscript materials in Quarto format, rendered with the `apaquarto` extension.

```
Manuscript/
├── manuscript.qmd                  # Main manuscript file
├── introduction.qmd                # Introduction
├── modelsframework.qmd             # Models and framework
├── present_research.qmd            # Present research
├── design.qmd                      # Simulation design
├── results.qmd                     # Results (loads .rds files, generates figures)
├── discussion.qmd                  # Discussion
├── appendix_a.qmd                  # Appendix A
├── appendix_b.qmd                  # Appendix B
├── appendix_c.qmd                  # Appendix C
├── supplemental.qmd                # Supplemental materials
├── supplemental_tables.qmd         # Supplemental tables
├── references.bib                  # Bibliography
└── _extensions/wjschne/apaquarto/  # APA Quarto extension
```

### 3. Debug

Helper scripts for inspecting the data generation mechanism. Not part of the main analysis pipeline.

- `Data_Retrieval(debug).R` — reproduce the dataset for any specific iteration/condition
- `GenerateData(debug).R` — generate data given specific parameters to verify correctness

## Workflow

```
1. Simulation(X).R  →  Data/  (per-condition .RData files with raw iteration output)
       ↓ depends on: GenerateData.R, Methods.R, Models(X).R

2. Results(X).R     →  results_study_X_*.rds  (performance metrics)

3. Manuscript (results.qmd)  loads .rds files  →  tables and figures
       ↓ depends on: Plots.R
```

## Cloning this repository

This repository uses Git LFS for `.RData` files. If you clone without Git LFS installed, you will get pointer files instead of the actual data.

Install Git LFS (see https://git-lfs.github.com/), then run:

```bash
git lfs pull
```

## Reproducing the simulation

To reproduce the `.RData` files in `Simulations/Study_X/Data/`, run `Simulation(1).R` and `Simulation(2).R` respectively. Each simulation script depends on `GenerateData.R`, `Methods.R`, and `Models(X).R`.

### Packages needed

| Package        | Version     | Citation                       |
|----------------|-------------|--------------------------------|
| lavaan         | 0.6.21.2400 | Rosseel (2012)                 |
| modsem         | 1.0.12      | Slupphaug et al. (2024)        |
| doParallel     | 1.0.17      | Corporation and Weston (2022)  |
| doRNG          | 1.8.6.2     | Gaujoux (2025)                 |
| covsim         | 1.1.0       | Gronneberg et al. (2022)       |
| rvinecopulib   | 0.7.3.1.0   | Nagler & Vatter (2025)         |

## Reproducing the results

Run `Results(1).R` and `Results(2).R` to compute performance metrics (convergence rates, bias, RMSE, SE/SD ratio, coverage, Type I error) from the raw simulation data. These scripts produce the `.rds` files loaded by the manuscript.

### Packages needed

| Package        | Version     | Citation                         |
|----------------|-------------|----------------------------------|
| simhelpers     | 0.3.1       | Joshi & Pustejovsky (2025)       |
| dplyr          | 1.1.4       | Wickham et al. (2023)            |
| ggplot2        | 3.5.2       | Wickham (2025)                   |
