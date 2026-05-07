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

**Results variants:** Each study produces four result files based on the combination of convergence checking:
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

The project follows a sequential pipeline. Each step produces output files that are consumed by the next step. The same structure applies to both Simulation 1 (simple 3-factor model) and Simulation 2 (complex 5-factor model).

### Step 1: Run the simulation

`Simulation(1).R` and `Simulation(2).R` are the main entry points. Each script:

1. Sources `Models(1).R` or `Models(2).R` to load the population model specifications (lavaan syntax with fixed parameter values for each reliability level and model type).
2. Sources `Methods.R`, which defines the estimation functions for each approach (UPI, LMS, QML, LSAM).
3. Sources `GenerateData.R`, which defines the VITA-based data generation function.
4. Sets up a parallel cluster with `doParallel` and runs 1,000 replications per condition using `%dorng%` for reproducible parallel random number streams.
5. Saves one `.RData` file per condition to `Simulations/Study_1/Data/` or `Simulations/Study_2/Data/`. Each `.RData` file contains the raw parameter tables, timing, warnings, and observed metrics for all replications within that condition.
6. After all conditions complete, reloads the per-condition files and combines them into a single `Data_Study_1_final.RData` or `Data_Study_2_final.RData`.

```
Models(1).R ──┐
Methods.R ────┼──> Simulation(1).R ──> Data/Data_Study_1_*.RData
GenerateData.R┘

Models(2).R ──┐
Methods.R ────┼──> Simulation(2).R ──> Data/Data_Study_2_*.RData
GenerateData.R┘
```

File dependency detail:
- `GenerateData.R` uses `lavaan` (to parse model syntax and extract parameters), `covsim` and `rvinecopulib` (to generate non-normal multivariate data via vine copulas).
- `Methods.R` uses `modsem` (for UPI, LMS, QML estimation) and `lavaan` (for LSAM estimation via `sam()`).
- `Models(1).R` and `Models(2).R` are plain R scripts that define lists of lavaan model strings. They have no package dependencies.

### Step 2: Process results

`Results(1).R` and `Results(2).R` each define two main functions:

1. `ExtractConvergenceOutliers()` — loads the `.RData` files, applies convergence checks (basic or strict) and outlier detection (IQR-based), and returns the filtered data.
2. `CalculatePerformanceMetrics()` — computes bias, RMSE, relative bias, SE/SD ratio, coverage, Type I error, power, and their MCSEs from the filtered data using `simhelpers`.

Running the commented-out code at the bottom of each results script produces four `.rds` files per study:

| File | Description |
|------|-------------|
| `results_study_1_basic.rds` | Basic convergence checks, standard SEs |
| `results_study_1_basic_robust.rds` | Basic convergence checks, robust SEs |
| `results_study_1_strict.rds` | Strict convergence checks, standard SEs |
| `results_study_1_strict_robust.rds` | Strict convergence checks, robust SEs |

The same naming pattern applies to Simulation 2.

```
Data/Data_Study_1_final.RData ──> Results(1).R ──> Results/results_study_1_*.rds
Data/Data_Study_2_final.RData ──> Results(2).R ──> Results/results_study_2_*.rds
```

### Step 3: Render the manuscript

The manuscript is written in Quarto (`.qmd`) using the `apaquarto` extension. The main file is `manuscript.qmd`, which includes child documents for each section.

- `results.qmd` loads the `.rds` files produced in Step 2, sources `Plots.R` for the plotting functions, and generates all figures and tables inline.
- `supplemental.qmd` and `supplemental_tables.qmd` similarly load `.rds` files for supplemental analyses.
- `design.qmd` retrieves package versions dynamically for the computational details section.

```
Results/results_study_1_*.rds ──┐
Results/results_study_2_*.rds ──┼──> results.qmd (sources Plots.R) ──> manuscript.qmd
Plots.R ────────────────────────┘
```

### Full pipeline summary

```
Models(1).R ─┐                                                  ┌─> results.qmd ─┐
Methods.R ───┼─> Simulation(1).R ─> .RData ─> Results(1).R ─> .rds ─┤               ├─> manuscript.qmd
GenerateData.R┘                                                 └─> supplemental.qmd
                                                                        ↑
Models(2).R ─┐                                                  ┌─> results.qmd ─┘
Methods.R ───┼─> Simulation(2).R ─> .RData ─> Results(2).R ─> .rds ─┤
GenerateData.R┘                                                 └─> supplemental.qmd
```

## Cloning this repository

This repository uses Git LFS for the Study 1 `.RData` files. If you clone without Git LFS installed, you will get pointer files instead of the actual data.

Install Git LFS (see https://git-lfs.github.com/), then run:

```bash
git lfs pull
```

The Study 2 `.RData` files are too large for Git LFS and are archived on Zenodo: <https://doi.org/10.5281/zenodo.20066656>. Download them and place them in `Simulations/Study_2/Data/`. The Study 1 `.RData` files are also available on Zenodo at the same DOI if you prefer not to use Git LFS.

## Reproducing the simulation

To reproduce the `.RData` files, run `Simulation(1).R` and `Simulation(2).R` from the project root directory. Each simulation script sources `GenerateData.R`, `Methods.R`, and the corresponding `Models` file. The scripts expect to be run with the working directory set to the repository root (i.e., where `SimulationStudy_Nonlinear.Rproj` is located).

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

Run `Results(1).R` and `Results(2).R` to compute performance metrics (convergence rates, bias, RMSE, SE/SD ratio, coverage, Type I error, power) from the `.RData` files. These scripts produce the `.rds` files loaded by the manuscript. The commented-out code at the bottom of each file shows how to run the extraction and metric computation for each variant (basic/strict, standard/robust SEs).

### Packages needed

| Package        | Version     | Citation                         |
|----------------|-------------|----------------------------------|
| simhelpers     | 0.3.1       | Joshi & Pustejovsky (2025)       |
| dplyr          | 1.1.4       | Wickham et al. (2023)            |
| ggplot2        | 3.5.2       | Wickham (2025)                   |
