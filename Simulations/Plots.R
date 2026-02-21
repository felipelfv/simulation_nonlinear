############################ 1. General Information ############################

# This file contains all the code for data wrangling and plots based on the 
# simulation results (with various performance metrics):
# Consistent "APA-style" formatting across all plots
# Data preparation functions for the different simulation designs
# Plotting functions for bias, RMSE, SE/SD ratio, coverage, Type I error, and 
# power
#
#' @note Dependencies:
#'   Required packages: ggplot2 and dplyr
#'   
# Structure:
#   2. Configuration     - Method order, visual specifications
#   3. Shared Components - Theme, helper functions
#   4. Data Preparation
#      4.1 Simulation 1  - Simple model (3 nonlinear terms), normal errors
#      4.2 Simulation 2  - Complex model (8 nonlinear terms), normal errors
#      4.3 Additional Distributional Conditions - Non-normal errors (dumbbell plots)
#   5. Plotting Functions
#      5.1 Main Plots    - Bias, RMSE, SE/SD, coverage, Type I, power
#      5.2 Additional Distributional Conditions - Dumbbell plots

############################### 2. Configuration ###############################

# METHOD SPECIFICATIONS 
METHOD_ORDER_4 <- c("LSAM", "LMS", "QML", "UPI")
METHOD_ORDER_3 <- c("LSAM", "QML", "UPI")

# DISTRIBUTION LABELS
DIST_LABS <- c(normal = "Normal", nonnormal = "Right-skewed", uniform = "Uniform")

# VISUAL SPECIFICATIONS 
SHAPES_4 <- c(LSAM = 16, LMS = 17, QML = 15, UPI = 18)
SHAPES_3 <- c(LSAM = 16, QML = 15, UPI = 18)

LTYS_4 <- c(LSAM = "solid", LMS = "dashed", QML = "dotdash", UPI = "twodash")
LTYS_3 <- c(LSAM = "solid", QML = "dotdash", UPI = "twodash")

GREYS_4 <- c(LSAM = "grey20", LMS = "grey45", QML = "grey65", UPI = "grey85")
GREYS_3 <- c(LSAM = "grey20", QML = "grey50", UPI = "grey80")

############################### 3. Shared Components ###########################

#' APA-style theme for ggplot2
#' @param base_size Base font size. Default is 11.
#' @return ggplot2 theme object
theme_apa_bw <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

#' Create condition labels from sample size and reliability
#' @param N Sample size values
#' @param Rel Reliability values
#' @return Factor with ordered condition labels
make_condition <- function(N, Rel) {
  factor(paste0("N=", N, ", Rel=", Rel),
         levels = c("N=400, Rel=0.4", "N=400, Rel=0.6", "N=400, Rel=0.8",
                    "N=1000, Rel=0.4", "N=1000, Rel=0.6", "N=1000, Rel=0.8"))
}

############################### 4. Data Preparation ############################

# 4.1 STUDY 1

#' Prepare Simulation 1 data; normal errors only
#' @param results_data Full results data frame from CalculatePerformanceMetrics
#' @return Formatted data frame for Simulation 1
prep_study1 <- function(results_data) {
  
  PARAM_ORDER <- c("eta1:eta2", "eta1:eta1", "eta2:eta2")
  PARAM_LABS <- c("eta1:eta2" = "beta[33]",
                  "eta1:eta1" = "beta[34]",
                  "eta2:eta2" = "beta[35]")
  
  results_data %>%
    filter(
      Parameter %in% PARAM_ORDER,
      Distr_Epsilon == "normal",
      Distr_Zeta == "normal"
    ) %>%
    mutate(
      Method = factor(Method, levels = METHOD_ORDER_4),
      Distribution = factor(Distr_Exo, levels = names(DIST_LABS), labels = DIST_LABS),
      Parameter = factor(Parameter, levels = PARAM_ORDER, labels = PARAM_LABS),
      Reliability = factor(Reliability),
      SampleSize = factor(SampleSize),
      Condition = make_condition(SampleSize, Reliability),
      Model = factor(Model, levels = c("Linear", "Full"))
    )
}

#' Prepare all data frames for Simulation 1 plotting
#' @param results_data Raw simulation results data frame
#' @return List of formatted data frames for each metric
prepare_study1_data <- function(results_data) {
  
  data_prep <- prep_study1(results_data)
  
  list(
    bias = data_prep %>%
      transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
                y = Bias_Mean, yerr = Bias_Mean_MCSE),
    
    bias_relative = data_prep %>%
      transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
                y = RelativeBias_Mean, yerr = RelativeBias_MCSE),
    
    rmse = data_prep %>%
      transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
                y = RMSE_Mean, yerr = RMSE_Mean_MCSE),
    
    rmse_relative = data_prep %>%
      transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
                y = Relative_RMSE, yerr = Relative_RMSE_MCSE),
    
    sesd = data_prep %>%
      transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
                ratio = SE_SD_Ratio),
    
    coverage = data_prep %>%
      transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
                y = CoverageRate, yerr = CoverageRate_MCSE),
    
    type1 = data_prep %>%
      filter(Model == "Linear") %>%
      transmute(Distribution, Parameter, Condition, Method, 
                y = TypeI_Error, yerr = TypeI_Error_MCSE),
    
    power = data_prep %>%
      filter(Model == "Full") %>%
      transmute(Distribution, Parameter, Condition, Method, 
                y = Power, yerr = Power_MCSE)
  )
}

# 4.2 STUDY 2

#' Prepare Simulation 2 subset by distribution and model
#' @param df Data frame with simulation results
#' @param distribution Distribution name ("normal", "nonnormal", or "uniform")
#' @param model Either "Full" or "Linear"
#' @return Formatted data frame
prep_study2_subset <- function(df, distribution, model = c("Full", "Linear")) {
  model <- match.arg(model)
  
  eta4_params <- c("eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2")
  eta5_params <- c("eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
  
  df %>%
    filter(
      Distr_Exo == distribution,
      Distr_Epsilon == "normal",
      Distr_Zeta == "normal",
      Model == model,
      (Equation == "eta4" & Parameter %in% eta4_params) |
        (Equation == "eta5" & Parameter %in% eta5_params)
    ) %>%
    mutate(
      Distribution = factor(Distr_Exo, levels = names(DIST_LABS), labels = DIST_LABS),
      Parameter_Label = case_when(
        Equation == "eta4" & Parameter == "eta1:eta2" ~ "beta[44]",
        Equation == "eta4" & Parameter == "eta1:eta3" ~ "beta[45]",
        Equation == "eta4" & Parameter == "eta1:eta1" ~ "beta[46]",
        Equation == "eta4" & Parameter == "eta2:eta2" ~ "beta[47]",
        Equation == "eta5" & Parameter == "eta1:eta4" ~ "beta[55]",
        Equation == "eta5" & Parameter == "eta2:eta4" ~ "beta[56]",
        Equation == "eta5" & Parameter == "eta1:eta1" ~ "beta[57]",
        Equation == "eta5" & Parameter == "eta3:eta3" ~ "beta[58]",
        TRUE ~ Parameter
      ),
      Method = factor(Method, levels = METHOD_ORDER_3),
      Reliability = factor(Reliability),
      SampleSize = factor(SampleSize),
      Condition = make_condition(SampleSize, Reliability),
      Parameter_Label = factor(Parameter_Label, 
                               levels = c("beta[44]", "beta[45]", "beta[46]", "beta[47]",
                                          "beta[55]", "beta[56]", "beta[57]", "beta[58]"))
    )
}

#' Prepare all data frames for Simulation 2 plotting (per distribution)
#' @param results_data Raw simulation results data frame
#' @param dist_name Distribution name ("normal", "nonnormal", or "uniform")
#' @return List of formatted data frames for each metric
prepare_study2_data <- function(results_data, dist_name) {
  
  data_combined <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")))
  
  eta4_params <- c("eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2")
  eta5_params <- c("eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
  
  list(
    bias = data_combined %>%
      transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
                Method, Model, y = Bias_Mean, yerr = Bias_Mean_MCSE),
    
    bias_relative = data_combined %>%
      transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
                Method, Model, y = RelativeBias_Mean, yerr = RelativeBias_MCSE),
    
    rmse = data_combined %>%
      transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
                Method, Model, y = RMSE_Mean, yerr = RMSE_Mean_MCSE),
    
    rmse_relative = data_combined %>%
      transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
                Method, Model, y = Relative_RMSE, yerr = Relative_RMSE_MCSE),
    
    sesd = data_combined %>%
      transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
                Method, Model, ratio = SE_SD_Ratio),
    
    coverage = data_combined %>%
      transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
                Method, Model, y = CoverageRate, yerr = CoverageRate_MCSE),
    
    type1 = results_data %>%
      filter(
        Model == "Linear",
        Distr_Exo == dist_name,
        Distr_Epsilon == "normal",
        Distr_Zeta == "normal",
        (Equation == "eta4" & Parameter %in% eta4_params) |
          (Equation == "eta5" & Parameter %in% eta5_params)
      ) %>%
      mutate(
        Distribution = DIST_LABS[dist_name],
        Parameter = case_when(
          Equation == "eta4" & Parameter == "eta1:eta2" ~ "beta[44]",
          Equation == "eta4" & Parameter == "eta1:eta3" ~ "beta[45]",
          Equation == "eta4" & Parameter == "eta1:eta1" ~ "beta[46]",
          Equation == "eta4" & Parameter == "eta2:eta2" ~ "beta[47]",
          Equation == "eta5" & Parameter == "eta1:eta4" ~ "beta[55]",
          Equation == "eta5" & Parameter == "eta2:eta4" ~ "beta[56]",
          Equation == "eta5" & Parameter == "eta1:eta1" ~ "beta[57]",
          Equation == "eta5" & Parameter == "eta3:eta3" ~ "beta[58]"
        ),
        Method = factor(Method, levels = METHOD_ORDER_3),
        Condition = make_condition(SampleSize, Reliability),
        y = TypeI_Error,
        yerr = TypeI_Error_MCSE
      ) %>%
      select(Distribution, Parameter, Condition, Method, y, yerr),
    
    power = results_data %>%
      filter(
        Model == "Full",
        Distr_Exo == dist_name,
        Distr_Epsilon == "normal",
        Distr_Zeta == "normal",
        (Equation == "eta4" & Parameter %in% eta4_params) |
          (Equation == "eta5" & Parameter %in% eta5_params)
      ) %>%
      mutate(
        Distribution = DIST_LABS[dist_name],
        Parameter = case_when(
          Equation == "eta4" & Parameter == "eta1:eta2" ~ "beta[44]",
          Equation == "eta4" & Parameter == "eta1:eta3" ~ "beta[45]",
          Equation == "eta4" & Parameter == "eta1:eta1" ~ "beta[46]",
          Equation == "eta4" & Parameter == "eta2:eta2" ~ "beta[47]",
          Equation == "eta5" & Parameter == "eta1:eta4" ~ "beta[55]",
          Equation == "eta5" & Parameter == "eta2:eta4" ~ "beta[56]",
          Equation == "eta5" & Parameter == "eta1:eta1" ~ "beta[57]",
          Equation == "eta5" & Parameter == "eta3:eta3" ~ "beta[58]"
        ),
        Method = factor(Method, levels = METHOD_ORDER_3),
        Condition = make_condition(SampleSize, Reliability),
        y = Power,
        yerr = Power_MCSE
      ) %>%
      select(Distribution, Parameter, Condition, Method, y, yerr)
  )
}

# 4.3 ADDITIONAL DISTRIBUTIONAL CONDITIONS: Non-normal errors (epsilon and zeta)

#' Prepare additional distributional conditions data (nonnormal errors)
#' @param results_data Full results data frame
#' @param study Which study: 1 or 2
#' @param include_low_rel Include Rel=0.4? Default FALSE
#' @param distribution For Simulation 2: distribution name ("normal", "nonnormal", or "uniform")
#' @return Formatted data frame for additional distributional conditions analysis
prep_sensitivity <- function(results_data, study = 1, include_low_rel = FALSE, 
                             distribution = NULL) {
  
  if (study == 1) {
    PARAM_ORDER <- c("eta1:eta2", "eta1:eta1", "eta2:eta2")
    PARAM_LABS <- c("eta1:eta2" = "beta[33]",
                    "eta1:eta1" = "beta[34]",
                    "eta2:eta2" = "beta[35]")
    method_order <- METHOD_ORDER_4
    
    df <- results_data %>%
      filter(
        Parameter %in% PARAM_ORDER,
        !(Distr_Epsilon == "normal" & Distr_Zeta == "normal")
      )
    
    if (!include_low_rel) {
      df <- df %>% filter(Reliability != 0.4)
    }
    
    df %>%
      mutate(
        Method = factor(Method, levels = method_order),
        Distribution = factor(Distr_Exo, levels = names(DIST_LABS), labels = DIST_LABS),
        Parameter = factor(Parameter, levels = PARAM_ORDER, labels = PARAM_LABS),
        Reliability = factor(Reliability),
        SampleSize = factor(SampleSize),
        Condition = make_condition(SampleSize, Reliability),
        Model = factor(Model, levels = c("Linear", "Full")),
        Error_Condition = case_when(
          Distr_Epsilon == "exp.rate1" & Distr_Zeta == "normal"    ~ "Epsilon",
          Distr_Epsilon == "normal"    & Distr_Zeta == "exp.rate1" ~ "Zeta",
          Distr_Epsilon == "exp.rate1" & Distr_Zeta == "exp.rate1" ~ "Epsilon & Zeta"
        ),
        Error_Condition = factor(Error_Condition,
                                 levels = c("Epsilon", "Zeta", "Epsilon & Zeta"))
      )

  } else if (study == 2) {
    
    if (is.null(distribution)) {
      stop("For Study 2, 'distribution' must be specified ('normal', 'nonnormal', or 'uniform')")
    }
    
    eta4_params <- c("eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2")
    eta5_params <- c("eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
    method_order <- METHOD_ORDER_3
    
    df <- results_data %>%
      filter(
        Distr_Exo == distribution,
        !(Distr_Epsilon == "normal" & Distr_Zeta == "normal"),
        (Equation == "eta4" & Parameter %in% eta4_params) |
          (Equation == "eta5" & Parameter %in% eta5_params)
      )
    
    if (!include_low_rel) {
      df <- df %>% filter(Reliability != 0.4)
    }
    
    df %>%
      mutate(
        Method = factor(Method, levels = method_order),
        Distribution = factor(Distr_Exo, levels = names(DIST_LABS), labels = DIST_LABS),
        Parameter = case_when(
          Equation == "eta4" & Parameter == "eta1:eta2" ~ "beta[44]",
          Equation == "eta4" & Parameter == "eta1:eta3" ~ "beta[45]",
          Equation == "eta4" & Parameter == "eta1:eta1" ~ "beta[46]",
          Equation == "eta4" & Parameter == "eta2:eta2" ~ "beta[47]",
          Equation == "eta5" & Parameter == "eta1:eta4" ~ "beta[55]",
          Equation == "eta5" & Parameter == "eta2:eta4" ~ "beta[56]",
          Equation == "eta5" & Parameter == "eta1:eta1" ~ "beta[57]",
          Equation == "eta5" & Parameter == "eta3:eta3" ~ "beta[58]",
          TRUE ~ Parameter
        ),
        Parameter = factor(Parameter, 
                           levels = c("beta[44]", "beta[45]", "beta[46]", "beta[47]",
                                      "beta[55]", "beta[56]", "beta[57]", "beta[58]")),
        Reliability = factor(Reliability),
        SampleSize = factor(SampleSize),
        Condition = make_condition(SampleSize, Reliability),
        Model = factor(Model, levels = c("Linear", "Full")),
        Error_Condition = case_when(
          Distr_Epsilon == "exp.rate1" & Distr_Zeta == "normal"    ~ "Epsilon",
          Distr_Epsilon == "normal"    & Distr_Zeta == "exp.rate1" ~ "Zeta",
          Distr_Epsilon == "exp.rate1" & Distr_Zeta == "exp.rate1" ~ "Epsilon & Zeta"
        ),
        Error_Condition = factor(Error_Condition,
                                 levels = c("Epsilon", "Zeta", "Epsilon & Zeta"))
      )
  } else {
    stop("study must be 1 or 2")
  }
}

#' Compute data for dumbbell plot (baseline and final values with MCSE)
#' @param data_baseline Prepared baseline data
#' @param data_sensitivity Prepared sensitivity data
#' @param metric Column name for the metric
#' @param mcse_col Optional: MCSE column name. Auto-detected if NULL.
#' @param model Which model to use: "Full" or "Linear". Default is "Full".
#' @return Data frame with baseline and final values
compute_dumbbell_data <- function(data_baseline, data_sensitivity, metric, mcse_col = NULL,
                                  model = "Full") {
  
  if (is.null(mcse_col)) {
    mcse_col <- switch(metric,
                       "Bias_Mean" = "Bias_Mean_MCSE",
                       "RelativeBias_Mean" = "RelativeBias_MCSE",
                       "PercentRelativeBias_Mean" = "RelativeBias_MCSE",
                       "RMSE_Mean" = "RMSE_Mean_MCSE",
                       "Relative_RMSE" = "Relative_RMSE_MCSE",
                       "CoverageRate" = "CoverageRate_MCSE",
                       "TypeI_Error" = "TypeI_Error_MCSE",
                       "Power" = "Power_MCSE",
                       "SE_SD_Ratio" = NA_character_,
                       paste0(metric, "_MCSE"))
  }
  
  has_mcse <- !is.na(mcse_col) && 
    mcse_col %in% names(data_baseline) && 
    mcse_col %in% names(data_sensitivity)
  
  # baseline values
  baseline <- data_baseline %>%
    filter(Model == model) %>%
    select(Distribution, Parameter, SampleSize, Reliability, Method,
           y_baseline = !!sym(metric),
           n_baseline = N_Final)
  
  if (has_mcse) {
    baseline <- baseline %>%
      left_join(
        data_baseline %>%
          filter(Model == model) %>%
          select(Distribution, Parameter, SampleSize, Reliability, Method,
                 mcse_baseline = !!sym(mcse_col)),
        by = c("Distribution", "Parameter", "SampleSize", "Reliability", "Method")
      )
  } else {
    baseline$mcse_baseline <- NA_real_
  }
  
  data_sens_prep <- data_sensitivity %>%
    filter(Model == model) %>%
    mutate(y_final = !!sym(metric), n_final = N_Final)
  
  if (has_mcse) {
    data_sens_prep <- data_sens_prep %>% mutate(mcse_final = !!sym(mcse_col))
  } else {
    data_sens_prep <- data_sens_prep %>% mutate(mcse_final = NA_real_)
  }
  
  # join baseline and final
  data_sens_prep %>%
    left_join(baseline, by = c("Distribution", "Parameter", "SampleSize", 
                               "Reliability", "Method")) %>%
    select(Distribution, Parameter, SampleSize, Reliability, Method, Error_Condition,
           y_baseline, y_final, mcse_baseline, mcse_final, n_baseline, n_final)
}

#' Prepare all dumbbell data frames for additional distributional conditions analysis
#' @param results_data Raw simulation results data frame
#' @param study Which study: 1 or 2
#' @param include_low_rel Include Rel=0.4? Default FALSE
#' @param distribution For Simulation 2: distribution name ("normal", "nonnormal", or "uniform")
#' @param model Which model for main metrics: "Full" or "Linear". Default is "Full".
#'              (Type I always uses Linear, Power always uses Full)
#' @return List of dumbbell data frames for each metric
prepare_sensitivity_dumbbell_data <- function(results_data, study = 1, include_low_rel = FALSE,
                                              distribution = NULL, model = "Full") {
  
  if (study == 1) {
    data_baseline <- prep_study1(results_data)
    
    if (!include_low_rel) {
      data_baseline <- data_baseline %>% filter(Reliability != 0.4)
    }
    
  } else if (study == 2) {
    
    if (is.null(distribution)) {
      stop("For Study 2, 'distribution' must be specified ('normal', 'nonnormal', or 'uniform')")
    }
    
    data_baseline <- prep_study2_subset(results_data, distribution, "Full") %>%
      bind_rows(prep_study2_subset(results_data, distribution, "Linear")) %>%
      mutate(Model = factor(Model, levels = c("Linear", "Full")))
    
    # rename Parameter_Label to Parameter for consistency with dumbbell function
    data_baseline <- data_baseline %>%
      mutate(Parameter = Parameter_Label) %>%
      select(-Parameter_Label)
    
    if (!include_low_rel) {
      data_baseline <- data_baseline %>% filter(Reliability != 0.4)
    }
    
  } else {
    stop("study must be 1 or 2")
  }
  
  data_sensitivity <- prep_sensitivity(results_data, study = study, 
                                       include_low_rel = include_low_rel,
                                       distribution = distribution)
  
  list(
    # user-selected model metrics
    bias = compute_dumbbell_data(data_baseline, data_sensitivity, "Bias_Mean", model = model),
    bias_relative = compute_dumbbell_data(data_baseline, data_sensitivity, "RelativeBias_Mean", model = model),
    rmse = compute_dumbbell_data(data_baseline, data_sensitivity, "RMSE_Mean", model = model),
    rmse_relative = compute_dumbbell_data(data_baseline, data_sensitivity, "Relative_RMSE", model = model),
    coverage = compute_dumbbell_data(data_baseline, data_sensitivity, "CoverageRate", model = model),
    sesd = compute_dumbbell_data(data_baseline, data_sensitivity, "SE_SD_Ratio", model = model),
    
    # fixed model metrics
    power = compute_dumbbell_data(data_baseline, data_sensitivity, "Power", model = "Full"),
    type1 = compute_dumbbell_data(data_baseline, data_sensitivity, "TypeI_Error", model = "Linear")
  )
}

############################### 5. Plotting Functions ##########################

# 5.1 MAIN PLOTS (Simulations 1 and 2)

#' Plot bias (absolute or relative)
#' @param data Prepared data frame with bias values
#' @param shapes Named vector of point shapes. Default is SHAPES_4
#' @param ltys Named vector of line types. Default is LTYS_4
#' @param facet_formula Custom faceting formula. Default is NULL
#' @param y_breaks Custom y-axis breaks. Default is NULL
#' @param y_limits Custom y-axis limits. Default is NULL
#' @param bias_type "absolute" or "relative". Default is "absolute"
#' @return ggplot2 object
plot_bias <- function(data, shapes = SHAPES_4, ltys = LTYS_4,
                      facet_formula = NULL, y_breaks = NULL, y_limits = NULL,
                      bias_type = c("absolute", "relative")) {
  
  bias_type <- match.arg(bias_type)
  
  p <- ggplot(data,
              aes(x = Reliability, y = y, shape = Method, linetype = Method, group = Method)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_pointrange(aes(ymin = y - yerr, ymax = y + yerr),
                    position = position_dodge(width = .6), color = "black") +
    geom_line(position = position_dodge(width = .6), color = "black") +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    labs(x = "Reliability", y = "Bias") +
    theme_apa_bw()
  
  if (bias_type == "relative") {
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.10, ymax = 0.10,
                      fill = "grey93", alpha = .9)
    p$layers <- c(p$layers[length(p$layers)], p$layers[-length(p$layers)])
  }
  
  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula, 
                        labeller = labeller(Parameter = label_parsed,
                                            Parameter_Label = label_parsed))
  } else {
    p <- p + facet_grid(Distribution + Parameter ~ Model + SampleSize,
                        labeller = labeller(SampleSize = label_value,
                                            Parameter = label_parsed))
  }
  
  if (!is.null(y_breaks) || !is.null(y_limits)) {
    p <- p + scale_y_continuous(breaks = y_breaks, limits = y_limits)
  }
  
  p
}

#' Plot RMSE
#' @param data Prepared data frame with RMSE values
#' @param shapes Named vector of point shapes. Default is SHAPES_4
#' @param ltys Named vector of line types. Default is LTYS_4
#' @param facet_formula Custom faceting formula. Default is NULL
#' @param y_breaks Custom y-axis breaks. Default is NULL
#' @param y_limits Custom y-axis limits. Default is NULL
#' @return ggplot2 object
plot_rmse <- function(data, shapes = SHAPES_4, ltys = LTYS_4,
                      facet_formula = NULL, y_breaks = NULL, y_limits = NULL) {
  
  p <- ggplot(data,
              aes(x = Reliability, y = y, shape = Method, linetype = Method, group = Method)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_pointrange(aes(ymin = y - yerr, ymax = y + yerr),
                    position = position_dodge(width = .6), color = "black") +
    geom_line(position = position_dodge(width = .6), color = "black") +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    labs(x = "Reliability", y = "RMSE") +
    theme_apa_bw()
  
  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula, 
                        labeller = labeller(Parameter = label_parsed,
                                            Parameter_Label = label_parsed))
  } else {
    p <- p + facet_grid(Distribution + Parameter ~ Model + SampleSize,
                        labeller = labeller(SampleSize = label_value,
                                            Parameter = label_parsed))
  }
  
  if (!is.null(y_breaks) || !is.null(y_limits)) {
    p <- p + scale_y_continuous(breaks = y_breaks, limits = y_limits)
  }
  
  p
}

#' Plot SE/SD ratio
#' @param data Prepared data frame with SE/SD ratio values
#' @param shapes Named vector of point shapes. Default is SHAPES_4
#' @param ltys Named vector of line types. Default is LTYS_4
#' @param facet_formula Custom faceting formula. Default is NULL
#' @return ggplot2 object
plot_sesd <- function(data, shapes = SHAPES_4, ltys = LTYS_4,
                      facet_formula = NULL) {
  
  p <- ggplot(data,
              aes(x = Reliability, y = ratio, 
                  group = Method, shape = Method, linetype = Method)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.90, ymax = 1.10,
             fill = "grey93", alpha = .9) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_line(color = "black", position = position_dodge(width = 0.6)) +
    geom_point(color = "black", size = 2, position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    coord_cartesian(ylim = c(0.80, 1.40)) +
    labs(x = "Reliability", y = "SE/SD") +
    theme_apa_bw()
  
  p <- p + 
    geom_text(data = filter(data, ratio > 1.40) %>%
                group_by(Distribution, Parameter, Model, SampleSize, Reliability) %>%
                arrange(desc(ratio)) %>%
                mutate(label_y = 1.39 - (row_number() - 1) * 0.03),
              aes(label = sprintf("%s: %.2f", Method, ratio), y = label_y),
              size = 1.5, hjust = 0, vjust = 1) +
    geom_text(data = filter(data, ratio < 0.80) %>%
                group_by(Distribution, Parameter, Model, SampleSize, Reliability) %>%
                arrange(ratio) %>%
                mutate(label_y = 0.81 + (row_number() - 1) * 0.03),
              aes(label = sprintf("%s: %.2f", Method, ratio), y = label_y),
              size = 1.5, hjust = 0, vjust = 0)
  
  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula, 
                        labeller = labeller(Parameter = label_parsed,
                                            Parameter_Label = label_parsed))
  } else {
    p <- p + facet_grid(Distribution + Parameter ~ Model + SampleSize,
                        labeller = labeller(SampleSize = label_value,
                                            Parameter = label_parsed))
  }
  
  p
}

#' Plot coverage rate
#' @param data Prepared data frame with coverage rates
#' @param shapes Named vector of point shapes. Default is SHAPES_4
#' @param ltys Named vector of line types. Default is LTYS_4
#' @param facet_formula Custom faceting formula. Default is NULL
#' @return ggplot2 object
plot_coverage <- function(data, shapes = SHAPES_4, ltys = LTYS_4,
                          facet_formula = NULL) {
  
  p <- ggplot(data,
              aes(x = Reliability, y = y, shape = Method, linetype = Method, group = Method)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 93, ymax = 97,
             fill = "grey93", alpha = .9) +
    geom_hline(yintercept = 95, linetype = "dotted") +
    geom_pointrange(aes(ymin = pmax(y - yerr, 80), 
                        ymax = pmin(y + yerr, 100)),
                    position = position_dodge(width = .6), color = "black") +
    geom_line(position = position_dodge(width = .6), color = "black") +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    coord_cartesian(ylim = c(80, 100)) +
    labs(x = "Reliability", y = "Coverage (%)") +
    theme_apa_bw()
  
  p <- p +
    geom_text(data = filter(data, y < 80) %>%
                group_by(Distribution, Parameter, Model, SampleSize, Reliability) %>%
                arrange(Method) %>%
                mutate(label_y = 81 + (row_number() - 1) * 2.5),
              aes(label = sprintf("%s: %.0f%%", Method, y), y = label_y),
              size = 1.5, hjust = 0, vjust = 0)
  
  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula, 
                        labeller = labeller(Parameter = label_parsed,
                                            Parameter_Label = label_parsed))
  } else {
    p <- p + facet_grid(Distribution + Parameter ~ Model + SampleSize,
                        labeller = labeller(SampleSize = label_value,
                                            Parameter = label_parsed))
  }
  
  p
}

#' Plot Type I error rate
#' @param data Prepared data frame with Type I error rates (includes yerr for MCSE)
#' @param greys Named vector of grey colors. Default is GREYS_4
#' @param facet_formula Custom faceting formula. Default is NULL
#' @param show_errorbars Logical. Show MCSE error bars. Default is TRUE
#' @return ggplot2 object
plot_type1 <- function(data, greys = GREYS_4, facet_formula = NULL,
                       show_errorbars = TRUE) {
  
  ymax <- max(data$y + ifelse(is.null(data$yerr), 0, data$yerr), na.rm = TRUE) * 1.15
  
  p <- ggplot(data, aes(x = Method, y = y, fill = Method)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3, ymax = 7,
             fill = "grey92", alpha = .9) +
    geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.3) +
    geom_col(position = position_dodge(width = 0.9), width = 0.75, color = "black")
  
  if (show_errorbars && "yerr" %in% names(data) && !all(is.na(data$yerr))) {
    p <- p + geom_errorbar(aes(ymin = y - yerr, ymax = y + yerr),
                           width = 0.2, position = position_dodge(width = 0.9))
  }
  
  p <- p +
    geom_text(aes(label = sprintf("%.1f", y)), vjust = -0.4, size = 3) +
    scale_fill_manual(values = greys) +
    coord_cartesian(ylim = c(0, ymax)) +
    labs(x = NULL, y = "Type I error (%)") +
    theme_apa_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula, 
                        labeller = labeller(Parameter = label_parsed,
                                            Parameter_Label = label_parsed))
  } else {
    p <- p + facet_grid(Distribution + Parameter ~ Condition,
                        labeller = labeller(Parameter = label_parsed))
  }
  
  p
}

#' Plot statistical power
#' @param data Prepared data frame with power values (includes yerr for MCSE)
#' @param greys Named vector of grey colors. Default is GREYS_4
#' @param facet_formula Custom faceting formula. Default is NULL
#' @param show_errorbars Logical. Show MCSE error bars. Default is TRUE
#' @return ggplot2 object
plot_power <- function(data, greys = GREYS_4, facet_formula = NULL,
                       show_errorbars = TRUE) {
  
  p <- ggplot(data, aes(x = Method, y = y, fill = Method)) +
    geom_hline(yintercept = 80, linetype = "dotted", linewidth = 0.3) +
    geom_col(position = position_dodge(width = 0.9), width = 0.75, color = "black")
  
  if (show_errorbars && "yerr" %in% names(data) && !all(is.na(data$yerr))) {
    p <- p + geom_errorbar(aes(ymin = y - yerr, ymax = y + yerr),
                           width = 0.2, position = position_dodge(width = 0.9))
  }
  
  p <- p +
    geom_text(aes(label = sprintf("%.0f", y)), vjust = -0.4, size = 3) +
    scale_fill_manual(values = greys) +
    coord_cartesian(ylim = c(0, 110)) +
    labs(x = NULL, y = "Power (%)") +
    theme_apa_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula, 
                        labeller = labeller(Parameter = label_parsed,
                                            Parameter_Label = label_parsed))
  } else {
    p <- p + facet_grid(Distribution + Parameter ~ Condition,
                        labeller = labeller(Parameter = label_parsed))
  }
  
  p
}

# 5.2 ADDITIONAL DISTRIBUTIONAL CONDITIONS PLOTS

#' Horizontal dumbbell plot for additional distributional conditions analysis
#' @param data Data from compute_dumbbell_data()
#' @param metric_name Name for x-axis label. Also used to determine default acceptable_range and reference_line.
#' @param show_errorbars Show 95% CI around the additional distributional conditions estimate (black dot). Default TRUE
#' @param acceptable_range Optional vector of length 2 for shaded region. If NULL, uses defaults based on metric_name.
#' @param reference_line Optional reference line value. If NULL, uses defaults based on metric_name.
#' @param errorbar_height Height of error bar caps. Default 0.5
#' @return ggplot2 object
plot_dumbbell <- function(data,
                          metric_name = "Relative Bias",
                          show_errorbars = TRUE,
                          acceptable_range = NULL,
                          reference_line = NULL,
                          errorbar_height = 0.5,
                          x_breaks = NULL,
                          x_limits = NULL) {

  # default acceptable ranges and reference lines based on metric name
  # (matching the main plots for consistency)
  metric_defaults <- list(
    "Relative Bias" = list(range = c(-0.10, 0.10), ref = 0),
    "Bias" = list(range = NULL, ref = 0),
    "RMSE" = list(range = NULL, ref = 0),
    "Relative RMSE" = list(range = NULL, ref = 0),
    "SE/SD Ratio" = list(range = c(0.90, 1.10), ref = 1),
    "SE/SD" = list(range = c(0.90, 1.10), ref = 1),
    "Coverage (%)" = list(range = c(93, 97), ref = 95),
    "Coverage" = list(range = c(93, 97), ref = 95),
    "Type I Error (%)" = list(range = c(3, 7), ref = 5),
    "Type I Error" = list(range = c(3, 7), ref = 5),
    "Power (%)" = list(range = NULL, ref = 80),
    "Power" = list(range = NULL, ref = 80)
  )

  # use defaults if not specified
  if (is.null(acceptable_range) && metric_name %in% names(metric_defaults)) {
    acceptable_range <- metric_defaults[[metric_name]]$range
  }
  if (is.null(reference_line) && metric_name %in% names(metric_defaults)) {
    reference_line <- metric_defaults[[metric_name]]$ref
  }

  p <- ggplot(data)

  # add acceptable range shading first (so it's behind everything)
  if (!is.null(acceptable_range)) {
    p <- p + annotate("rect",
                      ymin = -Inf, ymax = Inf,
                      xmin = acceptable_range[1], xmax = acceptable_range[2],
                      fill = "grey93", alpha = 0.9)
  }

  # add reference line
  if (!is.null(reference_line)) {
    p <- p + geom_vline(xintercept = reference_line, linetype = "dotted", color = "grey50")
  }

  # segment connecting baseline to additional distributional conditions
  p <- p +
    geom_segment(
      aes(y = Method, yend = Method, x = y_baseline, xend = y_final),
      color = "grey40", linewidth = 0.8
    )

  # error bars for the nonnormal estimate (black dot) using its own MCSE
  # note: MCSE for the difference is computed separately for tabular reporting
  if (show_errorbars) {
    p <- p +
      geom_errorbarh(
        aes(y = Method,
            xmin = y_final - mcse_final,
            xmax = y_final + mcse_final),
        height = errorbar_height, linewidth = 0.5, color = "black"
      )
  }

  p <- p +
    geom_point(aes(y = Method, x = y_baseline),
               shape = 21, size = 2.5, fill = "white", color = "grey50", stroke = 1) +
    geom_point(aes(y = Method, x = y_final),
               shape = 21, size = 2.5, fill = "black", color = "black") +
    facet_grid(Distribution + Parameter ~ Error_Condition + SampleSize + Reliability,
               labeller = labeller(Parameter = label_parsed)) +
    labs(y = NULL, x = metric_name) +
    theme_apa_bw()

  if (!is.null(x_breaks) || !is.null(x_limits)) {
    p <- p + scale_x_continuous(breaks = x_breaks, limits = x_limits)
  }

  p
}

############ 6. Additional Distributional Conditions Difference Functions ########

#' Compute difference between additional distributional conditions and baseline with 95% CI
#' @param dumbbell_data Data frame from compute_dumbbell_data()
#' @param metric_name Character string labeling the metric (e.g., "Bias", "RMSE")
#' @return Data frame with difference, MCSE of difference, and 95% CI
compute_sensitivity_difference <- function(dumbbell_data, metric_name) {
  dumbbell_data %>%
    mutate(
      Metric     = metric_name,
      Difference = y_final - y_baseline,
      MCSE_Diff  = sqrt(mcse_baseline^2 + mcse_final^2),
      CI_Lower   = Difference - 1.96 * MCSE_Diff,
      CI_Upper   = Difference + 1.96 * MCSE_Diff
    ) %>%
    select(Distribution, Parameter, SampleSize, Reliability, Method,
           Error_Condition, Metric, Difference, MCSE_Diff, CI_Lower, CI_Upper,
           n_baseline, n_final)
}

#' Compute additional distributional conditions differences for all metrics from prepare_sensitivity_dumbbell_data()
#' @param dumbbell_list Named list output of prepare_sensitivity_dumbbell_data()
#' @return Combined data frame with a Metric column
compute_all_sensitivity_differences <- function(dumbbell_list) {
  metric_labels <- c(
    bias          = "Bias",
    bias_relative = "Relative Bias",
    rmse          = "RMSE",
    rmse_relative = "Relative RMSE",
    coverage      = "Coverage Rate",
    sesd          = "SE/SD Ratio",
    power         = "Power",
    type1         = "Type I Error"
  )

  do.call(rbind, lapply(names(dumbbell_list), function(nm) {
    label <- if (nm %in% names(metric_labels)) metric_labels[[nm]] else nm
    compute_sensitivity_difference(dumbbell_list[[nm]], label)
  }))
}

#' Bar plot of additional distributional conditions differences (nonnormal - baseline) with 95% CI
#' @param data Output of compute_sensitivity_difference() filtered to one metric + one Error_Condition
#' @param metric_name Label for the y-axis (e.g., "Relative Bias")
#' @param greys Named vector of grey fill colours. Default is GREYS_4
#' @param show_errorbars Show 95% CI error bars. Default TRUE
#' @return ggplot2 object
plot_sensitivity_difference <- function(data, metric_name,
                                        greys = GREYS_4,
                                        show_errorbars = TRUE) {

  p <- ggplot(data, aes(x = Method, y = Difference, fill = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_col(position = position_dodge(width = 0.9), width = 0.75, color = "black")

  if (show_errorbars && !all(is.na(data$MCSE_Diff))) {
    p <- p + geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                           width = 0.2, position = position_dodge(width = 0.9))
  }

  p + scale_fill_manual(values = greys) +
    labs(x = NULL, y = paste0("\u0394 ", metric_name)) +
    theme_apa_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    facet_grid(Distribution + Parameter ~ SampleSize + Reliability,
               labeller = labeller(Parameter = label_parsed))
}