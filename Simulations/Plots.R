############################ 1. General Information ############################

# This file contains all the code for data wrangling and plots based on the 
# simulation results (with various performance metrics):
# Consistent "APA-style" formatting across all plots
# Data preparation functions for different sim designs
# Plotting functions for bias, RMSE, SE/SD ratio, coverage, Type I error, and power

# Packages needed for this script:
# library(ggplot2); library(dplyr)

############################### 2. Configuration ################################

# METHOD SPECIFICATIONS 
# method orderings and labels (study-specific)
METHOD_ORDER_3 <- c("LSAM","LMS","QML","UPI")  
METHOD_ORDER_5 <- c("LSAM","QML","UPI")        
# distribution labels for all studies
DIST_LABS <- c(normal="Normal", nonnormal="Right-skewed", uniform="Uniform")

# VISUAL SPECIFICATIONS 
# shape specifications for point markers
SHAPES_4 <- c(LSAM=16, LMS=17, QML=15, UPI=18)  
SHAPES_3 <- c(LSAM=16, QML=15, UPI=18)          
# line type specifications
LTYS_4 <- c(LSAM="solid", LMS="dashed", QML="dotdash", UPI="twodash")
LTYS_3 <- c(LSAM="solid", QML="dotdash", UPI="twodash")
# greyscale colors for bar plots
GREYS_4 <- c(LSAM="grey20", LMS="grey45", QML="grey65", UPI="grey85")
GREYS_3 <- c(LSAM="grey20", QML="grey50", UPI="grey80")

############################### 3. Shared Components ############################

# APA THEME FUNCTION 
#' @param base_size Base font size for all text elements. Default is 11.
#' @return A ggplot2 theme object with APA-style formatting
theme_apa_bw <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

# HELPER FUNCTION
#' condition labels from sample size and reliability
#' @param N Sample size values
#' @param Rel Reliability values
#' @return Factor with properly ordered condition labels
make_condition <- function(N, Rel) {
  factor(paste0("N=", N, ", Rel=", Rel),
         levels = c("N=400, Rel=0.4","N=400, Rel=0.6","N=400, Rel=0.8",
                    "N=1000, Rel=0.4","N=1000, Rel=0.6","N=1000, Rel=0.8"))
}

############################### 4. Data Preparation ##############################

#  SIMULATION 1 
#' @param df Data frame with simulation results
#' @param model Either "Full" or "Linear" model specification
#' @return Formatted data frame ready for plotting
prep_study1_subset <- function(df, model = c("Full","Linear")) {
  model <- match.arg(model)
  
  PARAM_ORDER <- c("eta1:eta2","eta1:eta1","eta2:eta2")
  PARAM_LABS <- c("eta1:eta2"="beta[33]",
                  "eta1:eta1"="beta[34]",
                  "eta2:eta2"="beta[35]")
  
  df %>%
    filter(Parameter %in% PARAM_ORDER, Model == model) %>%
    mutate(
      Method = factor(Method, levels = METHOD_ORDER_3),
      Distribution = factor(Distribution, levels = names(DIST_LABS), labels = DIST_LABS),
      Parameter = factor(Parameter, levels = PARAM_ORDER, labels = PARAM_LABS),
      Reliability = factor(Reliability),
      SampleSize = factor(SampleSize),
      Condition = make_condition(SampleSize, Reliability)
    )
}

# SIMULATION 2 
#' @param df Data frame with simulation results
#' @param distribution Distribution name ("normal", "nonnormal", or "uniform")
#' @param model Either "Full" or "Linear" model specification
#' @return Formatted data frame ready for plotting
prep_study2_subset <- function(df, distribution, model = c("Full","Linear")) {
  model <- match.arg(model)
  
  eta4_params <- c("eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2")
  eta5_params <- c("eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
  
  df %>%
    filter(Distribution == distribution,
           Model == model,
           (Equation == "eta4" & Parameter %in% eta4_params) |
             (Equation == "eta5" & Parameter %in% eta5_params)) %>%
    mutate(
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
      Method = factor(Method, levels = METHOD_ORDER_5),
      Reliability = factor(Reliability),
      SampleSize = factor(SampleSize),
      Condition = make_condition(SampleSize, Reliability)
    ) %>%
    mutate(Parameter_Label = factor(Parameter_Label, 
                                    levels = c("beta[44]", "beta[45]", "beta[46]", "beta[47]",
                                               "beta[55]", "beta[56]", "beta[57]", "beta[58]")))
}

# MAIN DATA 
# all data frames for Study 1 plotting
#' @param results_data Raw simulation results data frame
#' @return List containing formatted data frames for each metric type
prepare_study1_data <- function(results_data) {
  
  # bias data (absolute)
  dat_bias <- bind_rows(
    prep_study1_subset(results_data, "Linear") %>% mutate(Model = "Linear"),
    prep_study1_subset(results_data, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
    transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
              y = Bias_Mean, yerr = Bias_Mean_MCSE)
  
  # relative bias data
  dat_bias_relative <- bind_rows(
    prep_study1_subset(results_data, "Linear") %>% mutate(Model = "Linear"),
    prep_study1_subset(results_data, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
    transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
              y = RelativeBias_Mean, yerr = RelativeBias_MCSE)
  
  # RMSE data
  dat_rmse <- bind_rows(
    prep_study1_subset(results_data, "Linear") %>% mutate(Model = "Linear"),
    prep_study1_subset(results_data, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
    transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
              y = RMSE_Mean, yerr = RMSE_Mean_MCSE)
  
  # relative RMSE data
  dat_rmse_relative <- bind_rows(
    prep_study1_subset(results_data, "Linear") %>% mutate(Model = "Linear"),
    prep_study1_subset(results_data, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
    transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
              y = Relative_RMSE, yerr = Relative_RMSE_MCSE)
  
  # se/sd ratio data
  dat_sesd <- bind_rows(
    prep_study1_subset(results_data, "Linear") %>% mutate(Model = "Linear"),
    prep_study1_subset(results_data, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
    transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
              ratio = SE_SD_Ratio)
  
  # coverage rate data
  dat_cov <- bind_rows(
    prep_study1_subset(results_data, "Linear") %>% mutate(Model = "Linear"),
    prep_study1_subset(results_data, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
    transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
              y = CoverageRate, yerr = 100 * CoverageRate_MCSE)
  
  # Type I error data (linear model only)
  dat_t1 <- prep_study1_subset(results_data, "Linear") %>%
    transmute(Distribution, Parameter, Condition, Method, y = TypeI_Error)
  
  # statistical power data (full model only)
  dat_pow <- prep_study1_subset(results_data, "Full") %>%
    transmute(Distribution, Parameter, Condition, Method, y = Power)
  
  list(
    bias = dat_bias,
    bias_relative = dat_bias_relative,
    rmse = dat_rmse,
    rmse_relative = dat_rmse_relative,
    sesd = dat_sesd,
    coverage = dat_cov,
    type1 = dat_t1,
    power = dat_pow
  )
}

# all data frames for Study 2 plotting (per distribution)
#' @param results_data Raw simulation results data frame
#' @param dist_name Distribution name ("normal", "nonnormal", or "uniform")
#' @return List containing formatted data frames for each metric type
prepare_study2_data <- function(results_data, dist_name) {
  
  # bias data (absolute)
  dat_bias <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")),
           Distribution = DIST_LABS[dist_name]) %>%
    transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability, 
              Method, Model, y = Bias_Mean, yerr = Bias_Mean_MCSE)
  
  # relative bias data
  dat_bias_relative <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")),
           Distribution = DIST_LABS[dist_name]) %>%
    transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability, 
              Method, Model, y = RelativeBias_Mean, yerr = RelativeBias_MCSE)
  
  # RMSE data
  dat_rmse <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")),
           Distribution = DIST_LABS[dist_name]) %>%
    transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
              Method, Model, y = RMSE_Mean, yerr = RMSE_Mean_MCSE)
  
  # relative RMSE data
  dat_rmse_relative <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")),
           Distribution = DIST_LABS[dist_name]) %>%
    transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
              Method, Model, y = Relative_RMSE, yerr = Relative_RMSE_MCSE)
  
  # se/sd ratio data
  dat_sesd <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")),
           Distribution = DIST_LABS[dist_name]) %>%
    transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
              Method, Model, ratio = SE_SD_Ratio)
  
  # coverage rate data
  dat_cov <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")),
           Distribution = DIST_LABS[dist_name]) %>%
    transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
              Method, Model, y = CoverageRate, yerr = 100 * CoverageRate_MCSE)
  
  # Type I error and power data (special handling for eta4/eta5 parameters)
  eta4_params <- c("eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2")
  eta5_params <- c("eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
  
  # Type I error (linear model only)
  dat_t1 <- results_data %>%
    filter(Model == "Linear", Distribution == dist_name,
           (Equation == "eta4" & Parameter %in% eta4_params) |
             (Equation == "eta5" & Parameter %in% eta5_params)) %>%
    mutate(
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
      Method = factor(Method, levels = METHOD_ORDER_5),
      Condition = make_condition(SampleSize, Reliability),
      Distribution = DIST_LABS[dist_name],
      y = TypeI_Error
    ) %>%
    select(Distribution, Parameter, Condition, Method, y)
  
  # statistical power (full model only)
  dat_pow <- results_data %>%
    filter(Model == "Full", Distribution == dist_name,
           (Equation == "eta4" & Parameter %in% eta4_params) |
             (Equation == "eta5" & Parameter %in% eta5_params)) %>%
    mutate(
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
      Method = factor(Method, levels = METHOD_ORDER_5),
      Condition = make_condition(SampleSize, Reliability),
      Distribution = DIST_LABS[dist_name],
      y = Power
    ) %>%
    select(Distribution, Parameter, Condition, Method, y)
  
  list(
    bias = dat_bias,
    bias_relative = dat_bias_relative,
    rmse = dat_rmse,
    rmse_relative = dat_rmse_relative,
    sesd = dat_sesd,
    coverage = dat_cov,
    type1 = dat_t1,
    power = dat_pow
  )
}

############################### 5. Plotting Functions ############################

#' @param data Prepared data frame with bias values
#' @param shapes Named vector of point shapes for methods. Default is SHAPES_4
#' @param ltys Named vector of line types for methods. Default is LTYS_4
#' @param title Plot title. Default is ""
#' @param facet_formula Custom faceting formula. Default is NULL (uses default faceting)
#' @param y_breaks Custom y-axis breaks. Default is NULL
#' @param y_limits Custom y-axis limits. Default is NULL
#' @param bias_type Type of bias plot: "absolute" or "relative". Default is "absolute"
#' @return ggplot2 object
#' 
plot_bias <- function(data, shapes = SHAPES_4, ltys = LTYS_4, title = "", 
                      facet_formula = NULL, y_breaks = NULL, y_limits = NULL,
                      bias_type = c("absolute", "relative")) {
  
  bias_type <- match.arg(bias_type)
  
  p <- ggplot(data,
              aes(x = Reliability, y = y, shape = Method, linetype = Method, group = Method)) +
    geom_pointrange(aes(ymin = y - yerr, ymax = y + yerr),
                    position = position_dodge(width = .6), color = "black") +
    geom_line(position = position_dodge(width = .6), color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    labs(x = "Reliability", y = "Bias", title = title) +
    theme_apa_bw()
  
  # shaded region only for relative bias (±10% acceptable range)
  if (bias_type == "relative") {
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.10, ymax = 0.10,
                      fill = "grey93", alpha = .9)
    # shaded region to the back
    p$layers <- c(p$layers[length(p$layers)], p$layers[-length(p$layers)])
  }
  
  # faceting if provided
  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula, 
                        labeller = labeller(Parameter = label_parsed,
                                            Parameter_Label = label_parsed))
  } else {
    # faceting
    p <- p + facet_grid(Distribution + Parameter ~ Model + SampleSize,
                        labeller = labeller(SampleSize = label_value,
                                            Parameter = label_parsed))
  }
  
  # y-axis if provided
  if (!is.null(y_breaks) || !is.null(y_limits)) {
    p <- p + scale_y_continuous(breaks = y_breaks, limits = y_limits)
  }
  
  p
}

#' @param data Prepared data frame with RMSE values
#' @param shapes Named vector of point shapes for methods. Default is SHAPES_4
#' @param ltys Named vector of line types for methods. Default is LTYS_4
#' @param title Plot title. Default is ""
#' @param facet_formula Custom faceting formula. Default is NULL
#' @param y_breaks Custom y-axis breaks. Default is NULL
#' @param y_limits Custom y-axis limits. Default is NULL
#' @return ggplot2 object
#' 
plot_rmse <- function(data, shapes = SHAPES_4, ltys = LTYS_4, title = "", 
                      facet_formula = NULL, y_breaks = NULL, y_limits = NULL) {
  
  p <- ggplot(data,
              aes(x = Reliability, y = y, shape = Method, linetype = Method, group = Method)) +
    geom_pointrange(aes(ymin = y - yerr, ymax = y + yerr),
                    position = position_dodge(width = .6), color = "black") +
    geom_line(position = position_dodge(width = .6), color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    labs(x = "Reliability", y = "RMSE", title = title) +
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

#' @param data Prepared data frame with SE/SD ratio values
#' @param shapes Named vector of point shapes for methods. Default is SHAPES_4
#' @param ltys Named vector of line types for methods. Default is LTYS_4
#' @param title Plot title. Default is ""
#' @param facet_formula Custom faceting formula. Default is NULL
#' @return ggplot2 object
#' 
plot_sesd <- function(data, shapes = SHAPES_4, ltys = LTYS_4, title = "", 
                      facet_formula = NULL) {
  
  p <- ggplot(data,
              aes(x = Reliability, y = ratio, 
                  group = Method, shape = Method, linetype = Method)) +
    # shaded region for acceptable SE/SD ratios (0.90-1.10)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.90, ymax = 1.10,
             fill = "grey93", alpha = .9) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    # plot actual values, not truncated - let coord_cartesian handle the limits
    geom_line(color = "black", position = position_dodge(width = 0.6)) +
    geom_point(color = "black", size = 2, position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    coord_cartesian(ylim = c(0.80, 1.40)) +  # clips the view without changing data
    labs(x = "Reliability", y = "SE / SD", title = title) +
    theme_apa_bw()
  
  # labels for values outside the visible range
  p <- p + 
    # upper extreme values (above 1.40)
    geom_text(data = filter(data, ratio > 1.40) %>%
                group_by(Distribution, Parameter, Model, SampleSize, Reliability) %>%
                arrange(desc(ratio)) %>%
                mutate(label_y = 1.39 - (row_number() - 1) * 0.03),  
              aes(label = sprintf("%s: %.2f", Method, ratio), y = label_y),
              size = 1.5, hjust = 0, vjust = 1) +
    # lower extreme values (below 0.80)
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

#' @param data Prepared data frame with coverage rates
#' @param shapes Named vector of point shapes for methods. Default is SHAPES_4
#' @param ltys Named vector of line types for methods. Default is LTYS_4
#' @param title Plot title. Default is ""
#' @param facet_formula Custom faceting formula. Default is NULL
#' @return ggplot2 object
#' 
plot_coverage <- function(data, shapes = SHAPES_4, ltys = LTYS_4, title = "", 
                          facet_formula = NULL) {
  
  p <- ggplot(data,
              aes(x = Reliability, y = y, shape = Method, linetype = Method, group = Method)) +
    # acceptable coverage range (93-97%)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 93, ymax = 97,
             fill = "grey93", alpha = .9) +
    geom_pointrange(aes(ymin = pmax(y - yerr, 80), 
                        ymax = pmin(y + yerr, 100)),
                    position = position_dodge(width = .6), color = "black") +
    geom_line(position = position_dodge(width = .6), color = "black") +
    geom_hline(yintercept = 95, linetype = "dotted") +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    coord_cartesian(ylim = c(80, 100)) +
    labs(x = "Reliability", y = "Coverage (%)", title = title) +
    theme_apa_bw()
  
  # labels for values below 80%
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


#' @param data Prepared data frame with Type I error rates
#' @param greys Named vector of grey colors for methods. Default is GREYS_4
#' @param title Plot title. Default is ""
#' @param facet_formula Custom faceting formula. Default is NULL
#' @return ggplot2 object
#' 
plot_type1 <- function(data, greys = GREYS_4, title = "", facet_formula = NULL) {
  ymax <- max(data$y, na.rm = TRUE) * 1.15
  
  p <- ggplot(data, aes(x = Method, y = y, fill = Method)) +
    # acceptable Type I error range (3-7%)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3, ymax = 7,
             fill = "grey92", alpha = .9) +
    geom_col(position = position_dodge(width = 0.9), width = 0.75, color = "black") +
    geom_text(aes(label = sprintf("%.1f", y)), vjust = -0.4, size = 3) +
    geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.3) +
    scale_fill_manual(values = greys) +
    coord_cartesian(ylim = c(0, ymax)) +
    labs(x = NULL, y = "Type I error (%)", title = title) +
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

#' @param data Prepared data frame with power values
#' @param greys Named vector of grey colors for methods. Default is GREYS_4
#' @param title Plot title. Default is ""
#' @param facet_formula Custom faceting formula. Default is NULL
#' @return ggplot2 object
#' 
plot_power <- function(data, greys = GREYS_4, title = "", facet_formula = NULL) {
  
  p <- ggplot(data, aes(x = Method, y = y, fill = Method)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.75, color = "black") +
    geom_text(aes(label = sprintf("%.0f", y)), vjust = -0.4, size = 3) +
    geom_hline(yintercept = 80, linetype = "dotted", linewidth = 0.3) +
    scale_fill_manual(values = greys) +
    coord_cartesian(ylim = c(0, 110)) +
    labs(x = NULL, y = "Power (%)", title = title) +
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