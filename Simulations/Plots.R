############################ 1. General Information ############################

# This file contains a flexible plotting system for visualizing SEM simulation 
# study results. Supports both Study 1 and 2 with various performance metrics.

# Consistent APA-style formatting across all plots
# Data preparation functions for different study designs
# Plotting functions for bias, RMSE, SE/SD ratio, coverage, Type I error, and power

############################### 2. Configuration ################################

# METHOD SPECIFICATIONS 

# method orderings and labels (study-specific)
METHOD_ORDER_3 <- c("LSAM","LMS","QML","UPI")  # For 3-factor study (Study 1)
METHOD_ORDER_5 <- c("LSAM","QML","UPI")        # For 5-factor study (Study 2)

# distribution labels for all studies
DIST_LABS <- c(normal="Normal", nonnormal="Nonnormal", uniform="Uniform")

# ============================= VISUAL SPECIFICATIONS ===========================
# shape specifications for point markers
SHAPES_4 <- c(LSAM=16, LMS=17, QML=15, UPI=18)  # 4 methods
SHAPES_3 <- c(LSAM=16, QML=15, UPI=18)          # 3 methods

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

#  STUDY 1 

#' @param df Data frame with simulation results
#' @param model Either "Full" or "Linear" model specification
#' @return Formatted data frame ready for plotting
prep_study1_subset <- function(df, model = c("Full","Linear")) {
  model <- match.arg(model)
  
  PARAM_ORDER <- c("eta1:eta2","eta1:eta1","eta2:eta2")
  PARAM_LABS <- c("eta1:eta2"="beta[3]",
                  "eta1:eta1"="beta[4]",
                  "eta2:eta2"="beta[5]")
  
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

# STUDY 2 

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
        Equation == "eta4" & Parameter == "eta1:eta2" ~ "beta[14]",
        Equation == "eta4" & Parameter == "eta1:eta3" ~ "beta[15]",
        Equation == "eta4" & Parameter == "eta1:eta1" ~ "beta[16]",
        Equation == "eta4" & Parameter == "eta2:eta2" ~ "beta[17]",
        Equation == "eta5" & Parameter == "eta1:eta4" ~ "beta[25]",
        Equation == "eta5" & Parameter == "eta2:eta4" ~ "beta[26]",
        Equation == "eta5" & Parameter == "eta1:eta1" ~ "beta[27]",
        Equation == "eta5" & Parameter == "eta3:eta3" ~ "beta[28]",
        TRUE ~ Parameter
      ),
      Method = factor(Method, levels = METHOD_ORDER_5),
      Reliability = factor(Reliability),
      SampleSize = factor(SampleSize),
      Condition = make_condition(SampleSize, Reliability)
    ) %>%
    mutate(Parameter_Label = factor(Parameter_Label, 
                                    levels = c("beta[14]", "beta[15]", "beta[16]", "beta[17]",
                                               "beta[25]", "beta[26]", "beta[27]", "beta[28]")))
}

# MAIN DATA 

# all data frames for Study 1 plotting
#' @param results_data Raw simulation results data frame
#' @return List containing formatted data frames for each metric type
prepare_study1_data <- function(results_data) {
  
  # Bias data (absolute)
  dat_bias <- bind_rows(
    prep_study1_subset(results_data, "Linear") %>% mutate(Model = "Linear"),
    prep_study1_subset(results_data, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
    transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
              y = Bias_Mean, yerr = Bias_Mean_MCSE)
  
  # Relative bias data
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
  
  # Relative RMSE data
  dat_rmse_relative <- bind_rows(
    prep_study1_subset(results_data, "Linear") %>% mutate(Model = "Linear"),
    prep_study1_subset(results_data, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
    transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
              y = Relative_RMSE, yerr = Relative_RMSE_MCSE)
  
  # SE/SD ratio data
  dat_sesd <- bind_rows(
    prep_study1_subset(results_data, "Linear") %>% mutate(Model = "Linear"),
    prep_study1_subset(results_data, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
    transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
              ratio = SE_SD_Ratio)
  
  # Coverage rate data
  dat_cov <- bind_rows(
    prep_study1_subset(results_data, "Linear") %>% mutate(Model = "Linear"),
    prep_study1_subset(results_data, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
    transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
              y = CoverageRate, yerr = 100 * CoverageRate_MCSE)
  
  # Type I error data (Linear model only)
  dat_t1 <- prep_study1_subset(results_data, "Linear") %>%
    transmute(Distribution, Parameter, Condition, Method, y = TypeI_Error)
  
  # Statistical power data (Full model only)
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
  
  # Bias data (absolute)
  dat_bias <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")),
           Distribution = DIST_LABS[dist_name]) %>%
    transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability, 
              Method, Model, y = Bias_Mean, yerr = Bias_Mean_MCSE)
  
  # Relative bias data
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
  
  # Relative RMSE data
  dat_rmse_relative <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")),
           Distribution = DIST_LABS[dist_name]) %>%
    transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
              Method, Model, y = Relative_RMSE, yerr = Relative_RMSE_MCSE)
  
  # SE/SD ratio data
  dat_sesd <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")),
           Distribution = DIST_LABS[dist_name]) %>%
    transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
              Method, Model, ratio = SE_SD_Ratio)
  
  # Coverage rate data
  dat_cov <- bind_rows(
    prep_study2_subset(results_data, dist_name, "Linear") %>% mutate(Model = "Linear"),
    prep_study2_subset(results_data, dist_name, "Full") %>% mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")),
           Distribution = DIST_LABS[dist_name]) %>%
    transmute(Distribution, Parameter = Parameter_Label, SampleSize, Reliability,
              Method, Model, y = CoverageRate, yerr = 100 * CoverageRate_MCSE)
  
  # Type I error and Power data (special handling for eta4/eta5 parameters)
  eta4_params <- c("eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2")
  eta5_params <- c("eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
  
  # Type I error (Linear model only)
  dat_t1 <- results_data %>%
    filter(Model == "Linear", Distribution == dist_name,
           (Equation == "eta4" & Parameter %in% eta4_params) |
             (Equation == "eta5" & Parameter %in% eta5_params)) %>%
    mutate(
      Parameter = case_when(
        Equation == "eta4" & Parameter == "eta1:eta2" ~ "beta[14]",
        Equation == "eta4" & Parameter == "eta1:eta3" ~ "beta[15]",
        Equation == "eta4" & Parameter == "eta1:eta1" ~ "beta[16]",
        Equation == "eta4" & Parameter == "eta2:eta2" ~ "beta[17]",
        Equation == "eta5" & Parameter == "eta1:eta4" ~ "beta[25]",
        Equation == "eta5" & Parameter == "eta2:eta4" ~ "beta[26]",
        Equation == "eta5" & Parameter == "eta1:eta1" ~ "beta[27]",
        Equation == "eta5" & Parameter == "eta3:eta3" ~ "beta[28]"
      ),
      Method = factor(Method, levels = METHOD_ORDER_5),
      Condition = make_condition(SampleSize, Reliability),
      Distribution = DIST_LABS[dist_name],
      y = TypeI_Error
    ) %>%
    select(Distribution, Parameter, Condition, Method, y)
  
  # Statistical power (Full model only)
  dat_pow <- results_data %>%
    filter(Model == "Full", Distribution == dist_name,
           (Equation == "eta4" & Parameter %in% eta4_params) |
             (Equation == "eta5" & Parameter %in% eta5_params)) %>%
    mutate(
      Parameter = case_when(
        Equation == "eta4" & Parameter == "eta1:eta2" ~ "beta[14]",
        Equation == "eta4" & Parameter == "eta1:eta3" ~ "beta[15]",
        Equation == "eta4" & Parameter == "eta1:eta1" ~ "beta[16]",
        Equation == "eta4" & Parameter == "eta2:eta2" ~ "beta[17]",
        Equation == "eta5" & Parameter == "eta1:eta4" ~ "beta[25]",
        Equation == "eta5" & Parameter == "eta2:eta4" ~ "beta[26]",
        Equation == "eta5" & Parameter == "eta1:eta1" ~ "beta[27]",
        Equation == "eta5" & Parameter == "eta3:eta3" ~ "beta[28]"
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
    # shaded region for acceptable SE/SD ratios (0.85-1.15)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.85, ymax = 1.15,
             fill = "grey93", alpha = .9) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    # truncate extreme values for display
    geom_line(aes(y = pmax(pmin(ratio, 1.6), 0.85)), 
              color = "black", position = position_dodge(width = 0.0)) +
    geom_point(aes(y = pmax(pmin(ratio, 1.6), 0.85)),
               color = "black", size = 2, position = position_dodge(width = 0.0)) +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    coord_cartesian(ylim = c(0.80, 1.4)) + 
    labs(x = "Reliability", y = "SE / SD", title = title) +
    theme_apa_bw()
  
  # labels for extreme values (upper) - within visible range
  p <- p + 
    geom_text(data = filter(data, ratio > 1.15) %>%
                group_by(Distribution, Parameter, Model, SampleSize, Reliability) %>%
                arrange(desc(ratio)) %>%
                mutate(label_y = pmin(1.39 - (row_number() - 1) * 0.05, 1.39)), 
              aes(label = sprintf("%s: %.2f", Method, ratio), y = label_y),
              size = 1.5, hjust = 0, vjust = 1) +
    geom_text(data = filter(data, ratio < 0.85) %>%
                group_by(Distribution, Parameter, Model, SampleSize, Reliability) %>%
                arrange(ratio) %>%
                mutate(label_y = pmax(0.81 + (row_number() - 1) * 0.05, 0.81)), 
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
    # Acceptable coverage range (93-97%)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 93, ymax = 97,
             fill = "grey93", alpha = .9) +
    geom_pointrange(aes(ymin = pmax(y - yerr, 80), 
                        ymax = pmin(y + yerr, 100)),
                    position = position_dodge(width = .6), color = "black") +
    geom_line(aes(y = pmax(y, 80)), position = position_dodge(width = .6), color = "black") +
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

# ============================= USAGE EXAMPLES =============================

# # STUDY 1 USAGE:
# # --------------
# 
# # Prepare all data for Study 1
# study1_data <- prepare_study1_data(results_study_1)
# 
# # Access raw data frames for custom filtering
# dat_bias <- study1_data$bias
# dat_rmse <- study1_data$rmse
# dat_sesd <- study1_data$sesd
# dat_coverage <- study1_data$coverage
# dat_t1 <- study1_data$type1
# dat_pow <- study1_data$power
# 
# # Example: Filter bias data to Full model & N = 400
# dat_bias_400_full <- dat_bias %>%
#   filter(Model == "Full", SampleSize == "400") %>%
#   droplevels()
# 
# # Create custom plot with filtered data
# p_bias_400_full <- plot_bias(
#   dat_bias_400_full, 
#   shapes = SHAPES_4, 
#   ltys = LTYS_4,
#   facet_formula = Distribution + Parameter ~ .,  # Custom faceting
#   y_breaks = seq(-0.20, 0.20, by = 0.05)         # Custom y-axis
# )
# 
# # Or create standard plots
#p_bias_all <- plot_bias(dat_bias, shapes = SHAPES_4, ltys = LTYS_4, 
#                        title = "Study 1: Bias")
#p_rmse_all <- plot_rmse(dat_rmse, shapes = SHAPES_4, ltys = LTYS_4, 
#                         title = "Study 1: RMSE")
# 
# # Get specific values
# sam_bias_normal <- get_value(dat_bias, method = "SAM", dist = "Normal", 
#                              param = "beta[3]", n = "400", rel = "0.6", 
#                              model = "Full")
# 
# # STUDY 2 USAGE:
# # --------------
# 
# # Prepare data for normal distribution
#study2_normal_data <- prepare_study2_data(results_study_2$results, "normal")
# 
# # Access data frames
#dat_bias_normal <- study2_normal_data$bias
#dat_rmse_normal <- study2_normal_data$rmse
# 
# # Filter and create custom plots
#dat_bias_normal_400 <- dat_bias_normal %>%
#   filter(SampleSize == "400")
# 
#p_bias_normal_400 <- plot_bias(
#   dat_bias_normal_400, 
#   shapes = SHAPES_3, 
#   ltys = LTYS_3,
#   title = "Study 2: Normal Distribution (N=400)"
#)