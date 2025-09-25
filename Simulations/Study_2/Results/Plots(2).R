############################### Study 2 plots by distribution ########################################

library(tidyverse); library(ggplot2)

#load("Simulations/Study_2/Data/Results_Study_2_Performance.RData")

# similar settings as stdy 1

METHOD_ORDER <- c("SAM","QML","DBLCENT")
DIST_LABS <- c(normal="Normal", nonnormal="Nonnormal", uniform="Uniform")

# nonlinear parameters with equation labels
PARAM_ORDER_ALL <- c("eta1:eta2", "eta1:eta1",  
                     "eta2:eta4", "eta2:eta2",  
                     "eta3:eta5", "eta3:eta3")

PARAM_LABS_ALL <- c(
  "eta1:eta2" = "η1:η2→η4",
  "eta1:eta1" = "η1²→η4",
  "eta2:eta4" = "η2:η4→η5",
  "eta2:eta2" = "η2²→η5",
  "eta3:eta5" = "η3:η5→η6",
  "eta3:eta3" = "η3²→η6"
)

shapes <- c(SAM=16, QML=15, DBLCENT=18)
ltys <- c(SAM="solid", QML="dotdash", DBLCENT="twodash")
GREYS_METHOD <- c(SAM="grey20", QML="grey50", DBLCENT="grey80")

# HELPER FUNCTIONS

# theme function
theme_apa_bw <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

# condition factor helper
make_condition <- function(N, Rel) {
  factor(paste0("N=", N, ", Rel=", Rel),
         levels = c("N=400, Rel=0.4","N=400, Rel=0.6","N=400, Rel=0.8",
                    "N=1000, Rel=0.4","N=1000, Rel=0.6","N=1000, Rel=0.8"))
}

# data prep function for a specific distribution
prep_nonlinear_by_dist <- function(df, distribution, model = c("Full","Linear")) {
  model <- match.arg(model)
  
  eta4_params <- c("eta1:eta2", "eta1:eta1")
  eta5_params <- c("eta2:eta4", "eta2:eta2")
  eta6_params <- c("eta3:eta5", "eta3:eta3")
  
  df %>%
    filter(Distribution == distribution,
           Model == model,
           (Equation == "eta4" & Parameter %in% eta4_params) |
             (Equation == "eta5" & Parameter %in% eta5_params) |
             (Equation == "eta6" & Parameter %in% eta6_params)) %>%
    mutate(
      Method = factor(Method, levels = METHOD_ORDER),
      Parameter = factor(Parameter, levels = PARAM_ORDER_ALL, labels = PARAM_LABS_ALL),
      Reliability = factor(Reliability),
      SampleSize = factor(SampleSize)
    )
}

# ===================== MAIN PLOTTING FUNCTION =====================

create_distribution_plots <- function(dist_name) {
  
  # BIAS PLOT 
  dat_bias <- bind_rows(
    prep_nonlinear_by_dist(results_study_2$results, dist_name, "Linear") %>% 
      mutate(Model = "Linear"),
    prep_nonlinear_by_dist(results_study_2$results, dist_name, "Full") %>% 
      mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")))
  
  p_bias <- ggplot(dat_bias,
                   aes(x = Reliability, y = Bias_Mean, 
                       shape = Method, linetype = Method, group = Method)) +
    geom_pointrange(aes(ymin = Bias_Mean - Bias_Mean_MCSE, 
                        ymax = Bias_Mean + Bias_Mean_MCSE),
                    position = position_dodge(width = .4), color = "black") +
    geom_line(position = position_dodge(width = .4), color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    facet_grid(Parameter ~ Model + SampleSize,
               labeller = labeller(SampleSize = label_value)) +
    labs(x = "Reliability", y = "Bias", 
         title = paste("Bias -", DIST_LABS[dist_name], "Distribution")) +
    theme_apa_bw()
  
  # RMSE PLOT 
  dat_rmse <- bind_rows(
    prep_nonlinear_by_dist(results_study_2$results, dist_name, "Linear") %>% 
      mutate(Model = "Linear"),
    prep_nonlinear_by_dist(results_study_2$results, dist_name, "Full") %>% 
      mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")))
  
  p_rmse <- ggplot(dat_rmse,
                   aes(x = Reliability, y = RMSE_Mean, 
                       shape = Method, linetype = Method, group = Method)) +
    geom_pointrange(aes(ymin = RMSE_Mean - RMSE_Mean_MCSE, 
                        ymax = RMSE_Mean + RMSE_Mean_MCSE),
                    position = position_dodge(width = .4), color = "black") +
    geom_line(position = position_dodge(width = .4), color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    facet_grid(Parameter ~ Model + SampleSize,
               labeller = labeller(SampleSize = label_value)) +
    labs(x = "Reliability", y = "RMSE", 
         title = paste("RMSE -", DIST_LABS[dist_name], "Distribution")) +
    theme_apa_bw()
  
  # SE/SD RATIO PLOT
  dat_sesd <- bind_rows(
    prep_nonlinear_by_dist(results_study_2$results, dist_name, "Linear") %>% 
      mutate(Model = "Linear"),
    prep_nonlinear_by_dist(results_study_2$results, dist_name, "Full") %>% 
      mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")))
  
  p_sesd <- ggplot(dat_sesd,
                   aes(x = Reliability, y = SE_SD_Ratio, 
                       group = Method, shape = Method, linetype = Method)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_line(aes(y = pmax(pmin(SE_SD_Ratio, 1.6), 0.85)), 
              color = "black", position = position_dodge(width = .4)) +
    geom_point(aes(y = pmax(pmin(SE_SD_Ratio, 1.6), 0.85)),
               color = "black", size = 2, position = position_dodge(width = .4)) +
    geom_text(data = filter(dat_sesd, SE_SD_Ratio > 1.6) %>%
                group_by(Parameter, Model, SampleSize, Reliability) %>%
                arrange(desc(SE_SD_Ratio)) %>%
                mutate(label_y = 1.59 - (row_number() - 1) * 0.08), 
              aes(label = sprintf("%s: %.2f", Method, SE_SD_Ratio),
                  y = label_y),
              size = 2, hjust = 0, vjust = 1) +
    geom_text(data = filter(dat_sesd, SE_SD_Ratio < 0.85) %>%
                group_by(Parameter, Model, SampleSize, Reliability) %>%
                arrange(SE_SD_Ratio) %>%
                mutate(label_y = 0.86 + (row_number() - 1) * 0.08), 
              aes(label = sprintf("%s: %.2f", Method, SE_SD_Ratio),
                  y = label_y),
              size = 2, hjust = 0, vjust = 0) +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    facet_grid(Parameter ~ Model + SampleSize,
               labeller = labeller(SampleSize = label_value)) +
    coord_cartesian(ylim = c(0.85, 1.6)) + 
    labs(x = "Reliability", y = "SE / SD",
         title = paste("SE/SD Ratio -", DIST_LABS[dist_name], "Distribution")) +
    theme_apa_bw()
  
  # COVERAGE PLOT 
  dat_cov <- bind_rows(
    prep_nonlinear_by_dist(results_study_2$results, dist_name, "Linear") %>% 
      mutate(Model = "Linear"),
    prep_nonlinear_by_dist(results_study_2$results, dist_name, "Full") %>% 
      mutate(Model = "Full")
  ) %>%
    mutate(Model = factor(Model, levels = c("Linear", "Full")))
  
  p_cov <- ggplot(dat_cov,
                  aes(x = Reliability, y = CoverageRate, 
                      shape = Method, linetype = Method, group = Method)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 93, ymax = 97,
             fill = "grey93", alpha = .9) +
    geom_pointrange(aes(ymin = pmax(CoverageRate - 100*CoverageRate_MCSE, 80), 
                        ymax = pmin(CoverageRate + 100*CoverageRate_MCSE, 100)),
                    position = position_dodge(width = .4), color = "black") +
    geom_line(aes(y = pmax(CoverageRate, 80)), 
              position = position_dodge(width = .4), color = "black") +
    geom_text(data = filter(dat_cov, CoverageRate < 80) %>%
                group_by(Parameter, Model, SampleSize, Reliability) %>%
                arrange(Method) %>%
                mutate(label_y = 81 + (row_number() - 1) * 2.5),
              aes(label = sprintf("%s: %.0f%%", Method, CoverageRate),
                  y = label_y),
              size = 2, hjust = 0, vjust = 0) +
    geom_hline(yintercept = 95, linetype = "dotted") +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = ltys) +
    facet_grid(Parameter ~ Model + SampleSize,
               labeller = labeller(SampleSize = label_value)) +
    coord_cartesian(ylim = c(80, 100)) +
    labs(x = "Reliability", y = "Coverage (%)", 
         title = paste("Coverage -", DIST_LABS[dist_name], "Distribution")) +
    theme_apa_bw()
  
  # TYPE I ERROR PLOT 
  dat_t1 <- results_study_2$results %>%
    filter(Model == "Linear",
           Distribution == dist_name,
           (Equation == "eta4" & Parameter %in% c("eta1:eta2", "eta1:eta1")) |
             (Equation == "eta5" & Parameter %in% c("eta2:eta4", "eta2:eta2")) |
             (Equation == "eta6" & Parameter %in% c("eta3:eta5", "eta3:eta3"))) %>%
    mutate(
      Method = factor(Method, levels = METHOD_ORDER),
      Parameter = factor(Parameter, levels = PARAM_ORDER_ALL, labels = PARAM_LABS_ALL),
      Condition = make_condition(SampleSize, Reliability)
    )
  
  ymax_t1 <- max(dat_t1$TypeI_Error, na.rm = TRUE) * 1.15
  
  p_type1 <- ggplot(dat_t1, aes(x = Method, y = TypeI_Error, fill = Method)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3, ymax = 7,
             fill = "grey92", alpha = .9) +
    geom_col(position = position_dodge(width = 0.9), width = 0.75, color = "black") +
    geom_text(aes(label = sprintf("%.1f", TypeI_Error)), vjust = -0.4, size = 2.5) +
    geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.3) +
    facet_grid(Parameter ~ Condition) +
    scale_fill_manual(values = GREYS_METHOD) +
    coord_cartesian(ylim = c(0, ymax_t1)) +
    labs(x = NULL, y = "Type I error (%)",
         title = paste("Type I Error -", DIST_LABS[dist_name], "Distribution")) +
    theme_apa_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # POWER PLOT 
  dat_pow <- results_study_2$results %>%
    filter(Model == "Full",
           Distribution == dist_name,
           (Equation == "eta4" & Parameter %in% c("eta1:eta2", "eta1:eta1")) |
             (Equation == "eta5" & Parameter %in% c("eta2:eta4", "eta2:eta2")) |
             (Equation == "eta6" & Parameter %in% c("eta3:eta5", "eta3:eta3"))) %>%
    mutate(
      Method = factor(Method, levels = METHOD_ORDER),
      Parameter = factor(Parameter, levels = PARAM_ORDER_ALL, labels = PARAM_LABS_ALL),
      Condition = make_condition(SampleSize, Reliability)
    )
  
  p_power <- ggplot(dat_pow, aes(x = Method, y = Power, fill = Method)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.75, color = "black") +
    geom_text(aes(label = sprintf("%.0f", Power)), vjust = -0.4, size = 2.5) +
    geom_hline(yintercept = 80, linetype = "dotted", linewidth = 0.3) +
    facet_grid(Parameter ~ Condition) +
    scale_fill_manual(values = GREYS_METHOD) +
    coord_cartesian(ylim = c(0, 110)) +
    labs(x = NULL, y = "Power (%)",
         title = paste("Power -", DIST_LABS[dist_name], "Distribution")) +
    theme_apa_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # list of all plots
  list(
    bias = p_bias,
    rmse = p_rmse,
    sesd = p_sesd,
    coverage = p_cov,
    type1 = p_type1,
    power = p_power
  )
}