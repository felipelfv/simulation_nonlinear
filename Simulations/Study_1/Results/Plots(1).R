############################ 1. General Information ############################

# Organized but for now include plots together with get_functions for reporting in markdown

############################### 2. Documentation ################################

#' Plotting and Reporting Functions for Simulation Study Results
#' 
#' @section Global Settings:
#' @param METHOD_ORDER   Character vector. Display order for methods: c("SAM","LMS","QML","DBLCENT")
#' @param DIST_LABS      Named vector. Labels for distributions: normal="Normal", nonnormal="Nonnormal", uniform="Uniform"
#' @param PARAM_ORDER    Character vector. Display order for parameters: c("eta1:eta2","eta1:eta1","eta2:eta2")
#' @param PARAM_LABS     Named vector. LaTeX labels for parameters in plots
#' @param GREYS_METHOD   Named vector. Grayscale colors for each method (SAM darkest to DBLCENT lightest)
#' @param GREYS_PARAM    Named vector. Grayscale colors for each parameter
#' 
#' @section Plotting Theme:
#' theme_apa_bw()
#' @param base_size      Numeric. Base font size for plots (default = 11)
#' @return ggplot2 theme object with APA-style formatting (white background, minimal grid, bold facet labels)
#' 
#' @section Helper Functions:
#' make_condition()
#' @param N              Numeric. Sample size (400 or 1000)
#' @param Rel            Numeric. Reliability level (0.4, 0.6, or 0.8)
#' @return Factor with formatted condition labels "N=..., Rel=..."
#' 
#' prep_nonlinear()
#' @param df             Data.frame. Results data containing performance metrics
#' @param model          Character. "Full" or "Linear" model type
#' @return Tibble with filtered nonlinear parameters and formatted factors for plotting
#' 
#' @section Functions:
#' 
#' get_absolute_bias()
#' @param method         Character. Method name ("SAM", "LMS", "QML", or "DBLCENT")
#' @param dist           Character. Distribution type ("Normal", "Nonnormal", or "Uniform")
#' @param param          Character. Parameter ("beta[3]", "beta[4]", or "beta[5]")
#' @param n              Numeric. Sample size (400 or 1000)
#' @param rel            Numeric. Reliability (0.4, 0.6, or 0.8)
#' @param model          Character. Model type ("Full" or "Linear")
#' @return Numeric. Absolute bias value (Mean estimate - True value)
#' 
#' get_bias()
#' @param method, dist, param, n, rel, model  Same as get_absolute_bias()
#' @return Numeric. Relative bias value ((Mean estimate / True value) - 1)
#' 
#' get_rmse()
#' @param method, dist, param, n, rel, model  Same as get_absolute_bias()
#' @return Numeric. Root mean squared error
#' 
#' get_relative_rmse()
#' @param method, dist, param, n, rel, model  Same as get_absolute_bias()
#' @return Numeric. Relative root mean squared error
#' 
#' get_ratio()
#' @param method, dist, param, n, rel, model  Same as get_absolute_bias()
#' @return Numeric. SE/SD ratio (Mean SE / SD of estimates)
#' 
#' get_coverage()
#' @param method, dist, param, n, rel, model  Same as get_absolute_bias()
#' @return Numeric. Coverage rate percentage
#' 
#' get_type()
#' @param method         Character. Method name
#' @param dist           Character. Distribution type
#' @param param          Character. Parameter
#' @param cond           Character. Condition label "N=..., Rel=..."
#' @return Numeric. Type I error rate percentage (for Linear model)
#' 
#' get_power()
#' @param method, dist, param, cond  Same as get_type()
#' @return Numeric. Statistical power percentage (for Full model)
#' 
#' @section Plot Objects:
#' 
#' p_bias           - Absolute bias plot with error bars (MCSE)
#' p_rmse           - RMSE plot with error bars (MCSE)  
#' p_slope          - SE/SD ratio plot with truncation at [0.85, 1.6] and text labels for outliers
#' p_cov            - Coverage rate plot with 95% confidence band and text labels for <80%
#' p_t1             - Type I error bar chart with 5% reference line (Linear model only)
#' p_pow            - Statistical power bar chart with 80% reference line (Full model only)

############################### 3. PLOTS ########################################

# orderings and labels
METHOD_ORDER <- c("SAM","LMS","QML","DBLCENT")
DIST_LABS    <- c(normal="Normal", nonnormal="Nonnormal", uniform="Uniform")
PARAM_ORDER  <- c("eta1:eta2","eta1:eta1","eta2:eta2")
PARAM_LABS   <- c("eta1:eta2"="beta[3]",
                  "eta1:eta1"="beta[4]",
                  "eta2:eta2"="beta[5]")

GREYS_METHOD <- c(SAM="grey20", LMS="grey45", QML="grey65", DBLCENT="grey85")
GREYS_PARAM  <- c("eta1:eta2"="grey20", "eta1:eta1"="grey55", "eta2:eta2"="grey85")

# APA-based
theme_apa_bw <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

# condition factor 
make_condition <- function(N, Rel) {
  factor(paste0("N=", N, ", Rel=", Rel),
         levels = c("N=400, Rel=0.4","N=400, Rel=0.6","N=400, Rel=0.8",
                    "N=1000, Rel=0.4","N=1000, Rel=0.6","N=1000, Rel=0.8"))
}

# data prep for nonlinear parameters only 
prep_nonlinear <- function(df, model = c("Full","Linear")) {
  model <- match.arg(model)
  df %>%
    filter(Parameter %in% PARAM_ORDER, Model == model) %>%
    mutate(
      Method       = factor(Method, levels = METHOD_ORDER),
      Distribution = factor(Distribution, levels = names(DIST_LABS), labels = DIST_LABS),
      Parameter    = factor(Parameter, levels = PARAM_ORDER, 
                            labels = c("beta[3]", "beta[4]", "beta[5]")),  
      Reliability  = factor(Reliability),
      SampleSize   = factor(SampleSize),
      Condition    = make_condition(SampleSize, Reliability)
    )
}

shapes <- c(SAM=16, LMS=17, QML=15, DBLCENT=18)
ltys   <- c(SAM="solid", LMS="dashed", QML="dotdash", DBLCENT="twodash")

# BIAS functions and plots
dat_bias <- bind_rows(
  prep_nonlinear(results_study_1, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_study_1, model = "Full")   %>% mutate(Model = "Full")
) %>%
  mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
  transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
            y = Bias_Mean, yerr = Bias_Mean_MCSE)

get_absolute_bias <- function(method, dist, param, n, rel, model = c("Full","Linear")) {
  model <- match.arg(model)
  dat_bias %>%
    filter(Method == method,
           Distribution == dist,
           Parameter == param,
           SampleSize == n,
           Reliability == rel,
           Model == model) %>%
    select(y)
}

p_bias <- ggplot(
  dat_bias,
  aes(x = Reliability, y = y, shape = Method, linetype = Method, group = Method)
) +
  geom_pointrange(aes(ymin = y - yerr, ymax = y + yerr),
                  position = position_dodge(width = .6), color = "black") +
  geom_line(position = position_dodge(width = .6), color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = ltys) +
  facet_grid(Distribution + Parameter ~ Model + SampleSize,
             labeller = labeller(SampleSize = label_value,
                                 Parameter = label_parsed)) +
  labs(x = "Reliability", y = "Bias", title = "") +
  theme_apa_bw()

dat_bias_relative <- bind_rows(
  prep_nonlinear(results_study_1, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_study_1, model = "Full")   %>% mutate(Model = "Full")
) %>%
  mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
  transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
            y = RelativeBias_Mean, yerr = RelativeBias_MCSE)  # MCSE to percentage

get_bias <- function(method, dist, param, n, rel, model = c("Full","Linear")) {
  model <- match.arg(model)
  dat_bias_relative %>%
    filter(Method == method,
           Distribution == dist,
           Parameter == param,
           SampleSize == n,
           Reliability == rel,
           Model == model) %>%
    select(y)
}


# RMSE functions and plots
dat_rmse <- bind_rows(
  prep_nonlinear(results_study_1, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_study_1, model = "Full")   %>% mutate(Model = "Full")
) %>%
  mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
  transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
            y    = RMSE_Mean,
            yerr = RMSE_Mean_MCSE)  

p_rmse <- ggplot(
  dat_rmse,
  aes(x = Reliability, y = y, shape = Method, linetype = Method, group = Method)
) +
  geom_pointrange(aes(ymin = y - yerr, ymax = y + yerr),
                  position = position_dodge(width = .6), color = "black") +
  geom_line(position = position_dodge(width = .6), color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = ltys) +
  facet_grid(Distribution + Parameter ~ Model + SampleSize,
             labeller = labeller(SampleSize = label_value,
                                 Parameter = label_parsed)) +
  labs(x = "Reliability", y = "RMSE", title = "") +
  theme_apa_bw()

get_rmse <- function(method, dist, param, n, rel, model = c("Full","Linear")) {
  model <- match.arg(model)
  dat_rmse %>%
    filter(Method == method,
           Distribution == dist,
           Parameter == param,
           SampleSize == n,
           Reliability == rel,
           Model == model) %>%
    select(y)
}


# RELATIVE RMSE functions and plots
dat_rmse_relative <- bind_rows(
  prep_nonlinear(results_study_1, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_study_1, model = "Full")   %>% mutate(Model = "Full")
) %>%
  mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
  transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
            y    = Relative_RMSE,
            yerr = Relative_RMSE_MCSE)

get_relative_rmse <- function(method, dist, param, n, rel, model = c("Full","Linear")) {
  model <- match.arg(model)
  dat_rmse_relative %>%
    filter(Method == method,
           Distribution == dist,
           Parameter == param,
           SampleSize == n,
           Reliability == rel,
           Model == model) %>%
    select(y)
}

# SE/SD RATIO functions and plots
dat_sesd2 <- bind_rows(
  prep_nonlinear(results_study_1, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_study_1, model = "Full") %>% mutate(Model = "Full")
) %>%
  mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
  transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
            ratio = SE_SD_Ratio)

shapes <- c(SAM=16, LMS=17, QML=15, DBLCENT=18)
ltys <- c(SAM="solid", LMS="dashed", QML="dotdash", DBLCENT="twodash")

p_slope <- ggplot(dat_sesd2,
                  aes(x = Reliability, y = ratio, group = Method,
                      shape = Method, linetype = Method)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(aes(y = pmax(pmin(ratio, 1.6), 0.85)), 
            color = "black", position = position_dodge(width = 0.0)) +
  geom_point(aes(y = pmax(pmin(ratio, 1.6), 0.85)),
             color = "black", size = 2, position = position_dodge(width = 0.0)) +
  # text labels for values above 1.6
  geom_text(data = filter(dat_sesd2, ratio > 1.6) %>%
              group_by(Distribution, Parameter, Model, SampleSize, Reliability) %>%
              arrange(desc(ratio)) %>%
              mutate(label_y = 1.59 - (row_number() - 1) * 0.8), 
            aes(label = sprintf("%s: %.2f", Method, ratio),
                y = label_y),
            size = 1.5, hjust = 0, vjust = 1) +
  # text labels for values below 0.85
  geom_text(data = filter(dat_sesd2, ratio < 0.85) %>%
              group_by(Distribution, Parameter, Model, SampleSize, Reliability) %>%
              arrange(ratio) %>%
              mutate(label_y = 0.86 + (row_number() - 1) * 0.8), 
            aes(label = sprintf("%s: %.2f", Method, ratio),
                y = label_y),
            size = 1.5, hjust = 0, vjust = 0) +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = ltys) +
  facet_grid(Distribution + Parameter ~ Model + SampleSize,
             labeller = labeller(SampleSize = label_value,
                                 Parameter = label_parsed)) +
  coord_cartesian(ylim = c(0.85, 1.6)) + 
  labs(x = "Reliability", y = "SE / SD",
       title = "") +
  theme_apa_bw()

get_ratio <- function(method, dist, param, n, rel, model = c("Full", "Linear")) {
  model <- match.arg(model)
  dat_sesd2 %>%
    filter(Method == method,
           Distribution == dist,
           Parameter == param,
           SampleSize == n,
           Reliability == rel,
           Model == model) %>%
    pull(ratio)
}


# COVERAGE functions and plots
dat_cov <- bind_rows(
  prep_nonlinear(results_study_1, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_study_1, model = "Full")   %>% mutate(Model = "Full")
) %>%
  mutate(Model = factor(Model, levels = c("Linear", "Full"))) %>%
  transmute(Distribution, Parameter, SampleSize, Reliability, Method, Model,
            y = CoverageRate, yerr = 100*CoverageRate_MCSE)

p_cov <- ggplot(
  dat_cov,
  aes(x = Reliability, y = y, shape = Method, linetype = Method, group = Method)
) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 93, ymax = 97,
           fill = "grey93", alpha = .9) +
  geom_pointrange(aes(ymin = pmax(y - yerr, 80), 
                      ymax = pmin(y + yerr, 100)),
                  position = position_dodge(width = .6), color = "black") +
  geom_line(aes(y = pmax(y, 80)), position = position_dodge(width = .6), color = "black") +
  # text labels for values below 80%
  geom_text(data = filter(dat_cov, y < 80) %>%
              group_by(Distribution, Parameter, Model, SampleSize, Reliability) %>%
              arrange(Method) %>%
              mutate(label_y = 81 + (row_number() - 1) * 2.5), # labels vertically
            aes(label = sprintf("%s: %.0f%%", Method, y),
                y = label_y),
            size = 1.5, hjust = 0, vjust = 0) +
  geom_hline(yintercept = 95, linetype = "dotted") +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = ltys) +
  facet_grid(Distribution + Parameter ~ Model + SampleSize,
             labeller = labeller(SampleSize = label_value,
                                 Parameter = label_parsed)) +
  coord_cartesian(ylim = c(80, 100)) +
  labs(x = "Reliability", y = "Coverage (%)", title = "") +
  theme_apa_bw()

get_coverage <- function(method, dist, param, n, rel, model = c("Full", "Linear")) {
  model <- match.arg(model)
  dat_cov %>%
    filter(Method == method,
           Distribution == dist,
           Parameter == param,
           SampleSize == n,
           Reliability == rel,
           Model == model) %>%
    pull(y)
}

# TYPE I functions and plots
dat_t1 <- prep_nonlinear(results_study_1, model = "Linear") %>%
  transmute(Distribution, Parameter, Condition, Method,
            y = TypeI_Error)

ymax_t1 <- max(dat_t1$y, na.rm = TRUE) * 1.15  # 15% headroom for labels

p_t1 <- ggplot(dat_t1, aes(x = Method, y = y, fill = Method)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3, ymax = 7,
           fill = "grey92", alpha = .9) +
  geom_col(position = position_dodge(width = 0.9), width = 0.75, color = "black") +
  geom_text(aes(label = sprintf("%.1f", y)), vjust = -0.4, size = 3) +
  geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.3) +
  facet_grid(Distribution + Parameter ~ Condition,
             labeller = labeller(Parameter = label_parsed)) +
  scale_fill_manual(values = GREYS_METHOD) +
  coord_cartesian(ylim = c(0, ymax_t1)) +
  labs(x = NULL, y = "Type I error (%)",
       title = "") +
  theme_apa_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())


get_type <- function(method, dist, param, cond) {
  dat_t1 %>%
    filter(Method == method,
           Distribution == dist,
           Parameter == param,
           Condition == cond) %>%
    select(y)
}

# POWER functions and plots
dat_pow <- prep_nonlinear(results_study_1, model = "Full") %>%
  transmute(Distribution, Parameter, Condition, Method,
            y = Power)

# y-axis to go to 110 to accommodate labels at 100%
ymax_pow <- 110  

p_pow <- ggplot(dat_pow, aes(x = Method, y = y, fill = Method)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.75, color = "black") +
  geom_text(aes(label = sprintf("%.0f", y)), vjust = -0.4, size = 3) +
  geom_hline(yintercept = 80, linetype = "dotted", linewidth = 0.3) +
  facet_grid(Distribution + Parameter ~ Condition,
             labeller = labeller(Parameter = label_parsed)) +
  scale_fill_manual(values = GREYS_METHOD) +
  coord_cartesian(ylim = c(0, ymax_pow)) +
  labs(x = NULL, y = "Power (%)",
       title = "") +
  theme_apa_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())

get_power <- function(method, dist, param, cond) {
  dat_pow %>%
    filter(Method == method,
           Distribution == dist,
           Parameter == param,
           Condition == cond) %>%
    select(y)
}
