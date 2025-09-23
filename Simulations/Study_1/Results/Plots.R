# Organized but for now include plots together with get_functions for reporting in markdown

# PLOTS

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

# condition factor "N=..., Rel=..."
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

# BIAS
dat_bias <- bind_rows(
  prep_nonlinear(results_with_mcse, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_with_mcse, model = "Full")   %>% mutate(Model = "Full")
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
  prep_nonlinear(results_with_mcse, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_with_mcse, model = "Full")   %>% mutate(Model = "Full")
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






# RMSE
dat_rmse <- bind_rows(
  prep_nonlinear(results_with_mcse, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_with_mcse, model = "Full")   %>% mutate(Model = "Full")
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


## RELATIVE RMSE 
dat_rmse_relative <- bind_rows(
  prep_nonlinear(results_with_mcse, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_with_mcse, model = "Full")   %>% mutate(Model = "Full")
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





# SE/SD RATIO
dat_sesd2 <- bind_rows(
  prep_nonlinear(results_with_mcse, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_with_mcse, model = "Full") %>% mutate(Model = "Full")
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





# COVERAGE 
dat_cov <- bind_rows(
  prep_nonlinear(results_with_mcse, model = "Linear") %>% mutate(Model = "Linear"),
  prep_nonlinear(results_with_mcse, model = "Full")   %>% mutate(Model = "Full")
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





# TYPE I 
dat_t1 <- prep_nonlinear(results_with_mcse, model = "Linear") %>%
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






# POWER 
dat_pow <- prep_nonlinear(results_with_mcse, model = "Full") %>%
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
