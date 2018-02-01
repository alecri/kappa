## Ignoring covariances of log RRs for dose-response meta-analysis of quadratic trend

library(tidyverse)
library(dosresmeta)
library(grid)
library(gridExtra)
library(ggExtra)
library(MASS)
library(parallel)
library(cowplot)

load("data/sim_quadr/param_sim.rda")
k_quadr <- k
file_names <- list.files(path = "data/sim_quadr/", pattern = "dtab([[:digit:]])*.\\.rda") %>%
  sample()
invisible(lapply(file_names, function(f) 
  load(paste0("data/sim_quadr/", f), envir = .GlobalEnv)))

#dose-response model approximating and ignoring the covariance
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type = "FORK")
res_gl_quadr <- parLapply(cl, ls()[grep("dtab", x = ls())], function(f)
  dosresmeta(logrr ~ dose + I(dose^2), id = id, type = type, se = se,
             cases = cases, n = n, data = subset(eval(parse(text = f)), id == 2),
             method = "fixed"))
res_indep_quadr <- parLapply(cl, ls()[grep("dtab", x = ls())], function(f)
  dosresmeta(logrr ~ dose + I(dose^2), id = id, type = type, se = se,
             cases = cases, n = n, data = subset(eval(parse(text = f)), id == 2),
             method = "fixed", covariance = "indep"))
stopCluster(cl)
rm(list = ls()[grep("dtab", x = ls())])
beta_quadr <- beta
save(beta_quadr, res_gl_quadr, k_quadr, res_indep_quadr,
     file = "data/sim_quadr/res_cov_quadr.RData")


## Results
load("data/sim_quadr/res_cov_quadr.RData")

# predicted curves
xref <- 1
xmax <- 10
pred <- data_frame(
  dose = seq(xref, xmax, length.out = 99),
  true = exp(beta_quadr[2]*(dose - xref) + beta_quadr[3]*(dose^2 - xref^2))
) %>%
  bind_cols(
    lapply(res_gl_quadr, function(m) 
      data.frame(`gl_` = predict(m, newdata = ., expo = T)$pred)),
    lapply(res_indep_quadr, function(m) 
      data.frame(`indep_` = predict(m, newdata = ., expo = T)$pred))
  ) %>% 
  gather(sim, pred, -dose, -true) %>%
  separate(sim, c("cov", "sim"), extra = "drop", fill = "right") %>%
  mutate(sim = replace(sim, sim == "", 0), sim = as.numeric(sim) + 1) %>%
  arrange(sim)

ggplot(pred, aes(dose, pred, group = sim)) + 
  geom_line(col = "grey") +
  geom_line(aes(y = true), lwd = 1.5) +
  scale_y_continuous(trans = "log", breaks = c(.25, .5, .75, 1, 2, 4)) +
  labs(x = "Dose", y = "Odds ratio") +
  facet_grid(~ cov, labeller = as_labeller(c("gl" = "Greenland-Longnecker", 
                                             "indep" = "Independence"))) +
  theme_classic()

# comparing distribution of betas
betas_quadr <- data.frame(
  b.gl = do.call("rbind", lapply(res_gl_quadr, coef)),
  b.indep = do.call("rbind", lapply(res_indep_quadr, coef))
) %>%
  gather(sim, beta) %>%
  separate(sim, c("b", "cov", "coef", "int"), extra = "drop", fill = "right") %>%
  mutate(sim = rep(1:k_quadr, 4))

betas_quadr$xintercept <- beta_quadr[2]*(betas_quadr$coef == "dose") + 
  beta_quadr[3]*(betas_quadr$coef == "I")
ggplot(betas_quadr, aes(x = beta, col = cov)) +
  geom_line(stat = "density") +
  facet_wrap(~ coef, scales = "free", 
             labeller = as_labeller(c(dose = "beta[1]", I = "beta[2]"), label_parsed)) +
  labs(x = "", col = "Method") +
  geom_vline(aes(xintercept = xintercept), linetype = "dashed")

betas_quadr %>%
  dplyr::select(coef, beta, sim, cov) %>%
  spread(coef, beta) %>%
  ggplot(aes(x = dose, y = I, group = cov)) +
  stat_density2d(aes(fill = ..level..), geom = "polygon") +
  #geom_density2d() +
  labs(x = expression(beta[1]), y = expression(beta[2])) +
  facet_grid(~ cov) +
  theme_classic()

betas_quadr %>%
  dplyr::select(coef, beta, sim, cov) %>%
  spread(coef, beta) %>%
  ggplot(aes(x = dose, y = I, col = cov)) +
  #stat_density2d(aes(fill = ..level..), geom = "polygon") +
  geom_density2d() +
  labs(x = expression(beta[1]), y = expression(beta[2])) +
  scale_color_discrete(name = "Method") +
  theme_classic()
# ggsave("../figures/comp_beta_densities.pdf")


# adding the marginal distribution to the plot
p_gl_quadr <- betas_quadr %>%
  filter(cov == "gl") %>%
  dplyr::select(coef, beta, sim) %>%
  spread(coef, beta) %>%
  ggplot(aes(x = dose, y = I)) +
  #stat_density2d(aes(fill = ..level..), geom = "polygon") +
  geom_density2d() +
  #xlim(c(-.2, 0)) + #ylim(c(-0.005, 0.01)) + 
  theme_classic()
ggMarginal(p_gl_quadr)

p_indep_quadr <- betas_quadr %>%
  filter(cov == "indep") %>%
  dplyr::select(coef, beta, sim) %>%
  spread(coef, beta) %>%
  ggplot(aes(x = dose, y = I)) +
  #stat_density2d(aes(fill = ..level..), geom = "polygon") +
  geom_density2d() +
  #xlim(c(-.2, 0)) + #ylim(c(-0.005, 0.01)) + 
  theme_classic()
ggMarginal(p_indep_quadr)

grid.arrange(ggMarginal(p_gl_quadr), ggMarginal(p_indep_quadr), ncol = 2)


# comparing standard error (or variances)
bvcov_quadr <- data.frame(
  v.gl = do.call("rbind", lapply(res_gl_quadr, function(r) c(vcov(r))[-2])),
  v.indep = do.call("rbind", lapply(res_indep_quadr, function(r) c(vcov(r))[-2]))
) %>%
  gather(sim, vcov) %>%
  separate(sim, c("v", "cov", "coef")) %>% 
  mutate(sim = rep(1:k, 6))

ggplot(bvcov_quadr, aes(x = vcov, col = cov)) +
  geom_line(stat = "density") +
  facet_wrap(~ coef, scales = "free") +
  labs(x = "", col = "Method") +
  theme_classic()

betas_quadr %>% 
  mutate(coef = 1*(coef == "dose") + 3*(coef == "I")) %>%
  group_by(cov, coef) %>%
  summarise(var(beta)) %>%
  merge(
    bvcov_quadr %>% group_by(cov, coef) %>%
      summarise(mean(vcov)), all = T
  )

# save(betas_quadr, bvcov_quadr, beta_quadr, k_quadr, file = "data/res_cov_quadr.RData")
# load("data/res_cov_quadr.RData")

p_betas_quadr <- betas_quadr %>%
  dplyr::select(coef, beta, sim, cov) %>%
  spread(coef, beta) %>%
  ggplot(aes(x = dose, y = I, col = cov)) +
  #stat_density2d(aes(fill = ..level..), geom = "polygon") +
  geom_density2d() +
  geom_vline(xintercept = beta_quadr[2], linetype = "dotted") +
  geom_hline(yintercept = beta_quadr[3], linetype = "dotted") +
  scale_color_discrete(labels = c("Greenland-Longnecker", "Independence")) +
  labs(x = expression(hat(beta)[1]), y = expression(hat(beta)[2]), col = "Method") +
  theme_classic()
p_se_quadr <- bvcov_quadr %>%
  dplyr::select(coef, vcov, sim, cov) %>%
  spread(coef, vcov) %>%
  ggplot(aes(x = `1`^.5, y = `3`^.5, col = cov)) +
  geom_density2d() +
  scale_color_discrete(labels = c("Greenland-Longnecker", "Independence")) +
  labs(x = expression(paste("SE(", hat(beta)[1], ")")), col = "Method", 
       y = expression(paste("SE(", hat(beta)[2], ")"))) +
  theme_classic()
p_sim_quadr <- plot_grid( 
  p_betas_quadr + theme(legend.position="none"),
  p_se_quadr + theme(legend.position="none"),
  align = 'vh', labels = c("A", "B"), hjust = -1, nrow = 1)
legend_sim_quadr <- get_legend(p_betas_quadr + theme(legend.position="bottom"))
p_sim_quadr <- plot_grid(p_sim_quadr, legend_sim_quadr, ncol = 1, rel_heights = c(1, .1))
p_sim_quadr
