## Ignoring covariances of log RRs for dose-response meta-analysis of linear trend

library(tidyverse)
library(dosresmeta)
library(grid)
library(gridExtra)
library(ggExtra)
library(MASS)
library(parallel)
library(cowplot)

# Loading simulated data
load("data/sim_lin/param_sim.rda")
k_lin <- k 
file_names <- list.files(path = "data/sim_lin/", pattern = "dtab([[:digit:]])*.\\.rda") %>%
  sample()
invisible(lapply(file_names, function(f) 
  load(paste0("data/sim_lin/", f), envir = .GlobalEnv)))


#dose-response model approximating and ignoring the covariance
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type = "FORK")
res_gl_lin <- parLapply(cl, ls()[grep("dtab", x = ls())], function(f)
  dosresmeta(logrr ~ dose, id = id, type = type, se = se,
             cases = cases, n = n, data = subset(eval(parse(text = f)), id == 1),
             method = "fixed"))
res_indep_lin <- parLapply(cl, ls()[grep("dtab", x = ls())], function(f)
  dosresmeta(logrr ~ dose, id = id, type = type, se = se,
             cases = cases, n = n, data = subset(eval(parse(text = f)), id == 1),
             method = "fixed", covariance = "indep"))
p_val_gl_lin <- parSapply(cl, ls()[grep("dtab", x = ls())], function(f)
  summary(
    dosresmeta(logrr ~ dose + I(dose^2), id = id, type = type, se = se,
               cases = cases, n = n, data = subset(eval(parse(text = f)), id == 1),
               method = "fixed")
  )$coefficients[2, "Pr(>|z|)"])
p_val_indep_lin <- parSapply(cl, ls()[grep("dtab", x = ls())], function(f)
  summary(
    dosresmeta(logrr ~ dose + I(dose^2), id = id, type = type, se = se,
               cases = cases, n = n, data = subset(eval(parse(text = f)), id == 1),
               method = "fixed", covariance = "indep")
  )$coefficients[2, "Pr(>|z|)"])
stopCluster(cl)
rm(list = ls()[grep("dtab", x = ls())])
beta_lin <- beta
save(beta_lin, res_gl_lin, k_lin = k, res_indep_lin, p_val_gl_lin, p_val_indep_lin,
     file = "data/sim_lin/res_cov_lin.RData")



## Results
load("data/sim_lin/res_cov_lin.RData")

# predicted curves
xref <- 1
xmax <- 10
pred <- data_frame(
  dose = seq(xref, xmax, length.out = 99),
  true = exp(beta[2]*(dose - xref))
) %>%
  bind_cols(
    lapply(res_gl_lin, function(m) 
      data.frame(`gl_` = predict(m, newdata = ., expo = T)$pred)),
    lapply(res_indep_lin, function(m) 
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
betas_lin <- data_frame(
  b.gl = sapply(res_gl_lin, coef),
  b.indep = sapply(res_indep_lin, coef)
) %>%
  gather(cov, beta) %>%
  separate(cov, c("b", "cov")) %>%
  mutate(sim = rep(1:k_lin, 2))

ggplot(betas_lin, aes(x = beta, col = cov)) +
  geom_line(stat = "density")  +
  geom_vline(xintercept = beta_lin[2], linetype = "dotted") +
  scale_color_discrete(labels = c("Greenland-Longnecker", "Independence")) +
  labs(col = "Method", x = expression(hat(beta)))
# ggsave("figures/sim_lin_betadist.pdf")


# comparing standard error (or variances)
bvcov_lin <- data.frame(
  v.gl = sapply(res_gl_lin, vcov),
  v.indep = sapply(res_indep_lin, vcov)
) %>%
  gather(sim, vcov) %>%
  separate(sim, c("v", "cov")) %>% 
  mutate(sim = rep(1:k, 2))

ggplot(bvcov_lin, aes(x = cov, y = vcov)) +
  geom_boxplot() +
  labs(x = "Method") +
  theme_classic()

ggplot(bvcov_lin, aes(x = vcov, col = cov)) +
  geom_line(stat = "density") +
  labs(x = expression(paste("VAR(", hat(beta), ")")), col = "Method") +
  scale_color_discrete(labels = c("Greenland-Longnecker", "Independence"))
# ggsave("/figures/dist_se_method.pdf")

betas_lin %>% 
  group_by(cov) %>%
  summarise(sd(beta)) %>%
  bind_cols(
    bvcov_lin %>% group_by(cov) %>%
      summarise(mean(vcov)^.5)
  )

# p-values
mean(p_val_gl_lin < .05)
mean(p_val_indep_lin < .05)

# save(betas_lin, bvcov_lin, beta_lin, k_lin, p_val_gl_lin, p_val_indep_lin,
#      file = "data/res_cov_lin.RData")
# load("data/res_cov_lin.RData")

p_betas_lin <- ggplot(betas_lin, aes(x = beta, col = cov)) +
  geom_line(stat = "density")  +
  geom_vline(xintercept = beta_lin[2], linetype = "dotted") +
  scale_color_discrete(labels = c("Greenland-Longnecker", "Independence")) +
  labs(col = "Method", x = expression(hat(beta)))
p_se_lin <- ggplot(bvcov_lin, aes(x = vcov^.5, col = cov)) +
  geom_line(stat = "density") +
  labs(x = expression(paste("SE(", hat(beta), ")")), col = "Method") +
  scale_color_discrete(labels = c("Greenland-Longnecker", "Independence"))
p_sim_lin <- plot_grid( 
  p_betas_lin + theme(legend.position="none"),
  p_se_lin + theme(legend.position="none"),
  align = 'vh', labels = c("A", "B"), hjust = -1, nrow = 1)
legend_sim_lin <- get_legend(p_betas_lin + theme(legend.position="bottom"))
p_sim_lin <- plot_grid(p_sim_lin, legend_sim_lin, ncol = 1, rel_heights = c(1, .1))
p_sim_lin