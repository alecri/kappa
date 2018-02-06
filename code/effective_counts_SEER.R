## SEER data: effective counts

library(dosresmeta)
library(tidyverse)
library(Epi)
library(scales)

# data selected from paper IV
load("data/breast_1501.RData")

# linear model on IPD
lin_IPD <- glm(Y ~ age + histology + maritalStatus + race +
                 tumor_size + radiation + grade + ER + PR + node_pos, 
               data = breast_1501, family = "poisson", offset = log(surv_time))
summary(lin_IPD)

# categorical model on IPD
k_cat <- quantile(breast_1501$age, c(0, .25, .5, .75, 1))
breast_1501$age_cat <- cut(breast_1501$age, k_cat, include.lowest = T)

pos <- 1
cat_1501 <- glm(Y ~ relevel(age_cat, 1) + histology + maritalStatus + race +
                  tumor_size + radiation + grade + ER + PR + node_pos, 
                data = breast_1501, family = "poisson", offset = log(surv_time))
summary(cat_1501)

# creatig aggregated data
cat_res <- ci.exp(cat_1501)[1:4,] %>% 
  data.frame() %>%
  mutate(
    ci = paste0("(", round(X2.5., 2), ", ", round(X97.5., 2), ")"),
    logrr = coef(cat_1501)[1:4],
    se = diag(vcov(cat_1501))[1:4]^.5,
    refpos = c(pos, setdiff(1:4, pos))
  ) %>%
  arrange(refpos)
cat_res[pos, ] <- c(1, 1, 1, NA, 0, 0, pos)
tab_1501 <-breast_1501 %>%
  group_by(age_cat) %>%
  summarise(dose = mean(age),
            cases = sum(Y),
            PT = sum(surv_time)) %>%
  cbind(cat_res)


# linear model on AD using method by GL
lin_ad_gl <- dosresmeta(logrr ~ dose, type = "ir", se = se, cases = cases, n = PT,
                     data = tab_1501)
# linear model on AD using method by Hamling
lin_ad_h <- dosresmeta(logrr ~ dose, type = "ir", se = se, cases = cases, n = PT,
                        data = tab_1501, covariance = "h")

ec_h <- hamling(y = logrr, v = se^2, cases = cases, n = PT, type = "ir", data = tab_1501) %>%
  as.data.frame() %>%
  mutate(
    cases = round(A),
    n = round(N),
    controls = round(N - A),
    dose = tab_1501$dose
    )

lin_glm_cc <- glm(cbind(cases, controls) ~ dose, data = ec_h, family = "binomial")
lin_glm_ci <- glm(cases ~ dose, data = ec_h, family = "poisson", offset = log(n))

modi <- list(lin_ad_gl, lin_ad_h, lin_glm_cc, lin_glm_ci)
results <- data_frame(
  analysis = c("GL", "Hamling", "logit EC", "poisson EC"),
  coef = c(sapply(modi[1:2], coef), sapply(modi[3:4], function(x) coef(x)[2])),
  se = c(sapply(modi[1:2], vcov), sapply(modi[3:4], function(x) vcov(x)[4]))^.5,
  ci.lb = coef - 1.96*se, ci.ub = coef + 1.96*se
)
results

ggplot(results, aes(x = analysis, y = exp(coef))) +
  geom_point() +
  geom_errorbar(aes(ymin = exp(ci.lb), ymax = exp(ci.ub)), width = 0.2) +
  coord_flip() +
  scale_y_continuous(trans = "log", breaks = pretty_breaks(), limits = c(.98, 1.02)) +
  labs(y = "Relative risk (95% CI)", x = "")
