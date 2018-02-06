library(dosresmeta)
library(tidyverse)
library(rms)
library(scales)
library(cowplot)
library(knitr)
library(kableExtra)
library(Epi)
source("code/functions.R")

# reading data
data("coffee_mort")
data("coffee_mort_add")
coffee <- rbind(coffee_mort, coffee_mort_add) %>%
  arrange(id) %>%
  subset(area == "Europe")

Slist <- lapply(unique(coffee$id), function(i)
  with(subset(coffee, id == i), {
    if (any(is.na(cases) | is.na(n))){
      diag(se[se != 0 & !is.na(se)]^2, nrow = sum(se != 0 & !is.na(se)))
    }
    else {
      covar.logrr(y = logrr, v = I(se^2), cases = cases, n = n,
                  type = type)
    }
  }))
names(Slist) <- unique(coffee$id)
newd <- data.frame(dose = seq(0, 8, length.out = 100))
newd_tab <- data.frame(dose = c(0:6))


# With and without exclusion
twostage <- dosresmeta(logrr ~ dose + I(dose^2), id = id, se = se, 
                       data = subset(coffee, !(id %in% c(28, 29))), 
                       covariance = "user", Slist = Slist[!names(Slist) %in% c("28", "29")])
summary(twostage)
onestage <- dosresmeta(logrr ~ dose + I(dose^2), id = id, se = se, data = coffee, 
                  covariance = "user", Slist = Slist, proc = "1stage")
summary(onestage)
p_ts <- round(summary(twostage)$coefficients[2, 4], 3)
p_os <- round(summary(onestage)$coefficients[2, 4], 3)
se_ts <- round(summary(twostage)$coefficients[1:2, 2], 3)
se_os <- round(summary(onestage)$coefficients[1:2, 2], 3)
psi_ts <- round(c(twostage$Psi)[-2], 5)
psi_os <- round(c(onestage$Psi)[-2], 5)

# comparison of estimated quantities
cbind(twostage = coef(twostage), onestage = coef(onestage))
cbind(twostage = c(vcov(twostage))[-2], onestage = c(vcov(onestage))[-2]) %>% round(5)

# predicted curves
newd %>%
  mutate(
    twostage = predict(twostage, ., expo = T)$pred,
    onestage = predict(onestage, ., expo = T)$pred
  ) %>%
  gather(analysis, pred, -dose) %>%
  mutate(analysis = factor(analysis, labels = c("One-stage", "Two-stage"))) %>%
  ggplot(aes(dose, pred, linetype = analysis)) +
  geom_line() +
  scale_y_continuous(trans = "log", breaks = pretty_breaks()) + 
  labs(x = "Coffee consumption (cups/day)", y = "Relative risk", linetype = "Curve") +
  scale_linetype_manual(values = c(`One-stage` = "solid", `Two-stage` = "dashed"))

tab_coffee <- newd_tab %>%
  bind_cols(predict(twostage, ., expo = T)[, -c(1:2)]) %>%
  bind_cols(predict(onestage, ., expo = T)[, -c(1:2)]) %>%
  round(2)

# blup
bi_os <- data.frame(
  id = unique(coffee$id),
  xmin = tapply(coffee$dose, coffee$id, min),
  xref = tapply(coffee$dose, coffee$id, head, 1),
  xmax = tapply(coffee$dose, coffee$id, max),
  bi = t(apply(dosresmeta:::blup.dosresmeta(onestage), 1, function(x) x + rbind(coef(onestage))))
)
bi_ts <- data.frame(
  id = setdiff(unique(coffee$id), c(28, 29)),
  bi = t(apply(dosresmeta:::blup.dosresmeta(twostage), 1, function(x) x + rbind(coef(twostage))))
)
bi <- merge(bi_os, bi_ts, by = "id", all = T)
predi <- data.frame(
  x = newd$dose,
  y_os = pred_sq(coef =  bi[, c("bi.1.x", "bi.2.x")], newd$dose, bi$xref),
  y_ts = pred_sq(coef =  bi[, c("bi.1.y", "bi.2.y")], newd$dose, bi$xref)
)

p_coffee <- lapply(seq_along(unique(coffee$id)), function(i){
  idi <- unique(coffee$id)[i]
  ggplot(subset(coffee, id == idi), aes(x = dose, y = exp(logrr))) + 
    geom_errorbar(aes(ymin = exp(logrr - 1.96*se), ymax = exp(logrr + 1.96*se)), width = .3) +
    geom_line(data = predi, aes_string(x = "x", y = paste0("y_os.pred", i), lty = "'One-stage'"), lwd = 1) +
    geom_line(data = predi, aes_string(x = "x", y = paste0("y_ts.pred", i), lty = "'Two-stage'"), lwd = 1) +
    scale_y_continuous(trans = "log", breaks = pretty_breaks()) + 
    theme(legend.position = "none")  + xlim(c(0, 8)) +
    labs(title = paste0("Study ID ", idi), x = "Dose", y = "Relative Risk", lty = "Curvwe")
})
legend_grid_coffee <- get_legend(p_coffee[[1]] + theme(legend.position = "bottom"))
p_grid_coffee <- plot_grid(plotlist = c(p_coffee, list(legend_grid_coffee)), nrow = 5, ncol = 3,
                           rel_heights = c(1, 1, 1, 1, .1))


# model comparison
quadr_ml <- dosresmeta(logrr ~ dose + I(dose^2), id = id, se = se, data = coffee, 
                     covariance = "user", Slist = Slist, proc = "1stage", method = "ml")
spk_ml <- dosresmeta(logrr ~ I(1*(dose < 1)) + I((dose-1)*(dose >= 1)) +
                       I((dose-4)*(dose >= 4)), 
                        id = id, se = se, data = coffee, 
                        covariance = "user", proc = "1stage", method = "ml",
                        Slist = Slist, control = list(maxiter = 5000))
k2 <- c(0, 1, 3, 5, 7, 10)
categ <- dosresmeta(logrr ~ relevel(cut(dose, breaks = k2, include.lowest = T, right = F), 2), 
                    id = id, se = se, data = coffee, 
                    covariance = "user", proc = "1stage", method = "ml",
                    Slist = Slist, control = list(maxiter = 5000))
waldtest(vcov(categ), coef(categ), 2:4)
AIC_os <- round(sapply(list(quadr = quadr_ml, spk = spk_ml, categ = categ), AIC), 2)

xref <- 1
comp <- rbind(xref, newd) %>%
  mutate(
    quadr = predict(quadr_ml, ., expo = T)$pred,
    spike = predict(spk_ml, ., expo = T)$pred,
    categories = predict(categ, ., expo = T)$pred
  )
p_comp <- ggplot(comp, aes(dose, quadr, col = "Quadratic")) +
  geom_line() +
  geom_line(data = subset(comp, dose < 1), aes(y = spike, col = "Spike at 0")) +
  geom_line(data = subset(comp, dose >= 1), aes(y = spike, col = "Spike at 0")) +
  scale_y_continuous(trans = "log", breaks = pretty_breaks()) + 
  labs(x = "Coffee consumption (cups/day)", y = "Relative risk", col = "Model") +
  scale_color_manual(values = c(`Quadratic` = "blue", `Spike at 0` = "red",
                                `Categories` = "green"))
for (i in seq_along(k2[-1])){
  p_comp <- p_comp + geom_line(data = subset(comp, dose >= k2[i] & dose < k2[i+1]), 
                               aes(y = categories, col = "Categories"))
}
p_comp


# meta-regression
quadr_reg <- dosresmeta(logrr ~ dose + I(dose^2), id = id, se = se, data = coffee, 
                       covariance = "user", Slist = Slist, proc = "1stage",
                       mod = ~ gender, method = "reml")
summary(quadr_reg)
waldtest(vcov(quadr_reg), coef(quadr_reg), 3:6)

coffee_vpc <- coffee %>%
  filter(se != 0) %>%
  mutate(
    `Quadratic` = vpc(onestage),
    `Quadratic meta-regression` = vpc(quadr_reg)
  ) %>%
  gather(curve, vpc, Quadratic:`Quadratic meta-regression`)

p_coffee_vpc <- ggplot(coffee_vpc, aes(dose, vpc, group = curve)) +
  geom_point(aes(shape = curve, col = curve)) +
  geom_smooth(aes(col = curve), method = "loess", se = F) +
  labs(y = "VPC", x = "Coffee consumption (cups/day)", shape = "Curve", col = "Curve") +
  scale_color_manual(values = c(`Quadratic` = "blue", `Quadratic meta-regression` = "red")) +
  theme(legend.position = "top")
p_coffee_vpc
