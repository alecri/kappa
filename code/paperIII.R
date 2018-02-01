library(dosresmeta)
library(tidyverse)
library(metafor)
library(hetmeta)
library(knitr)
library(kableExtra)
library(rms)

source("code/functions.R")

# data(package = "dosresmeta")

# red meat
# processed meat

data("process_bc")
k_processed <- quantile(process_bc$dose, c(.10, .5, .90))
process <- process_bc %>%
  group_by(id) %>%
  nest() %>%
  mutate(
    Slist = map(data, ~ if (any(is.na(.x$cases) | is.na(.x$n))){
      diag(.x$se[.x$se != 0]^2, nrow = sum(.x$se != 0))
    } else {
      covar.logrr(cases = .x$cases, n = .x$n, y = .x$logrr, v = I(.x$se^2), 
                  type = .x$type)
    }),
    lin = map2(data, Slist, 
               ~ dosresmeta(logrr ~ dose, se = se, data = .x, covariance = "user",
                            Slist = .y)),
    quadr = map2(data, Slist, 
                 ~ dosresmeta(logrr ~ dose + I(dose^2), se = se, data = .x, covariance = "user",
                              Slist = .y)),
    spline = map2(data, Slist, 
                 ~ dosresmeta(logrr ~ rcs(dose, k_processed), se = se, data = .x, covariance = "user",
                              Slist = .y)),
    yi = 50*map_dbl(lin, coef),
    vi = 50^2*map_dbl(lin, vcov),
    author = map_chr(data, ~ as.character(.x$author[1])),
    year = map_chr(data, ~ as.character(.x$year[1])),
    type = map_chr(data, ~ as.character(.x$type[1]))
  )
p_pr_quadr <- summary(mvmeta(do.call("rbind", map(process$quadr, ~ coef(.x))) ~ 1, 
                             map(process$quadr, ~vcov(.x)), method = "mm"))$coefficients[2, 4]
p_pr_spline <- summary(mvmeta(do.call("rbind", map(process$spline, ~ coef(.x))) ~ 1, 
                              map(process$spline, ~vcov(.x)), method = "mm"))$coefficients[2, 4]
lin_pr <- rma.uni(yi = yi, vi = vi, data = process, method = "DL")
pred_pr <- pred.rma(lin_pr, transf = exp, digits = 2)
hetmes_pr <- confint(lin_het_pr)

# forest plot
xlim <- c(-1.25, 2.25)
old <- par()
par(bg = "white", cex = 1, font = 1, las = 1)
forest(lin_pr, slab = paste0(process$author, ", ", process$year), 
       atransf = exp, order = order(process$type), showweights = T,
       ylim = c(-1, 14), xlab = "",
       xlim = xlim, at=log(c(.6, 1, 1.5, 2, 3))
)
par(font = 2)
text(xlim[1], 12.5, "Author(s), Year", pos = 4)
text(xlim[2]-.05, 12.5, "RR [95% CI]", pos = 2)
text(xlim[2]-.55, 12.5, "Weight", pos = 2)
mtext("processed meat and bladder cancer", side = 3, line = 1, cex = 1.5)
mtext("for every 50 g per day increment", side = 3, line = -1, cex = 1.5)

ggplot(NULL, aes(lin_pr$vi)) + 
  geom_line(stat = "density", adjust = 2) +
  labs(x = "Within-study error variance")


tab_rb <- data_frame(
  exposure = "Processed meat",
  beta =  pred.rma(lin_pr, transf = exp, string = T),
  Qtest = paste(round(c(lin_pr$QE, lin_pr$QEp), c(1, 3)), collapse = ", "),
  CV_vi = lin_het_pr$cv_vi
) %>%
  cbind(
    rbind(apply(round(confint(lin_het_pr), c(0, 0, 0, 2)), 1, function(x) 
      paste0(x[1], " (", x[3], ", ", x[4], ")")))
  )



## red meat

data("red_bc")
k_red <- quantile(red_bc$dose, c(.10, .5, .90))
red <- red_bc %>%
  group_by(id) %>%
  nest() %>%
  mutate(
    Slist = map(data, ~ if (any(is.na(.x$cases) | is.na(.x$n))){
      diag(.x$se[.x$se != 0]^2, nrow = sum(.x$se != 0))
    } else {
      covar.logrr(cases = .x$cases, n = .x$n, y = .x$logrr, v = I(.x$se^2), 
                  type = .x$type)
    }),
    lin = map2(data, Slist, 
               ~ dosresmeta(logrr ~ dose, se = se, data = .x, covariance = "user",
                            Slist = .y)),
    quadr = map2(data, Slist, 
                 ~ dosresmeta(logrr ~ dose + I(dose^2), se = se, data = .x, covariance = "user",
                              Slist = .y)),
    spline = map2(data, Slist, 
                  ~ dosresmeta(logrr ~ rcs(dose, k_red), se = se, data = .x, covariance = "user",
                               Slist = .y)),
    yi = 100*map_dbl(lin, coef),
    vi = 100^2*map_dbl(lin, vcov),
    author = map_chr(data, ~ as.character(.x$author[1])),
    year = map_chr(data, ~ as.character(.x$year[1])),
    type = map_chr(data, ~ as.character(.x$type[1])),
    area = map_chr(data, ~ .x$area[1]),
    unit = map_chr(data, ~ as.character(.x$unit[1]))
  )

lin_red <- rma.uni(yi = yi, vi = vi, data = red, method = "DL")
lin_redcc <- rma.uni(yi = yi, vi = vi, data = subset(red, type == "cc"), method = "DL")
lin_redci <- rma.uni(yi = yi, vi = vi, data = subset(red, type != "cc"), method = "DL")
lin_het_red <- hetmeta(lin_red)
lin_het_redcc <- hetmeta(lin_redcc)
lin_het_redci <- hetmeta(lin_redci)
hetmes_red <- confint(lin_het_pr)

tab_rb <- rbind(tab_rb,
                data_frame(
                  exposure = "Red meat",
                  beta =  pred.rma(lin_red, transf = exp, string = T),
                  Qtest = paste(c(round(lin_red$QE, 1), "< 0.01"), collapse = ", "),
                  CV_vi = lin_het_red$cv_vi
                ) %>%
                  cbind(
                    rbind(apply(round(confint(lin_het_red), c(0, 0, 0, 2)), 1, function(x) 
                      paste0(x[1], " (", x[3], ", ", x[4], ")")))
                  ))


ggplot(NULL, aes(lin_pr$vi, linetype = "Processed meat")) +
  geom_line(stat = "density", adjust = 2.5) +
  geom_line(aes(lin_red$vi, linetype = "Read meat"), 
            stat = "density", adjust = 2.5) +
  labs(x = "Within-study error variance", color = "Exposure")


lin_red_reg <- rma.uni(yi = yi, vi = vi, data = red, method = "DL",
                       mods = ~ I(type == "cc"))
hetmeta(lin_red_reg)

forest(lin_pr)
forest(lin_red)
lin_red_cc <- rma.uni(yi = yi, vi = vi, data = filter(red, type == "cc"), method = "DL")
hetmeta(lin_red_cc)
lin_red_ci <- rma.uni(yi = yi, vi = vi, data = filter(red, type != "cc"), method = "DL")
hetmeta(lin_red_ci)
