library(dosresmeta)
library(tidyverse)
library(metafor)
library(hetmeta)
library(rms)

source("code/functions.R")

# processed meat: linear, quadratic and spline analysis
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
lin_pr <- rma.uni(yi = yi, vi = vi, data = process, method = "DL")
pred_pr <- pred.rma(lin_pr, transf = exp, digits = 2)

# heterogeneity measures with 95% confidence intervals
lin_het_pr <- hetmeta(lin_pr)
confint(lin_het_pr)

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

# disitribution for v_i
ggplot(NULL, aes(lin_pr$vi)) + 
  geom_line(stat = "density", adjust = 2) +
  labs(x = "Within-study error variance")

# saving results in a table
tab_rb <- data_frame(
  analysis = "Processed meat",
  beta =  pred.rma(lin_pr, transf = exp, string = T),
  Qtest = paste(round(c(lin_pr$QE, lin_pr$QEp), c(0, 2)), collapse = ", "),
  CV_vi = round(lin_het_pr$cv_vi, 2)
) %>%
  cbind(
    rbind(apply(round(confint(lin_het_pr)[-4, ]), 1, function(x) 
      paste0(x[1], " (", x[3], ", ", x[4], ")")))
  )


## red meat, similar analysis but seperate for study desing
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
confint(lin_het_pr)
confint(lin_het_redcc)
confint(lin_het_redci)

# forest plot
a <- table(red$type)["cc"]
b <- sum(table(red$type)[c("ir", "ci")])
xlim <- c(-1.75, 2.75)
old <- par()
par(bg = "white", cex = 1, font = 1, las = 1)
forest(lin_red, slab = paste0(red$author, ", ", red$year), 
       atransf = exp, order = order(red$type), showweights = T,
       rows = c(3:(a+2), (6+a):(5+a+b)), ylim = c(-3, a+b+9), xlab = "",
       mlab = Qlab(lin_red), xlim = xlim, at=log(c(.65, 1, 1.5, 2, 3.5))
)
text(xlim[1], c(6 + a + b, 3 + a), c("Cohort", "Case-control"), pos = 4)
par(font = 2)
text(xlim[1], a + b + 7.5, "Author(s), Year", pos = 4)
text(xlim[2]-.05, a + b + 7.5, "RR [95% CI]", pos = 2)
text(xlim[2]-.9, a + b + 7.5, "Weight", pos = 2)
par(font = 1)
addpoly(lin_redcc, row = 2, cex = 1, atransf = exp, mlab = Qlab(lin_redcc, F))
addpoly(lin_redci, row = 5+a, cex = 1, atransf = exp, mlab = Qlab(lin_redci, F))
par(font = 2)
mtext("Red meat and bladder cancer", side = 3, line = 2, cex = 1.25)
mtext("for every 100 g per day increment", side = 3, line = 1, cex = 1.25)


# adding results to the previous table
tab_rb <- rbind(tab_rb,
                data.frame(analysis = c("Red meat", "Red meat, Prospective", "Red meat, Case-control")) %>%
                  cbind(
                    do.call("rbind", 
                            map(list(lin_red, lin_redci, lin_redcc), ~
                                  cbind(pred.rma(.x, transf = exp, string = T),
                                        paste(c(round(.x$QE), ifelse(.x$QEp < 0.01, "< 0.01", round(.x$QEp, 1))),
                                              collapse = ", "),
                                        round(hetmeta(.x)$cv_vi, 2)
                                  )
                            ))  
                  ) %>%
                  cbind(
                    do.call("rbind", map(list(confint(lin_het_red), confint(lin_het_redci), confint(lin_het_redcc)), 
                                         ~ apply(round(.x[-4, ]), 1, function(x) 
                                           paste0(x[1], " (", x[3], ", ", x[4], ")"))
                    ))
                  ) %>%
                  `colnames<-`(names(tab_rb))
)
tab_rb

# disitribution for v_i
ggplot(NULL, aes(lin_pr$vi, linetype = "Processed meat")) +
  geom_line(stat = "density", adjust = 2.5) +
  geom_line(aes(lin_red$vi, linetype = "Read meat"), 
            stat = "density", adjust = 2.5) +
  labs(x = "Within-study error variance", color = "Exposure")
