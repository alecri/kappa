# Graphs to be included as a collage in my thesis' cover

library(metafor)
library(tidyverse)
library(dosresmeta)
library(rms)
library(hetmeta)

options(show.signif.stars = FALSE)
dat <- get(data(coffee_mort))

k <- quantile(dat$dose, c(.1, .5, .9))
spl <- dosresmeta(logrr ~ rcs(dose, k), id = id, se = se, type = type,
                  cases = cases, n = n, data = dat, proc = "1stage")

# plot 1 & 4: code output from dosresmeta, two parts:
#  1) summary from a one-stage model
#  2) results for goodness-of-fit

pdf <- pipe("a2ps -o - | ps2pdf - code.pdf", "w")
capture.output(
  summary(spl), cat("\n\n\n"), gof(spl), file = pdf)
close(pdf)


xref <- 0
newd <- data.frame(dose = c(xref, seq(0, 8, length.out = 100)))
pred_spl <- function(x, k, b, mod){
  X <- rcs(x, k, inclx = T)
  X <- (apply(X, 2, function(x) x - x[1]))
  exp(X %*% (b + coef(mod)))
}

predi <- newd %>% cbind(
  apply(blup(spl), 1, function(b) pred_spl(newd$dose, k, b, spl))
  )

# plot 2: individual curves + combined dose-response
predi %>%
  mutate(pooled = predict(spl, newd, expo = T)$pred) %>%
  gather(study, pred, -dose) %>%
  ggplot(aes(dose, pred, group = study)) +
  geom_line(aes(size = study == "pooled", col = study == "pooled")) +
  scale_y_continuous(trans = "log", breaks = c(.8, .85, .9, .95, 1)) +
  labs(x = "Exposure", y = "Relative Risk") + guides(col = FALSE, size = FALSE) + 
  scale_size_manual(values = c(.35, 1.5)) +
  scale_color_manual(values = c("grey", "black")) +
  theme_classic()
ggsave("predi.pdf", width = 5, height = 3)


lin <- dosresmeta(logrr ~ dose, id = id, se = se, type = type,
                  cases = cases, n = n, data = dat)
lin_rma <- rma.uni(yi = lin$bi, vi = unlist(lin$Si), method = "REML")

# Plot 3: forest plot of linear trend
pdf("forest.pdf")
par(font = 1)
forest(lin_rma, atransf = exp, at = log(c(.6, .7, .8, .9, 1, 1.2)),
       mlab = paste0("Overall (Rb = ", round(hetmeta(lin_rma)$Rb), "%, p < 0.01)"), 
       alim = c(-.7, .3), xlim = c(-.7, .3),
       xlab = "", slab = NA, annotate = F)
par(font = 2)
text(-.35, 24, "Linear trend", pos = 4)
dev.off()


# additional plot: combined curve
data.frame(dose = c(2, seq(0, 8, length.out = 100))) %>%
  cbind(
    predict(spl, ., expo = T)
  ) %>%
  ggplot(aes(dose, pred)) +
  geom_line() +
  geom_segment(x = 2, y = -1, xend = 2, yend = 0, linetype = "dashed") +
  geom_segment(x = 2, y = 0, xend = 0, yend = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = .2) +
  scale_y_continuous(trans = "log", breaks = c(.95, 1, 1.05, 1.1, 1.15, 1.2)) +
  labs(x = "Exposure", y = "Relative Risk") +
  theme_classic()
