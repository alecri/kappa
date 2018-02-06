library(dosresmeta)
library(tidyverse)
library(gridExtra)
library(rms)
library(metafor)
library(hetmeta)
library(cowplot)
library(knitr)
library(kableExtra)

## example with processed meat

data("red_bc")

chr_red <- red_bc %>%
  group_by(id) %>%
  summarise(xref = dose[se == 0], min = min(dose), p25 = quantile(dose, .25),
            median = median(dose), p75 = quantile(dose, .75), max = max(dose))
kable(chr_red, "latex", align = "c", booktabs = T, escape = F, digits = 1, 
      col.names = c("ID", "Referent", "Min", "P25", "Median", "P75", "Max"), 
      caption = "Descriptive statistics of the assigned dose levels for 13 studies 
      included in a dose--response meta-analysis between red meat consumption and bladder
      cancer risk")

ggplot(chr_red, aes(x = min, y = seq_along(id))) +
  geom_point(aes(x = xref), shape = 4, size = 1.5) +
  geom_segment(aes(xend = max, yend = seq_along(id)), linetype = "dotted") +
  geom_point(data = filter(red_bc, se != 0) %>%
               mutate(id = replace(id, id > 12, id[id > 12] - 1)), aes(x = dose, y = id)) +
  labs(x = "Red meat consumption (g per day)", y = "Study ID") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank())

fpgrid <- as.matrix(fpgrid())
#max(chr_red$ref)
shift <- 10
scale <- 10

Slist_red <- do.call("list", by(red_bc, red_bc$id, function(x){
  if (any(is.na(x$n))){
    diag(x$se[x$se!=0]^2, ncol = length(x$se[x$se!=0]))
  } else {
    covar.logrr(y = logrr, v = I(se^2), cases = cases, n = n,
                type = type, data = x)
  }
}))

# individual fp2 model: 13 lists of 36
mod_pi_red <- map2(split(red_bc, red_bc$id), Slist_red,
                  ~ lapply(array_branch(fpgrid, 1), function(p)
                    dosresmeta(logrr ~ fracpol(dose, p = p, shift = shift, scale = scale), 
                               se = se, data = .x, covariance = "user", Slist = .y)
                  ))
AIC_pi_red <- lapply(mod_pi_red, function(x) sapply(x, AIC))
pb_red <- lapply(AIC_pi_red, which.min)
fpgrid[unlist(pb_red), ]
# 13 lists of 1 model 
mod_pb_red <- Map(function(x, y) x[[y]], mod_pi_red, pb_red)

# combined fp2 model: list of 36
mod_pi1_red <- map(array_branch(fpgrid, 1),
                  ~ dosresmeta(logrr ~ fracpol(dose, p = .x, shift = shift, scale = scale),
                               id = id, se = se, data = red_bc, covariance = "user", 
                               Slist = Slist_red, proc = "1stage"))
# final model 
mod_pb1_red <- mod_pi1_red[[which.min(sapply(mod_pi1_red, AIC))]]

pbi_plot_red <- Map(function(d, m, m1){
  newd <- data.frame(dose = seq(min(d$dose), max(red_bc$dose), length.out = 100)) %>%
    cbind(predict(m, newdata = ., expo = T)) %>%
    mutate(pred1 = predict(m1, newdata = ., expo = T)$pred)
  ggplot(subset(newd, dose <= max(d$dose)), aes(dose, pred)) +
    geom_line() +
    geom_line(aes(y = pred1), col = "red") +
    geom_line(data = subset(newd, dose > max(d$dose)), linetype = "dashed") +
    geom_line(data = subset(newd, dose >= max(d$dose)-1), aes(y = pred1), col = "red", linetype = "dashed") +
    geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = .15) +
    geom_ribbon(data = subset(newd, dose >= max(d$dose)-.5), aes(ymin = ci.lb, ymax = ci.ub), alpha = .05) +
    #geom_errorbar(data = d, aes(x = dose, y = rr, ymin = lrr, ymax = urr)) +
    scale_y_continuous(trans = "log", limits = c(max(.2, min(newd$ci.lb)), min(25, max(newd$ci.ub))), 
                       breaks = scales::pretty_breaks()) +
    labs(x = "Red meat (g per day)", y = "Relative risk", 
         title = paste0("Study ID ", d$id[1]))
}, split(red_bc, red_bc$id), mod_pb_red, Map(function(x) x[[which.min(sapply(mod_pi1_red, AIC))]], mod_pi_red))
layout <- matrix(c(1:13, NA, NA), nrow = 5, ncol = 3, byrow = TRUE)
grid.arrange(grobs = pbi_plot_red, nrow = 5, ncol = 3, as.table = FALSE, 
             layout_matrix = layout)

# pointwise strategy
xref <- 85
newd <- data.frame(
  dose = c(xref, seq(min(red_bc$dose), max(red_bc$dose), 5))
  ) %>%
  filter(dose >= 5, dose <= 300)
newd <- newd[-which(newd$dose == xref)[2], , drop = FALSE]

predi <- cbind(x = newd[-1, ], do.call("cbind", Map(function(m){
  predict(m, newdata = newd, exp = F, se = T)[-1, c("pred", "se")]
}, mod_pb_red)))
metamodi <- lapply(split(predi, seq_len(nrow(predi))), function(m)
  rma.uni(y = unlist(m)[grep("pred", names(m))], 
          sei = unlist(m)[grep("se", names(m))], method = "DL")
)
predavg <- cbind(x = newd[-1, ], as.data.frame(do.call("rbind", lapply(metamodi, function(x)
  predict(x, transf = exp)))))

extr <- do.call("cbind", lapply(split(red_bc, red_bc$id), function(d){
  newd$dose[-1] >= min(d$dose) & newd$dose[-1] <= max(d$dose)
}))
predi_extr <- predi
predi_extr[, grep("pred", names(predi))][!extr] <- NA
predi_extr[, grep("se", names(predi))][!extr] <- NA
metamodi_extr <- suppressWarnings(
  lapply(split(predi_extr, seq_len(nrow(predi_extr))), function(m)
    rma.uni(y = as.double(m[grep("pred", names(m))]), 
            sei = as.double(m[grep("se", names(m))]), method = "DL")
  ))
predavg_extr <- rbind(c(xref, 1, 1, 1),
                      cbind(x = newd[-1, ], do.call("rbind", lapply(metamodi_extr, function(x)
                        do.call("data.frame", predict(x, transf = exp))[c("pred", "ci.lb", "ci.ub")])))) %>%
  arrange(x)

# matplot(newd$dose[-1], exp(predi_extr[grep("pred", names(predi_extr))]),
#         log = "y", type = "l", col = "black", lty = 1, lwd = 3,
#         las = 1, bty = "n", xlab = "Red meat consumption (g per day)", ylab = "Risk ratio")
# matlines(newd$dose[-1], exp(predi[grep("pred", names(predi))]),
#          col = "black", lty = 2, lwd = 1.5)

n <- sapply(metamodi_extr, function(m) m$k)
i <- which(predavg_extr$x %in% seq(0, 300, 25))
predfp1 <- cbind(dose = newd$dose, predict(mod_pb1_red, newd, expo = T))[-1, ]
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for the new axis
with(predavg_extr, errbar(x[i], pred[i], ci.lb[i], ci.ub[i],
                          log = "y", las = 1, bty = "n", col = c("black"), 
                          lwd = c(1.25), pch = 19, cap = 0.01, xlim = c(0, 300),
                          ylim = c(.5, 3),
                          xlab = "Red meat consumption (g per day)", ylab = "Risk ratio"))
abline(h = 1)
with(subset(predfp1, dose %in% seq(0, 350, 25)),
     errbar(dose + 5, pred, ci.lb, ci.ub, col = c("black"), add = T,
            lty = 2, lwd = c(1.25), cap = 0.01))
par(new = TRUE)
plot(newd$dose[-1], n, axes = F, ylim = c(3, 120), xlab = "", cex = .2)
axis(side = 4, at = c(3, 7, 12), las = 1)
mtext("Number of studies                                    ", 
      side = 4, line = 3)
legend("topleft", exp(4), c("Point-wise average", "Two-stage"), lty = c("solid", "dashed"), 
       bty = "n", cex = 1.1, lwd = 2)


# # alternatevely
# with(predavg_extr, plot(c(x[i]), pred[i], type = "b",
#                           log = "y", las = 1, bty = "n", col = c("black"),
#                           xlab = "Red meat consumption (g per day)", ylab = "Risk ratio"))
# with(subset(predfp1, dose %in% seq(0, 350, 25)) %>% arrange(dose),
#      points(dose, pred, col = c("black"), type = "b", pch = 2))

# pointwise results
metamodi_res <- newd %>%
  filter(dose != xref) %>%
  mutate(
    tau2 = map_dbl(metamodi_extr, ~ .x$tau2),
    QEp = map_dbl(metamodi_extr, ~ .x$QEp),
    Rb = map_dbl(metamodi_extr, ~ hetmeta(.x)$Rb),
    I2 =  map_dbl(metamodi_extr, ~ .x$I2)
    # ,
    # vi = map(metamodi_extr, ~ .x$vi),
    # wi = map(metamodi_extr, ~ ((.x$vi + .x$tau2)^-1)/sum((.x$vi + .x$tau2)^-1))
  )
id_i <- c(1, 5, 10)
vi <- (predi_extr[grep("se", names(predi_extr))]^2)
colnames(vi) <- sub("se", "vi", colnames(vi))
vi_id <- vi[, paste0(id_i, ".vi")]
wi <- do.call("rbind", map2(array_branch(vi, 1), metamodi_res$tau2,
     ~ ((.x + .y)^-1)/sum((.x + .y)^-1, na.rm = T)))
colnames(wi) <- sub("vi", "wi", colnames(vi))
wi_id <- wi[, paste0(id_i, ".wi")]
metamodi_res <- cbind(metamodi_res, vi_id, wi_id)

xlim <- 250
p_pwa1 <- subset(metamodi_res, dose < xlim) %>%
  ggplot(aes(dose, tau2)) +
  geom_point() +
  geom_line() +
  labs(x = "Red meat consumption (g per day)", y = expression(tau^2))
p_pwa2 <- metamodi_res %>%
  subset(dose < 125) %>%
  gather(study, wi, ends_with("wi")) %>%
  mutate(study = factor(study, labels = id_i)) %>%
  filter(dose < xlim) %>%
  ggplot(aes(dose, wi, group = study)) +
  geom_point(aes(shape = study)) +
  geom_line(aes(linetype = study)) +
  labs(x = "Red meat consumption (g per day)", 
       y = expression(w[i]), shape = "study ID", linetype = "study ID") +
  theme(legend.position = c(.8, 0.6))
p_pwa3 <- subset(metamodi_res, dose < xlim & dose > 5) %>%
  ggplot(aes(dose, QEp)) +
  geom_point() +
  geom_line() +
  labs(x = "Red meat consumption (g per day)", 
       y = "p-value for Q test")
p_pwa4 <- subset(metamodi_res, dose < xlim) %>%
  ggplot(aes(dose, Rb)) +
  geom_point() +
  geom_line() +
  labs(x = "Red meat consumption (g per day)", 
       y = expression(R[b]))

#grid.arrange(p_pwa1, p_pwa2, p_pwa3, p_pwa4, ncol = 2, nrow = 2)
plot_grid(p_pwa1, p_pwa2, p_pwa3, p_pwa4, align = 'vh', labels = c("A", "B", "C", "D"))


# table for appendix
lapply(AIC_pi_red, round, 2) %>%
  do.call("cbind", .) %>%
  `rownames<-`(apply(fpgrid, 1, function(f) paste0("(", paste(f, collapse = ", "), ")"))) %>%
  rbind(
    `Best p` = apply(fpgrid[unlist(pb_red), ], 1, function(f) paste(f, collapse = ", "))
  ) %>%
  as.data.frame() %>%
  rownames_to_column(var = "p") %>%
  kable("latex", align = "c", booktabs = T, escape = F, longtable = T,
        caption = "AIC for the study-specific second-degree fractional polynomials with power terms specified by $p$ in a dose-response meta-analysis between red mead consumption (g per day) and bladder cancer risk. The last row reports the power term corresponding to the lowest AIC.") %>%
  add_header_above(c(" " = 1, "Study ID" = 13)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  landscape()

# old backup
lapply(AIC_pi_red, round, 2) %>%
  do.call("cbind", .) %>%
  `rownames<-`(apply(fpgrid, 1, function(f) paste0("(", paste(f, collapse = ", "), ")"))) %>%
  rbind(
    `best p` = apply(fpgrid[unlist(pb_red), ], 1, function(f) paste(f, collapse = ", "))
  ) %>%
  as.data.frame() %>%
  rownames_to_column(var = "p") %>%
  kable("latex", align = "c", booktabs = T, escape = F, 
        caption = "AIC for the study-specific second-degree fractional polynomials with power terms specified by $p$ in a dose-response meta-analysis between red mead consumption (g per day) and bladder cancer risk. The last row reports the power term corresponding to the lowest AIC.") %>%
  kable_styling(font_size = 8)

