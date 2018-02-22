## ----eval=FALSE----------------------------------------------------------
#  library(devtools)
#  install_github('sinhrks/ggfortify')
#  library(ggfortify)

## ----global_options, include=FALSE---------------------------------------
library(knitr)
opts_chunk$set(fig.width=6, fig.height=3, fig.path='figures/pca-', warning=FALSE)

## ---- message = FALSE, eval=FALSE----------------------------------------
#  library(ggfortify)
#  df <- iris[c(1, 2, 3, 4)]
#  autoplot(prcomp(df))

## ---- message = FALSE, eval=FALSE----------------------------------------
#  autoplot(prcomp(df), data = iris, colour = 'Species')

## ---- message = FALSE, eval=FALSE----------------------------------------
#  autoplot(prcomp(df), data = iris, colour = 'Species', label = TRUE, label.size = 3)

## ---- message = FALSE, eval=FALSE----------------------------------------
#  autoplot(prcomp(df), data = iris, colour = 'Species', shape = FALSE, label.size = 3)

## ---- message = FALSE, eval=FALSE----------------------------------------
#  autoplot(prcomp(df), data = iris, colour = 'Species', loadings = TRUE)

## ---- message = FALSE, eval=FALSE----------------------------------------
#  autoplot(prcomp(df), data = iris, colour = 'Species',
#           loadings = TRUE, loadings.colour = 'blue',
#           loadings.label = TRUE, loadings.label.size = 3)

## ---- message = FALSE, eval=FALSE----------------------------------------
#  d.factanal <- factanal(state.x77, factors = 3, scores = 'regression')
#  autoplot(d.factanal, data = state.x77, colour = 'Income')
#  autoplot(d.factanal, label = TRUE, label.size = 3,
#           loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 3)

## ---- message = FALSE, eval=FALSE----------------------------------------
#  
#  autoplot(kmeans(USArrests, 3), data = USArrests)
#  
#  autoplot(kmeans(USArrests, 3), data = USArrests, label = TRUE, label.size = 3)

## ---- message = FALSE, eval=FALSE----------------------------------------
#  library(cluster)
#  autoplot(clara(iris[-5], 3))

## ---- eval=FALSE---------------------------------------------------------
#  autoplot(fanny(iris[-5], 3), frame = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  autoplot(pam(iris[-5], 3), frame = TRUE, frame.type = 'norm')

## ---- message = FALSE, eval = FALSE--------------------------------------
#  library(lfda)
#  
#  # Fisher局部判别分析 (LFDA)
#  model <- lfda(iris[-5], iris[, 5], 4, metric="plain")
#  autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')

## ---- message = FALSE, eval = FALSE--------------------------------------
#  # 半监督Fisher局部判别分析 (SELF)
#  model <- self(iris[-5], iris[, 5], beta = 0.1, r = 3, metric="plain")
#  autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')

## ---- message = FALSE, eval=FALSE----------------------------------------
#  library(ggfortify)
#  autoplot(AirPassengers)

## ---- eval=FALSE---------------------------------------------------------
#  autoplot(AirPassengers, ts.colour = 'red', ts.linetype = 'dashed')

## ----message = FALSE, eval=FALSE-----------------------------------------
#  library(vars)
#  data(Canada)
#  autoplot(Canada)

## ---- eval=FALSE---------------------------------------------------------
#  autoplot(Canada, facets = FALSE)

## ---- message = FALSE, eval=FALSE----------------------------------------
#  library(xts)
#  autoplot(as.xts(AirPassengers), ts.colour = 'green')
#  
#  library(timeSeries, eval=FALSE)
#  autoplot(as.timeSeries(AirPassengers), ts.colour = ('dodgerblue3'))

## ---- message = FALSE, eval=FALSE----------------------------------------
#  autoplot(AirPassengers, ts.geom = 'bar', fill = 'blue')

## ---- message = FALSE, eval=FALSE----------------------------------------
#  autoplot(AirPassengers, ts.geom = 'point', shape = 3)

## ---- message = FALSE, eval=FALSE----------------------------------------
#  library(forecast)
#  d.arima <- auto.arima(AirPassengers)
#  d.forecast <- forecast(d.arima, level = c(95), h = 50)
#  autoplot(d.forecast)

## ---- eval=FALSE---------------------------------------------------------
#  autoplot(d.forecast, ts.colour = 'firebrick1', predict.colour = 'red',
#           predict.linetype = 'dashed', conf.int = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  library(vars)
#  d.vselect <- VARselect(Canada, lag.max = 5, type = 'const')$selection[1]
#  d.var <- VAR(Canada, p = d.vselect, type = 'const')

## ---- eval=FALSE---------------------------------------------------------
#  autoplot(predict(d.var, n.ahead = 50), ts.colour = 'dodgerblue4',
#           predict.colour = 'blue', predict.linetype = 'dashed')

## ----message = FALSE, eval=FALSE-----------------------------------------
#  library(changepoint)
#  autoplot(cpt.meanvar(AirPassengers))

## ---- eval=FALSE---------------------------------------------------------
#  autoplot(cpt.meanvar(AirPassengers), cpt.colour = 'blue', cpt.linetype = 'solid')

## ----message = FALSE, eval=FALSE-----------------------------------------
#  library(strucchange)
#  autoplot(breakpoints(Nile ~ 1), ts.colour = 'blue', ts.linetype = 'dashed',
#           cpt.colour = 'dodgerblue3', cpt.linetype = 'solid')

## ---- message=FALSE, warning=FALSE, eval=FALSE---------------------------
#  library(dlm)
#  form <- function(theta){
#    dlmModPoly(order = 1, dV = exp(theta[1]), dW = exp(theta[2]))
#  }
#  
#  model <- form(dlmMLE(Nile, parm = c(1, 1), form)$par)
#  filtered <- dlmFilter(Nile, model)
#  
#  autoplot(filtered)

## ---- eval=FALSE---------------------------------------------------------
#  autoplot(filtered, ts.linetype = 'dashed', fitted.colour = 'blue')

## ---- eval=FALSE---------------------------------------------------------
#  smoothed <- dlmSmooth(filtered)
#  class(smoothed)
#  autoplot(smoothed)

## ---- eval=FALSE---------------------------------------------------------
#  p <- autoplot(filtered)
#  autoplot(smoothed, ts.colour = 'blue', p = p)

## ---- message=FALSE, warning=FALSE, eval=FALSE---------------------------
#  library(KFAS)
#  model <- SSModel(
#    Nile ~ SSMtrend(degree=1, Q=matrix(NA)), H=matrix(NA)
#  )
#  
#  fit <- fitSSM(model=model, inits=c(log(var(Nile)),log(var(Nile))), method="BFGS")
#  smoothed <- KFS(fit$model)
#  autoplot(smoothed)

## ---- eval=FALSE---------------------------------------------------------
#  filtered <- KFS(fit$model, filtering="mean", smoothing='none')
#  autoplot(filtered)

## ---- eval=FALSE---------------------------------------------------------
#  trend <- signal(smoothed, states="trend")
#  class(trend)

## ---- eval=FALSE---------------------------------------------------------
#  p <- autoplot(filtered)
#  autoplot(trend, ts.colour = 'blue', p = p)

## ---- message = FALSE, eval=FALSE----------------------------------------
#  autoplot(stl(AirPassengers, s.window = 'periodic'), ts.colour = 'blue')

## ---- fig.width = 4, fig.height = 2, eval=FALSE--------------------------
#  autoplot(acf(AirPassengers, plot = FALSE))

## ---- fig.width = 4, fig.height = 2, eval=FALSE--------------------------
#  autoplot(acf(AirPassengers, plot = FALSE), conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma')

## ---- fig.width = 4, fig.height = 2, eval=FALSE--------------------------
#  autoplot(spec.ar(AirPassengers, plot = FALSE))

## ---- message = FALSE, fig.width=3, fig.height=3, eval=FALSE-------------
#  ggcpgram(arima.sim(list(ar = c(0.7, -0.5)), n = 50))

## ---- message = FALSE, fig.width = 5, fig.height = 5, eval=FALSE---------
#  library(forecast)
#  ggtsdiag(auto.arima(AirPassengers))

## ---- message = FALSE, fig.width = 4, fig.height = 4, eval=FALSE---------
#  gglagplot(AirPassengers, lags = 4)

