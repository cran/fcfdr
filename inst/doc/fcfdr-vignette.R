## ---- include = FALSE,eval=FALSE----------------------------------------------
#  knitr::opts_chunk$set(
#    collapse = TRUE,
#    comment = "#>",
#    dev = 'png'
#  )
#  Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "false")

## ----eval=FALSE---------------------------------------------------------------
#  library(fcfdr)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  n = 50000
#  n1p = 500 # associated variants
#  zp = c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1)) # z-scores
#  p = 2*pnorm(-abs(zp)) # convert to p-values
#  hist(p)

## ----eval=FALSE---------------------------------------------------------------
#  mixture_comp1 <- function(x) rnorm(x, mean = -0.5, sd = 0.5)
#  mixture_comp2 <- function(x) rnorm(x, mean = 2, sd = 1)
#  n = length(p)
#  z = runif(n)
#  
#  q <- c(mixture_comp1(n1p), mixture_comp2(n-n1p))
#  hist(q)

## ---- fig.width = 6, fig.height = 5,eval=FALSE--------------------------------
#  corr_plot(p, q)

## ---- fig.width = 6, fig.height = 5,eval=FALSE--------------------------------
#  stratified_qqplot(data_frame = data.frame(p, q), prin_value_label = "p", cond_value_label = "q", thresholds = quantile(q)[-1])

## ----eval=FALSE---------------------------------------------------------------
#  res <- flexible_cfdr(p, q, indep_index = seq(1, n, 1))

## ----eval=FALSE---------------------------------------------------------------
#  str(res)
#  
#  p = res[[1]]$p
#  q = res[[1]]$q
#  v = res[[1]]$v

## ----eval=FALSE---------------------------------------------------------------
#  pv_plot(p = p, q = q, v = v)
#  log10pv_plot(p = p, q = q, v = v,
#               axis_lim = c(0, 10)) # zoom in to interesting region

## ----eval=FALSE---------------------------------------------------------------
#  hit = which(p.adjust(v, method = "BH") <= 0.05)

## ----eval=FALSE---------------------------------------------------------------
#  hit_p = which(p.adjust(p, method = "BH") <= 0.05)

## ----eval=FALSE---------------------------------------------------------------
#  # cFDR
#  1 - (length(intersect(hit,c(1:500)))/length(hit))
#  
#  # p-value
#  1 - (length(intersect(hit_p,c(1:500)))/length(hit_p))

## ----eval=FALSE---------------------------------------------------------------
#  # number of extra true associations identified by flexible cFDR
#  length(which(hit[!hit %in% hit_p] <= 500))

