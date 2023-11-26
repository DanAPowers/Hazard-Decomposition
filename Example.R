# Example.R 
setwd("~....")
library(tidyverse)
library(sandwich)

df <- read.csv(file='hazdata.csv')

##########################################################################
# Notes:
#  data should be constructed as person-epsiod/counting-process format
#  model formulas should include offset()
#   -set to 0 or log interval width in discrete-time models
#    (if NULL will be set to 0)
#   -set to log exposure: log(t[i] - t[i-1]) in continuous-time models
#   -account for frailty
#             offset(logexp + logv)
#   where v is a known relative risk at the 
#   person level or higher (i.e., persons nested in clusters)
#   (example uses sibling data identified by famid) 
#   details in Powers and Yun (2008) Soc. Methodology
#   decompose differences in hazard of premarital birth betweeen 
#   Black and White women in the NLSY79 data.
############################################################################
#                     D A T A   D E S C R I P T I O N
# data are structured as split-episodes based on age at risk of a first 
# premarital birth: [12,16) [16,18) [18,20) [20,22) [22,24) [24,36) 
# (mean-centered except for magebir) : 
#  pctsmom  = proportion of years from age 0-18 spent in single mother household
#  nfamtran = number of family-structure changes from 0-18
#  medu     = respondent mother's years of schooling)
#  inc1000  = adjusted family income in thousands of $
#  nosibs   = number of older siblings
#  magebir  = respondent mother's age at respondent's birth)
#  age      = factor for age-intervals 
#  race     = factor for race
#  famid    = family identifier for sibling clusters
#  iid      = individual indentifier
#############################################################################
# USER INPUT: model formula
# continuous time formula (supplied by user)
uformC <- devnt ~ age +
  pctsmom + nfamtran + medu + inc1000 + nosibs + magebir + offset(logexp) - 1 

# discrete time formula (supplied by user)
uformD <- devnt ~ a1 + a2 + a3 + a4 + a5 + a6 +
  pctsmom + nfamtran + medu + inc1000 + nosibs + magebir + offset(fake) - 1 

############################################################################
# function to assemble model results for input to decomp functions
# (update as needed for additional models)
############################################################################

mod.out <- function(m, gval, df, group, cluster=NULL,...) {
  b   <- coef(m)
  x   <- m$x
  # robust types also possible
  # 
  if(is.null(cluster))  v <- summary(m)$cov.unscaled
  if(!is.null(cluster)) v <- vcovCL(m, cluster = df$"cluster"[df$"group"== gval]) 
  w      <- m$prior.weights
  off    <- m$offset
  if (is.null(off)) off <- rep(0,nrow(x))
  id     <- with(df, iid[`group` == gval])
  mget(ls())
}  

# USER INPUT 
######################################################
# define binary grouping variable   
group  <- with(df, ifelse(race=="Black",1,0)) 
gv     <- group==1 # logical for glm subset arg
gval   <- c(0,1)   # index for post model processing
# set cluster variable if needed
cluster <- with(df, famid)
#  weights, etc. for glm can be included as needed
#  df$w is a random person weight
######################################################

# Fit Models by Subset

# Blacks (gv = TRUE = blk=1)
## continuous time (piecewise exponential)
mod.1a <- glm(terms(uformC), data = df, sub = gv, family = poisson, x = TRUE)
## discrete time (logit)
#mod.1a <- glm(terms(uformD), data = df, sub = gv, family = binomial(link=logit), x = TRUE)
## discrete time (cloglog)
# mod.1a <- glm(terms(uformD), data=df, sub=group==1, family = binomial(link=cloglog), x = TRUE)

mean(fitted(mod.1a))

# Whites (gv = FALSE = blk=0)
## continuous time (piecewise exponential)
mod.0a <- glm(terms(uformC), data=df, sub = !gv, family='poisson', x=TRUE)
## discrete time (logit)
#mod.0a <- glm(terms(uformD), data=df, sub = !gv,  family=binomial(link=logit), x=TRUE)
# mean(fitted(mod.1a))
## discrete time (cloglog)
# mod.0a <- glm(terms(uformD), data=df, sub =  !gv, family=binomial(link=cloglog), x=TRUE)

mean(fitted(mod.0a))

####################################
# gval: identifies blacks=1 whites=0 
####################################
outm0 <- mod.out(mod.0a, gval=0 , df, group, cluster=cluster)  # low outcome model
# check that n's are equal for all args
cbind(length(outm0$id), nrow(outm0$x), length(outm0$w))
outm1 <- mod.out(mod.1a, gval=1,  df, group, cluster=cluster)  # high outcome model
# check that n's are equal for all args
cbind(length(outm1$id), nrow(outm1$x), length(outm1$w))
## to pass to decomp subroutines
m0stuff <- list(b=outm0$b,x=outm0$x, off=outm0$off, v=outm0$v, w=outm0$w, id=outm0$id)
m1stuff <- list(b=outm1$b,x=outm1$x, off=outm1$off, v=outm1$v, w=outm1$w, id=outm1$id)

source('~/documents/R_progs/hazardratedecomp/model_setup.R')
source('~/documents/R_progs/hazardratedecomp/decomp_functions.R')

# two decomposions: one from each perspective (to be averaged)
C1 <- decomp.cthaz(m1stuff, m0stuff, printit=TRUE, scale=1000)
C2 <- decomp.cthaz(m0stuff, m1stuff, printit=TRUE, scale=1000)
# D1 <- decomp.dtlogit(m1stuff, m0stuff, printit=TRUE, scale=100)
# D2 <- decomp.dtlogit(m0stuff, m1stuff, printit=TRUE, scale=100)
# Dc <- decomp.dtcloglog(m1stuff, m0stuff, printit=TRUE, scale=100)

# average
ave.decompE <- function(D1, D2) {
  s  <- 1000
  vn <- D1$vnames
  b1 <- D1$bE
  b2 <- D2$bE
  v1 <- D1$varbE
  v2 <- D2$varbE
  b <- s * (b1 - b2)/2
  v <- s^2 * (v1 + v2)/4
  se.b <- sqrt(diag(v))
  #colnames(b, vn)
  return(data.frame(term=vn, b=b, se=se.b, z=b/se.b, 
                    b.lower = b-1.96*se.b,
                    b.upper = b+1.96*se.b))
}

ave.decompC <- function(D1, D2) {
  s  <- 1000
  vn <- D1$vnames
  b1 <- D1$bC
  b2 <- D2$bC
  v1 <- D1$varbC
  v2 <- D2$varbC
  b <- s * (b1 - b2)/2
  v <- s^2 * (v1 + v2)/4
  se.b <- sqrt(diag(v))
  #colnames(b,vn)
  return(data.frame(term = vn, b=b, se=se.b, z=b/se.b, 
                    b.lower = b-1.96*se.b,
                    b.upper = b+1.96*se.b))
}
E.ave <- ave.decompE(C1,C2)
C.ave <- ave.decompC(C1,C2)

# simple tables
E.ave[, -c(1)]
C.ave[, -c(1)]
# coef plots
pE <- ggplot(E.ave, aes(b, term)) +
         geom_point(color="seagreen") +
         geom_errorbarh(aes(xmin=b.lower, xmax=b.upper), 
                 linewidth=.3,
                 height=.3,
                 linetype=1,
                 color="turquoise4") +
  # add in a dotted line at zero
  geom_vline(xintercept = 0, lty = 2, color="green") +
  theme_gray() +
  labs(
    x = "E-coefs",
    y = NULL,
    title = "Plot of Characteristic Effects (averaged)"
  ) 

pC <- ggplot(C.ave, aes(b, term)) +
  geom_point(color="seagreen") +
  geom_errorbarh(aes(xmin=b.lower, xmax=b.upper), 
                 linewidth=.3,
                 height=.3,
                 linetype=1,
                 color="turquoise4") +
  # add in a dotted line at zero
  geom_vline(xintercept = 0, lty = 2, color="green") +
  theme_gray() +
  labs(
    x = "C-coefs",
    y = NULL,
    title = "Plot of Coefficient Effects (averaged)"
  ) 

gridExtra::grid.arrange(pE,pC, nrow=1, ncol=2)


