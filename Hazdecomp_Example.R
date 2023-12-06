# HazdecompExample.R 
rm(list=ls())
setwd("/Users/dpowers/Documents/R_progs/hazardratedecomp")

require(tidyverse)
require(survey)


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
# (mean-centered) : 
#  pctsmom  = proportion of years from age 0-18 spent in single mother household
#  nfamtran = number of family-structure changes from 0-18
#  medu     = respondent mother's years of schooling)
#  inc1000  = adjusted family income in thousands of $
#  nosibs   = number of older siblings
#  magebir  = respondent mother's age at respondent's birth)
#  age      = factor for age-intervals 
#  race     = factor for race
#  famid    = family identifier for sibling clusters
#  iid      = individual identifier
#############################################################################
df <- read.csv(file='hazdata.csv')


# libraries
require(survey)

# set survey design: required regardless of actual design: 
# clusters, no weights
dfsvy <- svydesign(ids=~iid, weights=~1, cluster=~famid, nest=TRUE, data=df)
# no weights, no clusters 
dfsvy <- svydesign(ids=~iid, weights=~1, nest=TRUE, data=df)
# OK (ignore warning)
dfsvy <- svydesign(ids=~iid, data=df)
# make subsets etc.
# check on values of grouping variable
table(df$race) #OK
# define subsets here
Asub <- subset(dfsvy, race == "Black")
Bsub <- subset(dfsvy, race == "White")
#
# source main subroutines
source('hazard_decomp_functions_svy.R')
# rate decomposition
#
# test calls
m1 <- decomp.pwcexp(devnt ~ age + pctsmom + nfamtran + medu + 
                      inc1000 + nosibs + magebir + offset(logexp) - 1,
                    Asub, Bsub, scale=1000, printit=TRUE)

m2 <- decomp.dtlogit(devnt ~ age + pctsmom + nfamtran + medu + 
                       inc1000 + nosibs + magebir  - 1,
                     Asub, Bsub, scale=100, printit=TRUE)

m2 <- decomp.dtcloglog(devnt ~ age + pctsmom + nfamtran + medu + 
                         inc1000 + nosibs + magebir  - 1,
                       Asub, Bsub, scale=100, printit=TRUE)



require(tidyverse)
# average two decompositions: A-B and B-A (reversed)
m1a <- decomp.pwcexp(devnt ~ age + pctsmom + nfamtran + medu + 
                      inc1000 + nosibs + magebir + offset(logexp) - 1,
                    Asub, Bsub, scale=1000, reverse=FALSE)
m1b <- decomp.pwcexp(devnt ~ age + pctsmom + nfamtran + medu + 
                       inc1000 + nosibs + magebir + offset(logexp) - 1,
                     Asub, Bsub, scale=1000, reverse=TRUE)

E.ave <- ave.decompE(m1a, m1b, scale=1000)
C.ave <- ave.decompC(m1a, m1b, scale=1000)

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


# source all below here

ave.decompE <- function(D1, D2, scale) {
  s  <- scale
  vn <- D1$vnames
  b1 <- D1$bE
  b2 <- D2$bE
  v1 <- D1$varbE
  v2 <- D2$varbE
  b <- s * (b1 - b2)/2
  v <-     (v1 + v2)/4
  se.b <- sqrt(diag(v)) * s
  return(data.frame(term=vn, b=b, se=se.b, z=b/se.b, 
                    b.lower = b-1.96*se.b,
                    b.upper = b+1.96*se.b))
}

ave.decompC <- function(D1, D2, scale) {
  s  <- scale
  vn <- D1$vnames
  b1 <- D1$bC
  b2 <- D2$bC
  v1 <- D1$varbC
  v2 <- D2$varbC
  b <- s * (b1 - b2)/2
  v <- s^2 * (v1 + v2)/4
  se.b <- sqrt(diag(v))
  return(data.frame(term = vn, b=b, se=se.b, z=b/se.b, 
                    b.lower = b-1.96*se.b,
                    b.upper = b+1.96*se.b))
}


