require(tidyverse)
require(survey)
############################################
# process model calls:
############################################
mod.out <- function(m,sub) {
  id   <- unlist(sub[[1]])
  # cat(length(id),"\n")
  w   <- weights(sub)
  # cat(length(w),"\n")
  b   <- coef(m)
  x   <- m$x
  v   <- m$cov.unscaled
  off <- m$offset
  if (is.null(off)) off <- rep(0,nrow(x))
  sub <- sub
  mget(ls())
} 

decomp.pwcexp  <- function(form, designA, designB, scale=NULL, printit=FALSE, reverse=FALSE) {
  A <- substitute(svyglm(terms(form), family=quasipoisson(), design = designA, 
                         x=TRUE))
  B <- substitute(svyglm(terms(form), family=quasipoisson(), design = designB,
                         x=TRUE))
  
  if (reverse) outA <- mod.out(eval(B), Bsub)
      else outA     <- mod.out(eval(A), Asub)
  if (reverse) outB <- mod.out(eval(A), Asub)
      else outB     <- mod.out(eval(B), Bsub)
 
  Alist <- list(id=outA$id, b=outA$b,x=outA$x, v=outA$v, w=outA$w, sub=outA$sub, off=outA$off)
  Blist <- list(id=outB$id, b=outB$b,x=outB$x, v=outB$v, w=outB$w, sub=outB$sub, off=outB$off)

decomp_cthaz(Alist, Blist, scale, printit)
}



decomp.dtlogit  <- function(form, designA, designB, scale=NULL, printit=FALSE, reverse=FALSE) {
  A <- substitute(svyglm(terms(form), family=quasibinomial(link=logit), design = designA, 
                         x=TRUE))
  B <- substitute(svyglm(terms(form), family=quasibinomial(link=logit), design = designB,
                         x=TRUE))

  if (reverse) outA <- mod.out(eval(B), Bsub)
      else outA     <- mod.out(eval(A), Asub)
  if (reverse) outB <- mod.out(eval(A), Asub)
      else outB     <- mod.out(eval(B), Bsub)
  
  Alist <- list(id=outA$id, b=outA$b,x=outA$x, v=outA$v, w=outA$w, sub=outA$sub, off=outA$off)
  Blist <- list(id=outB$id, b=outB$b,x=outB$x, v=outB$v, w=outB$w, sub=outB$sub, off=outB$off)
  
  decomp_dtlogit(Alist, Blist, scale, printit)
}

decomp.dtcloglog  <- function(form, designA, designB, scale=NULL, printit=FALSE, reverse=FALSE) {
  A <- substitute(svyglm(terms(form), family=quasibinomial(link=cloglog), design = designA, 
                         x=TRUE))
  B <- substitute(svyglm(terms(form), family=quasibinomial(link=cloglog), design = designB,
                         x=TRUE))

  if (reverse) outA <- mod.out(eval(B),Bsub)
      else outA     <- mod.out(eval(A),Asub)
  if (reverse) outB <- mod.out(eval(A),Asub)
      else outB     <- mod.out(eval(B), Bsub)
  
  Alist <- list(id=outA$id, b=outA$b,x=outA$x, v=outA$v, w=outA$w, sub=outA$sub, off=outA$off)
  Blist <- list(id=outB$id, b=outB$b,x=outB$x, v=outB$v, w=outB$w, sub=outB$sub, off=outB$off)
 
  decomp_dtcloglog(Alist, Blist, scale, printit)
}

############################################################
#               specific model processing
############################################################
############################################################
# functions for piecewise-constant exponential hazard model
############################################################
decomp_cthaz <- function(A, B, scale=NULL, printit) {
  
  if (is.null(scale)) scale <- 1  
  
  mean.fix <- function(df) {
    df %>% group_by(id) %>%
      summarize(across(starts_with('x'),
                       list(m = ~weighted.mean(., weight))),
                .groups = 'drop') %>%  #-> m.i
      summarize(across(starts_with('x'), list(m = ~mean(.,))),
                .groups = 'drop') -> m
    return(unlist(m))
  }
  
  dfA <- data.frame(id=A$id, x=A$x, weight=A$w)
  dfB <- data.frame(id=B$id, x=B$x, weight=B$w)
  
  meanA <- mean.fix(dfA)
  meanB <- mean.fix(dfB)
  
  #
  Asub  <- A$sub
  Bsub  <- B$sub
  bA    <- A$b
  bB    <- B$b
  xA    <- A$x
  xB    <- B$x
  offA  <- A$off
  offB  <- B$off
  varbA <- A$v
  varbB <- B$v
  wA    <- A$w
  wB    <- B$w
  xpA   <- exp(offA) # define exposure for for denominators of rates.
  xpB   <- exp(offB)
  mxpA  <- svymean(~xpA, Asub, na=TRUE) # apply svy weighting
  mxpB  <- svymean(~xpB, Bsub, na=TRUE)
  
  
  
  # F <-  function(b,x,off,wt) {
  #   F <- exp(x%*%b + off)*wt
  #   return(F)
  # }
  # f <-  function(b,x,off,wt) {
  #   f <-  exp(x%*%b + off)*wt
  #   return(f)
  # }
  # 
  # decomp(bA, bB, xA, xB, offA, offB, varbA, varbB, 
  #        meanA, meanB, wA, wB, printit, scale, 
  #        eval(F), eval(f) )
  # decomp(bA, bB, xA, xB, offA, offB, varbA, varbB, 
  #        meanA, meanB, wA, wB, A, B, printit, scale, 
  #        eval(F), eval(f) )
  
  # svy opts
  F <-  function(b,x,off) {
    F <- exp(x%*%b + off)
    return(F)
  }
  f <-  function(b,x,off) {
    f <-  exp(x%*%b + off)
    return(f)
  }
  decomp(bA, bB, xA, xB, offA, offB, varbA, varbB, 
         mxpA, mxpB, meanA, meanB, Asub, Bsub, printit, scale, 
         eval(F), eval(f) )
  
}

#
decomp_cthazB <- function(A, B, scale=NULL, printit) {
  
  if (is.null(scale)) scale <- 1  
  
  mean.fix <- function(df) {
    df %>% group_by(id) %>%
      summarize(across(starts_with('x'),
                       list(m = ~weighted.mean(., weight))),
                .groups = 'drop') %>%  #-> m.i
      summarize(across(starts_with('x'), list(m = ~mean(.,))),
                .groups = 'drop') -> m
    return(unlist(m))
  }
  
  dfA <- data.frame(id=A$id, x=A$x, weight=A$w)
  dfB <- data.frame(id=B$id, x=B$x, weight=B$w)
  
  meanA <- mean.fix(dfA)
  meanB <- mean.fix(dfB)
  
  #
  Asub  <- A$sub
  Bsub  <- B$sub
  bA    <- A$b
  bB    <- B$b
  xA    <- A$x
  xB    <- B$x
  offA  <- A$off
  offB  <- B$off
  varbA <- A$v
  varbB <- B$v
  wA    <- A$w
  wB    <- B$w
  xpA   <- exp(offA) # define exposure for for denominators of rates.
  xpB   <- exp(offB)
  mxpA  <- svymean(~xpA, Asub, na=TRUE) # apply svy weighting
  mxpB  <- svymean(~xpB, Bsub, na=TRUE)
  
  
  
  F <-  function(b,x,off,wt) {
    F <- exp(x%*%b + off)*wt
    return(F)
  }
  f <-  function(b,x,off,wt) {
    f <-  exp(x%*%b + off)*wt
    return(f)
  }

  decompB(bA, bB, xA, xB, offA, offB, varbA, varbB,
         meanA, meanB, wA, wB, printit, scale,
         eval(F), eval(f) )
 
  # # svy opts
  # F <-  function(b,x,off) {
  #   F <- exp(x%*%b + off)
  #   return(F)
  # }
  # f <-  function(b,x,off) {
  #   f <-  exp(x%*%b + off)
  #   return(f)
  # }
  # decomp(bA, bB, xA, xB, offA, offB, varbA, varbB, 
  #        mxpA, mxpB, meanA, meanB, Asub, Bsub, printit, scale, 
  #        eval(F), eval(f) )
  
}
################################################
# functions for discrete-time logit hazard model
################################################
decomp_dtlogit <- function(A,  B, scale=NULL, printit) {
  
  if (is.null(scale)) scale <- 1
  
  mean.fix <- function(df) {
    df %>% group_by(id) %>%
      summarize(across(starts_with('x'),
                       list(m = ~weighted.mean(., weight))),
                .groups = 'drop') %>%  #-> m.i
      summarize(across(starts_with('x'), list(m = ~mean(.,))),
                .groups = 'drop') -> m
    return(unlist(m))
  }
  
  dfA <- data.frame(id=A$id, x=A$x, weight=A$w)
  dfB <- data.frame(id=B$id, x=B$x, weight=B$w)
  
  meanA <- mean.fix(dfA)
  meanB <- mean.fix(dfB)
  
  #
  bA    <- A$b
  bB    <- B$b
  xA    <- A$x
  xB    <- B$x
  offA  <- A$off
  offB  <- B$off
  varbA <- A$v
  varbB <- B$v
  wA    <- A$w
  wB    <- B$w
  
  F <-  function(b,x,off,wt) {
    F <- wt * exp(x%*%b + off)/(1  + exp(x%*%b + off))
    return(F)
  }
  f <-  function(b,x,off,wt) {
    f <-  wt * (exp(x%*%b + off)/(1 + exp(x%*%b + off)))^2
    return(f)
  }
  
  decompB(bA, bB, xA, xB, offA, offB, varbA, varbB, 
         meanA, meanB, wA, wB, printit, scale, 
         eval(F), eval(f) )
  
  # decomp(bA, bB, xA, xB, offA, offB, varbA, varbB, mxpA, mxpB,
  #        meanA, meanB, Asub, Bsub, printit, scale,
  #        eval(F), eval(f) )
  
  
}

##################################################
# functions for discrete time cloglog hazard model
##################################################
decomp_dtcloglog <- function(A,  B, scale=NULL, printit) {
  
  mean.fix <- function(df) {
    df %>% group_by(id) %>%
      summarize(across(starts_with('x'),
                       list(m = ~weighted.mean(., weight))), 
                .groups = 'drop') %>%  #-> m.i
      summarize(across(starts_with('x'), list(m = ~mean(.,))),
                .groups = 'drop') -> m
    return(unlist(m))
  }
  
  dfA <- data.frame(id=A$id, x=A$x, weight=A$w)
  dfB <- data.frame(id=B$id, x=B$x, weight=B$w)
  
  meanA <- mean.fix(dfA)
  meanB <- mean.fix(dfB)
  
  #
  bA    <- A$b
  bB    <- B$b
  xA    <- A$x
  xB    <- B$x
  offA  <- A$off
  offB  <- B$off
  varbA <- A$v
  varbB <- B$v
  wA    <- A$w
  wB    <- B$w
  
  
  if (is.null(scale)) scale <- 1  
  
  # distribution functions
  
  F <-  function(b,x,off,wt) {
    F <- wt * (1 - exp(-exp(x%*%b + off)))
    return(F)
  }
  f <-  function(b,x,off,wt) {
    f <-  wt * exp(-exp(x%*%b + off))
    return(f)
  }
  
  
  decompB(bA, bB, xA, xB, offA, offB, varbA, varbB,
         meanA, meanB, wA, wB, printit, scale,
         eval(F), eval(f) )
  
  # decomp(bA, bB, xA, xB, offA, offB, varbA, varbB, mxpA, mxpB,
  #        meanA, meanB, Asub, Bsub, printit, scale,
  #        eval(F), eval(f) )
  
  
}

####################   M A I N  ##################################
#  at present these exploit svy and give samme results as decompB
##################################################################

decomp <- function(bA, bB, xA, xB, offA, offB, varbA, varbB, 
                     mxpA, mxpB, meanA, meanB, Asub, Bsub, printit, scale, F, f) {
    
F <- eval(F)
f <- eval(f)

mA <- meanA
mB <- meanB
 
# Yun weight function (composition)
Wdx.F  <- function(b,x1,x2){
  A <- (x1-x2)%*%b
  Wdx <- NULL
  for (i in 1:length(b)){
    Wdx[i] <- (x1[i] - x2[i])*b[i] / A
  }
  return(Wdx)
  
}
# Yun weight function (coefficient)
Wdb.F  <- function(b1,b2,x){
  A <- x%*%(b1-b2)
  Wdb <- NULL
  for (i in 1:length(x)){
    Wdb[i] <- (x[i]*(b1[i] - b2[i])) / A
  }
  return(Wdb)
}

dW.F <- function(b,x1,x2) {
  dW <- array(NA,c(length(b), length(b)))
  A <- (x1-x2)%*%b   
  for (k in 1:length(b)) {
    for (l in 1:length(b)) {
      dW[k,l]  <-  as.numeric(k==l) * ((x1[k] - x2[k])/A) - 
        (b[k] * ( (x1[k] - x2[k])*(x1[l] - x2[l]) )/A^2 )
    }
  }
  return(dW)
}     

dwA.F <- function(b1,b2,x2) {
  # derivative of Wdb part A #  NEW  (this is now a K x K result) 
  dwA1 <- array(NA, c(length(b1), length(b2)))
  A <- x2%*%(b1-b2)
  for (k in 1:length(b1)){
    for (l in 1:length(b2)) {
      dwA1[k,l] <- as.numeric(k==l) * x2[k]/A - (x2[k]*x2[l]*(b1[k]-b2[k]))/A^2 
    } 
  }
  return(dwA1)
}

dwB.F <- function(b1,b2,x2) { 
  # derivative of Wdb part B
  dwB1 <- array(NA, c(length(b1), length(b2)))
  A <- x2%*%(b1-b2)
  for (k in 1:length(b1)){
    for (l in 1:length(b2)) {
      dwB1[k,l] <- (x2[k]*x2[l]*(b1[k]-b2[k]))/A^2 - as.numeric(k==l) * (x2[k]/A) 
    } 
  }
  return(dwB1)
}

wbA  <-  dwA.F(bA,bB,mB)
wbB  <-  dwB.F(bA,bB,mB)

# weights
Wdx  <-  Wdx.F(bA, mA, mB)
Wdb  <-  Wdb.F(bA, bB, mB)

# i=x, j=b i,k {0,1}
# f00 <- f(bA,xA,offA,wA)
# f01 <- f(bA,xB,offB,wB)
# f11 <- f(bB,xB,offB,wB)
# F00 <- sum(F(bA, xA, offA, wA))/sum(exp(offA)*wA)
# F01 <- sum(F(bA, xB, offB, wB))/sum(exp(offB)*wB)
# F11 <- sum(F(bB, xB, offB, wB))/sum(exp(offB)*wB)

f00 <- f(bA,xA,offA)
f01 <- f(bA,xB,offB)
f11 <- f(bB,xB,offB)
F00 <- svymean(~F(bA, xA, offA),Asub)/mxpA
F01 <- svymean(~F(bA, xB, offB),Bsub)/mxpB
F11 <- svymean(~F(bB, xB, offB),Bsub)/mxpB

#Convention: Yhi - Ylo = 1st moment higher group - 1st moment lower group
# decomp composition: E = F(Bhi,Xhi) - F(Bhi,Xlo) 
E   <-  F00 - F01

#Coefficients + Unexplained: C = F(Bhi,Xlo) - F(Blo,Xlo) 
C   <-  F01 - F11

# get dEdb for variance estimator
dWx  <- dW.F(bA, mA, mB)

# dEdb <- array(NA, c(length(bA), length(bA)))
# for (k in 1:length(bA)){
#   for (l in 1:length(bA)) {
#     dEdb[k,l] <- Wdx[k]*(sum(f00*xA[,l])/sum(exp(offA)*wA) -  
#                          sum(f01*xB[,l])/sum(exp(offB)*wB))  + 
#       dWx[k,l]*E
#   }
# }

dEdb <- array(NA, c(length(bA), length(bA)))
for (k in 1:length(bA)) {
  for (l in 1:length(bA)) {
    # 
    df <- svymean(~as.vector(f00*xA[,l]), Asub)/mxpA - 
          svymean(~as.vector(f01*xB[,l]), Bsub)/mxpB
    dEdb[k,l] <- Wdx[k] * df + dWx[k,l] * E
  }
}
# 
# dCdbA <- array(NA, c(length(bA), length(bA)))
# dCdbB <- array(NA, c(length(bB), length(bB)))
# for (k in 1:length(bA)){
#   for (l in 1:length(bB)) { 
#     dCdbA[k,l] <- Wdb[k]*sum(f01*xB[,l])/sum(exp(offA)*wA) + wbA[k,l]*C
#     dCdbB[k,l] <- wbB[k,l]*C - Wdb[k]*sum(f11*xB[,l])/sum(exp(offB)*wB) 
#   }
# } 

dCdbA   <- array(NA, c(length(bA), length(bA)))
dCdbB   <- array(NA, c(length(bB), length(bB)))
for (k in 1:length(bA)){
  for (l in 1:length(bB)) { 
    dCdbA[k,l] <- Wdb[k] * svymean(~as.vector(f01*xB[,l]), Bsub)/mxpB + wbA[k,l]*C
    dCdbB[k,l] <- wbB[k,l] * C - Wdb[k]*svymean(~as.vector(f11*xB[,l]),Bsub)/mxpB
  }
} 


##   Variances
#    Composition (E.k)
K <- length(bA)

Var.E.k <- dEdb%*%varbA%*%t(dEdb)
seWdx <- sqrt(diag(Var.E.k))
##  Variances
#   Coefficients (C.k)

Var.C.k <- dCdbA%*%varbA%*%t(dCdbA) + dCdbB%*%varbB%*%t(dCdbB)
seWdb <- sqrt(diag(Var.C.k))

# get asymptotic variance of overall E and C
# E
# dEdb.0 <- NULL
# for (k in 1:length(bA)){
#   dEdb.0[k] <- sum(f00*xA[,k])/sum(exp(offA)*wA) - 
#     sum(f01*xB[,k])/sum(exp(offB)*wB)
# }

dEdb.0 <- NULL
for (k in 1:length(bA)){
  dEdb.0[k] <- svymean(~as.vector(f00*xA[,k]), Asub)/mxpA - 
               svymean(~as.vector(f01*xB[,k]), Bsub)/mxpB
}

# C
# dCdb.0A <- NULL
# for (k in 1:length(bA)){
#   dCdb.0A[k] <- sum(f01*xB[,k])/sum(exp(offB)*wB)
# }
# dCdb.0B <- NULL
# for (k in 1:length(bB)){
#   dCdb.0B[k] <- sum(f11*xB[,k])/sum(exp(offB)*wB)
# }
#C_k
dCdb.0A <- NULL
for (k in 1:length(bA)){
  dCdb.0A[k] <- svymean(~as.vector(f01*xB[,k]),Bsub)/mxpB
}
dCdb.0B <- NULL
for (k in 1:length(bB)){
  dCdb.0B[k] <- svymean(~as.vector(f11*xB[,k]),Bsub)/mxpB
}

#         1 x k   k x k     k x 1

var.E0 <- t(dEdb.0)%*%varbA%*%dEdb.0
s.E0   <- sqrt(var.E0)

var.C0 <- t(dCdb.0A)%*%varbA%*%dCdb.0A + t(dCdb.0B)%*%varbB%*%dCdb.0B 
s.C0   <- sqrt(var.C0)

R <- E + C
v.names <- names(bA)
if (is.null(v.names)) {
  v.names <- paste("Var", 1:length(bA), sep = "") # if names not available
}
ymA <- F00 * scale
ymB <- F11 * scale

if (printit==TRUE) {
  cat("","\n")
  cat("Hazard Rate Decomposition","\n")
  cat("","\n")
  cat("Rate Group A           ", formatC(ymA, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Rate Group B           ", formatC(ymB, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("__________________________________________________________________", "\n")
  cat("Total Difference (E + C)", formatC((E+C)*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Raw Contrib......(E,  C)", formatC(E*scale, flag=" ", dig=4, wid=10, format="f"), 
      formatC(C*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Std. Errors......(E,  C)", formatC(s.E0*scale, flag=" ", dig=4, wid=10, format="f"),
      formatC(s.C0*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Z-values.........(E,  C)", formatC(E/s.E0, flag=" ", dig=4, wid=10, format="f"),
      formatC(C/s.C0, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Pct Contrib......(E,  C)", formatC(100*E/(E+C), flag=" ", dig=4, wid=10, format="f"),
      formatC(100*C/(E+C), flag=" ", dig=4, wid=10, format="f"), "\n" )
  cat("","\n")
  cat("SUMMARY MEASURES","\n")
  cat("        Comparison Group", "\n")
  cat("        ----------------", "\n")
  cat("         Variable  ", "      Mean  ",   "Coefficient", "Std. Error", "Z-value", "\n")
  cat("__________________________________________________________________", "\n")
  for (i in 1:length(mA)){
    cat(formatC(v.names[i],format="s", wid=20),
        formatC(mA[i],flag=" ",dig=4,wid=10, format="f"),
        formatC(bA[i],wid=10,dig=4, format="f"),
        formatC(sqrt(diag(varbA))[i],wid=10,dig=4,format="f"), 
        formatC(bA[i]/sqrt(diag(varbA))[i],wid=10,dig=4,format="f"), "\n")
  }
  cat("","\n")
  cat("        Reference Group", "\n")
  cat("        ---------------", "\n")
  cat("         Variable  ", "      Mean  ",   "Coefficient", "Std. Error", "Z-value", "\n")
  cat("__________________________________________________________________", "\n")
  for (i in 1:length(mB)){
    cat(formatC(v.names[i],format="s", wid=20),
        formatC(mB[i],flag=" ",   dig=4,wid=10, format="f"),
        formatC(bB[i],wid=10,     dig=4, format="f"),
        formatC(sqrt(diag(varbB))[i],wid=10,dig=4,format="f"), 
        formatC(bB[i]/sqrt(diag(varbB))[i],wid=10,dig=4,format="f"), "\n")
  }
  cat("","\n")
  cat("DUE TO DIFFERENCE IN COMPOSITION (E)", "\n")
  cat("","\n")
  cat("          Variable  ", " Estimate",   "  Std. Error",  "  Z-Val",  "        PCT",  "\n")
  cat("___________________________________________________________________", "\n")
  for (i in 1:length(bA)){
    cat(formatC(v.names[i],format="s", wid=20),
        formatC(E*Wdx[i]*scale,       flag=" ", dig=7, wid=10, format="f"),
        formatC(seWdx[i]*scale,       flag=" ", dig=7, wid=10, format="f"),
        formatC(E*Wdx[i]/seWdx[i],    flag=" ", dig=4, wid=10, format="f"),
        formatC(100*(E*Wdx[i]/(E+C)), flag="+", dig=1, wid=10, format="f"), "\n")
  }
  cat("","\n")
  cat("DUE TO DIFFERENCE IN COEFFICIENTS (C)", "\n")
  cat("","\n")
  cat("          Variable  ", " Estimate",   "  Std. Error",  "  Z-Val",  "         PCT"  ,"\n")
  cat("___________________________________________________________________", "\n")
  for (i in 1:length(bA)){
    cat(formatC(v.names[i],format="s", wid=20),
        formatC(C*Wdb[i]*scale,       flag=" ",  dig=7, wid=10, format="f"),
        formatC(seWdb[i]*scale,       flag=" ",  dig=7, wid=10, format="f"),
        formatC(C*Wdb[i]/seWdb[i],    flag=" ",  dig=4, wid=10, format="f"), 
        formatC(100*(C*Wdb[i]/(E+C)), flag="+",  dig=1, wid=10, format="f"),"\n")
    
  }
} # end printit

K        <- length(bA)
# store these
bE       <- matrix(E*Wdx,K,1) 
varbE    <- matrix(Var.E.k,K,K)
bC       <- matrix(C*Wdb,K,1)
varbC    <- matrix(Var.C.k,K,K)
E        <- E
varE     <- var.E0
C        <- C
varC     <- var.C0
rownames(bE) <- v.names
rownames(bC) <- v.names
colnames(varbE) <- rownames(varbE) <- v.names
colnames(varbC) <- rownames(varbC) <- v.names

# returns 
return(list(vnames=v.names, E=E, varE=varE, C=C, varC=varC,
            bE=bE, varbE=varbE, 
            bC=bC, varbC=varbC))

  } # end decomp



#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
######################   M A I N  ####################################
#  at present these do not exploit svy but give same results as decomp    
######################################################################
decompB <- function(bA, bB, xA, xB, offA, offB, varbA, varbB, 
                   meanA, meanB, wA, wB, printit, scale, F, f) {
  
  F <- eval(F)
  f <- eval(f)
  
  mA <- meanA
  mB <- meanB
  
  # Yun weight function (composition)
  Wdx.F  <- function(b,x1,x2){
    A <- (x1-x2)%*%b
    Wdx <- NULL
    for (i in 1:length(b)){
      Wdx[i] <- (x1[i] - x2[i])*b[i] / A
    }
    return(Wdx)
    
  }
  # Yun weight function (coefficient)
  Wdb.F  <- function(b1,b2,x){
    A <- x%*%(b1-b2)
    Wdb <- NULL
    for (i in 1:length(x)){
      Wdb[i] <- (x[i]*(b1[i] - b2[i])) / A
    }
    return(Wdb)
  }
  
  dW.F <- function(b,x1,x2) {
    dW <- array(NA,c(length(b), length(b)))
    A <- (x1-x2)%*%b   
    for (k in 1:length(b)) {
      for (l in 1:length(b)) {
        dW[k,l]  <-  as.numeric(k==l) * ((x1[k] - x2[k])/A) - 
          (b[k] * ( (x1[k] - x2[k])*(x1[l] - x2[l]) )/A^2 )
      }
    }
    return(dW)
  }     
  
  dwA.F <- function(b1,b2,x2) {
    # derivative of Wdb part A #  NEW  (this is now a K x K result) 
    dwA1 <- array(NA, c(length(b1), length(b2)))
    A <- x2%*%(b1-b2)
    for (k in 1:length(b1)){
      for (l in 1:length(b2)) {
        dwA1[k,l] <- as.numeric(k==l) * x2[k]/A - (x2[k]*x2[l]*(b1[k]-b2[k]))/A^2 
      } 
    }
    return(dwA1)
  }
  
  dwB.F <- function(b1,b2,x2) { 
    # derivative of Wdb part B
    dwB1 <- array(NA, c(length(b1), length(b2)))
    A <- x2%*%(b1-b2)
    for (k in 1:length(b1)){
      for (l in 1:length(b2)) {
        dwB1[k,l] <- (x2[k]*x2[l]*(b1[k]-b2[k]))/A^2 - as.numeric(k==l) * (x2[k]/A) 
      } 
    }
    return(dwB1)
  }
  
  wbA  <-  dwA.F(bA,bB,mB)
  wbB  <-  dwB.F(bA,bB,mB)
  
  # weights
  Wdx  <-  Wdx.F(bA, mA, mB)
  Wdb  <-  Wdb.F(bA, bB, mB)
  
  # i=x, j=b i,k {0,1}
  f00 <- f(bA,xA,offA,wA)
  f01 <- f(bA,xB,offB,wB)
  f11 <- f(bB,xB,offB,wB)
  F00 <- sum(F(bA, xA, offA, wA))/sum(exp(offA)*wA)
  F01 <- sum(F(bA, xB, offB, wB))/sum(exp(offB)*wB)
  F11 <- sum(F(bB, xB, offB, wB))/sum(exp(offB)*wB)


  #Convention: Yhi - Ylo = 1st moment higher group - 1st moment lower group
  # decomp composition: E = F(Bhi,Xhi) - F(Bhi,Xlo) 
  E   <-  F00 - F01
  
  #Coefficients + Unexplained: C = F(Bhi,Xlo) - F(Blo,Xlo) 
  C   <-  F01 - F11
  
  # get dEdb for variance estimator
  dWx  <- dW.F(bA, mA, mB)
  
  dEdb <- array(NA, c(length(bA), length(bA)))
   for (k in 1:length(bA)){
     for (l in 1:length(bA)) {
       dEdb[k,l] <- Wdx[k]*(sum(f00*xA[,l])/sum(exp(offA)*wA) -
                            sum(f01*xB[,l])/sum(exp(offB)*wB))  +
         dWx[k,l]*E
     }
   }


  #
   dCdbA <- array(NA, c(length(bA), length(bA)))
   dCdbB <- array(NA, c(length(bB), length(bB)))
   for (k in 1:length(bA)){
     for (l in 1:length(bB)) {
       dCdbA[k,l] <- Wdb[k]*sum(f01*xB[,l])/sum(exp(offA)*wA) + wbA[k,l]*C
       dCdbB[k,l] <- wbB[k,l]*C - Wdb[k]*sum(f11*xB[,l])/sum(exp(offB)*wB)
     }
   }
  
  
  
  ##   Variances
  #    Composition (E.k)
  K <- length(bA)
  
  Var.E.k <- dEdb%*%varbA%*%t(dEdb)
  seWdx <- sqrt(diag(Var.E.k))
  ##  Variances
  #   Coefficients (C.k)
  
  Var.C.k <- dCdbA%*%varbA%*%t(dCdbA) + dCdbB%*%varbB%*%t(dCdbB)
  seWdb <- sqrt(diag(Var.C.k))
  
  # get asymptotic variance of overall E and C
  # E
  dEdb.0 <- NULL
   for (k in 1:length(bA)){
     dEdb.0[k] <- sum(f00*xA[,k])/sum(exp(offA)*wA) - 
       sum(f01*xB[,k])/sum(exp(offB)*wB)
   }
  
  
  
  # C
   dCdb.0A <- NULL
    for (k in 1:length(bA)){
      dCdb.0A[k] <- sum(f01*xB[,k])/sum(exp(offB)*wB)
    }
    dCdb.0B <- NULL
    for (k in 1:length(bB)){
      dCdb.0B[k] <- sum(f11*xB[,k])/sum(exp(offB)*wB)
    }
   
  #         1 x k   k x k     k x 1
  
  var.E0 <- t(dEdb.0)%*%varbA%*%dEdb.0
  s.E0   <- sqrt(var.E0)
  
  var.C0 <- t(dCdb.0A)%*%varbA%*%dCdb.0A + t(dCdb.0B)%*%varbB%*%dCdb.0B 
  s.C0   <- sqrt(var.C0)
  
  R <- E + C
  v.names <- names(bA)
  if (is.null(v.names)) {
    v.names <- paste("Var", 1:length(bA), sep = "") # if names not available
  }
  ymA <- F00 * scale
  ymB <- F11 * scale
  
  if (printit==TRUE) {
    cat("","\n")
    cat("Hazard Rate Decomposition","\n")
    cat("","\n")
    cat("Rate Group A           ", formatC(ymA, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Rate Group B           ", formatC(ymB, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("__________________________________________________________________", "\n")
    cat("Total Difference (E + C)", formatC((E+C)*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Raw Contrib......(E,  C)", formatC(E*scale, flag=" ", dig=4, wid=10, format="f"), 
        formatC(C*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Std. Errors......(E,  C)", formatC(s.E0*scale, flag=" ", dig=4, wid=10, format="f"),
        formatC(s.C0*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Z-values.........(E,  C)", formatC(E/s.E0, flag=" ", dig=4, wid=10, format="f"),
        formatC(C/s.C0, flag=" ", dig=4, wid=10, format="f"), "\n")
    cat("Pct Contrib......(E,  C)", formatC(100*E/(E+C), flag=" ", dig=4, wid=10, format="f"),
        formatC(100*C/(E+C), flag=" ", dig=4, wid=10, format="f"), "\n" )
    cat("","\n")
    cat("SUMMARY MEASURES","\n")
    cat("        Comparison Group", "\n")
    cat("        ----------------", "\n")
    cat("         Variable  ", "      Mean  ",   "Coefficient", "Std. Error", "Z-value", "\n")
    cat("__________________________________________________________________", "\n")
    for (i in 1:length(mA)){
      cat(formatC(v.names[i],format="s", wid=20),
          formatC(mA[i],flag=" ",dig=4,wid=10, format="f"),
          formatC(bA[i],wid=10,dig=4, format="f"),
          formatC(sqrt(diag(varbA))[i],wid=10,dig=4,format="f"), 
          formatC(bA[i]/sqrt(diag(varbA))[i],wid=10,dig=4,format="f"), "\n")
    }
    cat("","\n")
    cat("        Reference Group", "\n")
    cat("        ---------------", "\n")
    cat("         Variable  ", "      Mean  ",   "Coefficient", "Std. Error", "Z-value", "\n")
    cat("__________________________________________________________________", "\n")
    for (i in 1:length(mB)){
      cat(formatC(v.names[i],format="s", wid=20),
          formatC(mB[i],flag=" ",   dig=4,wid=10, format="f"),
          formatC(bB[i],wid=10,     dig=4, format="f"),
          formatC(sqrt(diag(varbB))[i],wid=10,dig=4,format="f"), 
          formatC(bB[i]/sqrt(diag(varbB))[i],wid=10,dig=4,format="f"), "\n")
    }
    cat("","\n")
    cat("DUE TO DIFFERENCE IN COMPOSITION (E)", "\n")
    cat("","\n")
    cat("          Variable  ", " Estimate",   "  Std. Error",  "  Z-Val",  "        PCT",  "\n")
    cat("___________________________________________________________________", "\n")
    for (i in 1:length(bA)){
      cat(formatC(v.names[i],format="s", wid=20),
          formatC(E*Wdx[i]*scale,       flag=" ", dig=7, wid=10, format="f"),
          formatC(seWdx[i]*scale,       flag=" ", dig=7, wid=10, format="f"),
          formatC(E*Wdx[i]/seWdx[i],    flag=" ", dig=4, wid=10, format="f"),
          formatC(100*(E*Wdx[i]/(E+C)), flag="+", dig=1, wid=10, format="f"), "\n")
    }
    cat("","\n")
    cat("DUE TO DIFFERENCE IN COEFFICIENTS (C)", "\n")
    cat("","\n")
    cat("          Variable  ", " Estimate",   "  Std. Error",  "  Z-Val",  "         PCT"  ,"\n")
    cat("___________________________________________________________________", "\n")
    for (i in 1:length(bA)){
      cat(formatC(v.names[i],format="s", wid=20),
          formatC(C*Wdb[i]*scale,       flag=" ",  dig=7, wid=10, format="f"),
          formatC(seWdb[i]*scale,       flag=" ",  dig=7, wid=10, format="f"),
          formatC(C*Wdb[i]/seWdb[i],    flag=" ",  dig=4, wid=10, format="f"), 
          formatC(100*(C*Wdb[i]/(E+C)), flag="+",  dig=1, wid=10, format="f"),"\n")
      
    }
  } # end printit
  
  K        <- length(bA)
  # store these
  bE       <- matrix(E*Wdx,K,1) 
  varbE    <- matrix(Var.E.k,K,K)
  bC       <- matrix(C*Wdb,K,1)
  varbC    <- matrix(Var.C.k,K,K)
  E        <- E
  varE     <- var.E0
  C        <- C
  varC     <- var.C0
  rownames(bE) <- v.names
  rownames(bC) <- v.names
  colnames(varbE) <- rownames(varbE) <- v.names
  colnames(varbC) <- rownames(varbC) <- v.names
  
  # returns 
  return(list(vnames=v.names, E=E, varE=varE, C=C, varC=varC,
              bE=bE, varbE=varbE, 
              bC=bC, varbC=varbC))
  
} # end decompB
