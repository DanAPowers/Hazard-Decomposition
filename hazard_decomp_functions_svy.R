############################################################################
# Author: Dan Powers (11/22/23)
#       (updated 10/15/26)
# functions for multivariate rate decomposition
# of difference in hazard rates---E[r(XB)]--between group a and group b: 
#
#  D = E[r(XaBa)] - E[r(XbBb)] 
#    = {E[r(XaBa)]  - E[r(XbBa)]} + {E[r(XbBa)] - E[F(XbBb)]}

#  X -> characteristics (observed)
#  B -> effects         (estimated)
# 
# use svyglm for continuous-time and discrete-time regression models
# compare to Stata -mvdcmp-
######################################################################################################################
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

###########################################
# models
###########################################
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
#           model-specific processing functions
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
         F, f)
  
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
  xpA   <- exp(offA)
  xpB   <- exp(offB)
  mxpA  <- svymean(~xpA, Asub, na=TRUE) # apply svy weighting
  mxpB  <- svymean(~xpB, Bsub, na=TRUE)
  
  F <-  function(b,x,off) {
    F <- exp(x%*%b + off)/(1 + exp(x%*%b + off))
    return(F)
  }
  f <-  function(b,x,off) {
    f <- exp(x%*%b + off)/(1 + exp(x%*%b + off))^2
    return(f)
  }
  
  decomp(bA, bB, xA, xB, offA, offB, varbA, varbB, 
         mxpA, mxpB, meanA, meanB, Asub, Bsub, printit, scale, 
         F, f)
  
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
  xpA   <- exp(offA)
  xpB   <- exp(offB)
  mxpA  <- svymean(~xpA, Asub, na=TRUE) # apply svy weighting
  mxpB  <- svymean(~xpB, Bsub, na=TRUE)
  
  
  if (is.null(scale)) scale <- 1  
  
  # distribution functions
  
  F <-  function(b,x,off) {
    F <- 1 - exp(-exp(x%*%b + off))
    return(F)
  }
  f <-  function(b,x,off) {
    f <-  exp(-exp(x%*%b + off))
    return(f)
  }
  
  
  decomp(bA, bB, xA, xB, offA, offB, varbA, varbB,
         mxpA, mxpB, meanA, meanB, Asub, Bsub, printit, scale,
         F, f)
}

####################   M A I N  ######################################
#        at present these functions exploit svy
######################################################################

decomp <- function(bA, bB, xA, xB, offA, offB, varbA, varbB, 
                     mxpA, mxpB, meanA, meanB, Asub, Bsub, printit, scale, F, f) {
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

dEdb <- array(NA, c(length(bA), length(bA)))
for (k in 1:length(bA)) {
  for (l in 1:length(bA)) {
    # 
    df <- svymean(~as.vector(f00*xA[,l]), Asub)/mxpA - 
          svymean(~as.vector(f01*xB[,l]), Bsub)/mxpB
    dEdb[k,l] <- Wdx[k] * df + dWx[k,l] * E
  }
}


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


dEdb.0 <- NULL
for (k in 1:length(bA)){
  dEdb.0[k] <- svymean(~as.vector(f00*xA[,k]), Asub)/mxpA - 
               svymean(~as.vector(f01*xB[,k]), Bsub)/mxpB
}

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

# collect things to return
K        <- length(bA)
bE       <- matrix(E*Wdx,K,1) 
varbE    <- matrix(Var.E.k,K,K)
bC       <- matrix(C*Wdb,K,1)
varbC    <- matrix(Var.C.k,K,K)
E        <- E[[1]]
varE     <- var.E0
C        <- C[[1]]
varC     <- var.C0
rownames(bE) <- v.names
rownames(bC) <- v.names
colnames(varbE) <- rownames(varbE) <- v.names
colnames(varbC) <- rownames(varbC) <- v.names

# returns to model object
return(list( vnames=v.names, 
             scale=scale,
             coefE=bE, 
             vcovE=varbE, 
             coefC=bC, 
             vcovC=varbC,
             E=E,
             C=C,    
             varE=varE,
             varC=varC )
)



} # end decomp

# helper functions:

# get.decomp processes results from a single decomposition
# ave.decomp averages results from two decompositions
# they return data frames that can be processed easily
# tabFun helps in making tables

get.decomp <- function(m, scale=1, component = NULL) {
  if (is.null(scale)) {
    s = 1
  } else {
    s = scale
  }
  if (is.null(component)) {
    print("component E or C must be specifed")
  }
  
  gTot = m$E + m$C
  vgTot = m$varE + m$varC
  
  if (component == "E") {
    vn = c(m$vnames, "Total", "gTotal")
    b = m$coefE
    v = diag(m$vcovE)
    tot = m$E
    vtot = m$varE
  }
  else {
    vn = c(m$vnames, "Total", "gTotal")
    b = m$coefC
    v = diag(m$vcovC)
    tot = m$C
    vtot = m$varC
  }
  
  b    <- c(b, tot, gTot) * s
  se.b <- sqrt(c(v, vtot, vgTot)) * s
  return(data.frame(term=vn,  b=b, se=se.b, z=b/se.b, 
                    b.lower = b - 1.96*se.b,
                    b.upper = b + 1.96*se.b))
}

ave.decomp <- function(m1, m2, scale=1, component=NULL) {
  if (is.null(scale)) {
    s=1
  } else {
    s=scale
  }
  if (is.null(component)) {
    print("component E or C must be specified") 
    exit
  }
  gAveE  = (m1$E - m2$E)/2
  gAveC  = (m1$C - m2$C)/2
  vgAveE = (m1$varE + m2$varE)/4
  vgAveC = (m1$varC + m2$varC)/4
  gTot   =  gAveE + gAveC
  vgTot  =  vgAveE + vgAveC
  
  if (component == "E") {
    vn = c(m1$vnames, "Total", "gTotal")
    b  = (m1$coefE - m2$coefE)/2
    v  = diag((m1$vcovE + m1$vcovE)/4)
    tot  = gAveE
    vtot  = vgAveE
  }
  else {
    vn = c(m1$vnames, "Total", "gTotal")
    b  = (m1$coefC  - m2$coefC)/2
    v  = diag((m1$vcovC  + m2$vcovC)/4)
    tot  = gAveC
    vtot  = vgAveC
  }
  
  b    <- c(b, tot, gTot) * s
  se.b <- sqrt(c(v, vtot, vgTot)) * s
  return(data.frame(term=vn,  b=b, se=se.b, z=b/se.b, 
                    b.lower = b  - 1.96*se.b,
                    b.upper = b  + 1.96*se.b))
}

# table builder
# takes E and C components (in that order) from get.decomp
# data frame returns a list
tabFun <- function(df1, df2) {
  dfE    <-  dplyr::select(df1, c(b, b.lower, b.upper))
  dfC    <-  dplyr::select(df2, c(b, b.lower, b.upper))
  indx   <- nrow(dfE)
  rownames(dfE) <-  df1$term
  rownames(dfC) <-  df2$term
  dfEp   <-  mutate(dfE,  pct = 100 * b/b[indx]) %>%
    select(c(pct)) %>% 
    filter(row_number() < indx)
  
  dfCp    <- mutate(dfC,  pct = 100 * b/b[indx]) %>%
    select(c(pct)) %>%
    filter(row_number() < indx)
  
  dfE <- dfE %>% filter(row_number() < indx)
  dfC <- dfC %>% filter(row_number() < indx)
  
  dfEc = bind_cols(dfE, dfEp)
  dfCc = bind_cols(dfC, dfCp)
  
  return(list(dfEc, dfCc))
}

# function to summarize output from decomp object

summary.decomp <- function (object,...) 
  
{
  s  <- object$scale
  cnE <- paste0(rownames(object$coefE))
  cnC <- paste0(rownames(object$coefC))
  coefE <- object$coefE  * s
  coefC <- object$coefC  * s
  coefEtot <- object$E * s
  coefCtot <- object$C * s
  covmatE <- object$vcovE
  covmatC <- object$vcovC
  var.cfE <- diag(covmatE)
  var.cfC <- diag(covmatC)
  s.errE  <- sqrt(var.cfE) * s
  s.errC  <- sqrt(var.cfC) * s
  s.errEtot <- sqrt(object$varE) * s
  s.errCtot <- sqrt(object$varC) * s
  tvalueE <- coefE/s.errE 
  tvalueC <- coefC/s.errC 
  tvalueEtot <- coefEtot/s.errEtot
  tvalueCtot <- coefCtot/s.errCtot
  dn      <- c("Estimate", "Std. Error")
  pvalueE <- 2 * pnorm(-abs(tvalueE))
  pvalueC <- 2 * pnorm(-abs(tvalueC))
  pvalueEtot <- 2 * pnorm(-abs(tvalueEtot))
  pvalueCtot <- 2 * pnorm(-abs(tvalueCtot))
  total <- sum(coefE) + sum(coefC)
  pctE  <- coefE/total * 100
  pctC  <- coefC/total * 100
  pctEtot  <- coefEtot/total * 100
  pctCtot  <- coefCtot/total * 100
  coeftot  <- rbind(coefEtot, coefCtot, total)
  s.errtot <- rbind(s.errEtot, s.errCtot, s.errEtot + s.errCtot)
  tvaluetot <- rbind(tvalueEtot, tvalueCtot, 
                     total/(s.errEtot + s.errCtot))
  pvaluetot <- rbind(pvalueEtot, pvalueCtot, 
                     2 * pnorm(-abs(total/(s.errEtot + s.errCtot))))
  pcttot    <- rbind(pctEtot, pctCtot, pctEtot+pctCtot)
  rownames(coeftot) <- c("E", "C", "E + C")
  cnTot     <- paste0(rownames(coeftot))
  coef.tableE <- cbind(coefE, s.errE, tvalueE, pvalueE, pctE)
  dimnames(coef.tableE) <- list(cnE, c(dn, "z value", "Pr(>|z|)", "pct"))
  coef.tableC <- cbind(coefC, s.errC, tvalueC, pvalueC, pctC)
  dimnames(coef.tableC) <- list(cnC, c(dn, "z value", "Pr(>|z|)", "pct"))
  coef.tabletot <- cbind(coeftot, s.errtot, tvaluetot, pvaluetot, pcttot)
  dimnames(coef.tabletot) <- list(cnTot, c(dn, "z value", "Pr(>|z|)", "pct"))
  keep <- match(c("call", "terms"), names(object), 0L)
  # ans <- c(object[keep], list(coefE=coefE, coefC=coefC,
  #                          summaryE=coef.tableE,   
  #                          summaryC=coef.tableC, 
  #       covmatE=covmatE, covmatC=covmatC ))
  ans <- c(object[keep], list(summaryE=coef.tableE,   
                              summaryC=coef.tableC, 
                              summaryTotal = coef.tabletot))
  class(ans) <- "summary.decomp"
  #  return(as.data.frame(rbind(coef.tableE, coef.tableC)))
  return(ans)
}

