# 
# decomp.cthaz <- function(bA, bB, xA, xB, offA, offB, varbA, varbB, 
#                          meanA, meanB, wA, wB, scale=NULL, printit) {

decomp.cthaz <- function(A, B, 
                          meanA, meanB, scale=NULL, printit) {

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
        F <- exp(x%*%b + off)*wt
        return(F)
}
  f <-  function(b,x,off,wt) {
        f <-  exp(x%*%b + off)*wt
        return(f)
}
  
decomp(bA, bB, xA, xB, offA, offB, varbA, varbB, 
       meanA, meanB, wA, wB, printit, scale, 
       eval(F), eval(f) )

}
##################################
# 
# decomp.dtlogit <- function(bA, bB, xA, xB, offA, offB, varbA, varbB,
#                          meanA, meanB, wA, wB, scale=NULL, printit) {
decomp.dtlogit <- function(A,  B,
                          meanA, meanB, scale=NULL, printit) {

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
  
  decomp(bA, bB, xA, xB, offA, offB, varbA, varbB, 
         meanA, meanB, wA, wB, printit, scale, 
         eval(F), eval(f) )
  
}

############# cloglog ############
#m0stuff <- list(b0=outm0$b,x0=outm0$x, off0=outm0$off, v0=outm0$v, w0=outm0$w)
#m1stuff <- list(b1=outm1$b,x1=outm0$x, off1=outm1$off, v1=outm1$v, w1=outm1$w)

# decomp.dtcloglog <- function(bA, xA, offA, varA, wA,  bB, xB, offB, varB, 
#                            meanA, meanB, scale=NULL, printit) 
decomp.dtcloglog <- function(A,  B, 
                              scale=NULL, printit) {
  
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
  
  
  decomp(bA, bB, xA, xB, offA, offB, varbA, varbB, 
         meanA, meanB, wA, wB, printit, scale, 
         eval(F), eval(f) )
  
}

