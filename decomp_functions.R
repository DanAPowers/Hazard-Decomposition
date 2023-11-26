# C.R (from B)
decomp <- function(bA, bB, xA, xB, offA, offB, varbA, varbB, 
                   meanA, meanB, wA, wB, printit, scale, F, f) {
  
F <- eval(F)
f <- eval(f)

mA <- meanA
mB <- meanB

# mA <- NULL
# for (i in 1:ncol(xA)){
#   mA[i] <- mean(xA[,i])
# }
# 
# mB <- NULL
# for (i in 1:ncol(xB)){
#   mB[i] <- mean(xB[,i])
# }
# 
# if (length(meanA) == ncol(xA)) mA <- meanA
# if (length(meanB) == ncol(xB)) mB <- meanB
# #  

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

#varb.b1  <- varb.bA <- matrix(varbA,length(bA),length(bA))
#varb.b2  <- varb.bB <- matrix(varbB,length(bB),length(bB))

dEdb <- array(NA, c(length(bA), length(bA)))
for (k in 1:length(bA)){
  for (l in 1:length(bA)) {
    dEdb[k,l] <- Wdx[k]*(sum(f00*xA[,l])/sum(exp(offA)*wA) -  
                         sum(f01*xB[,l])/sum(exp(offB)*wB))  + 
      dWx[k,l]*E
  }
}

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
#Var.E.k <- dEdb%*%varb.bA%*%t(dEdb)
Var.E.k <- dEdb%*%varbA%*%t(dEdb)
seWdx <- sqrt(diag(Var.E.k))
##  Variances
#   Coefficients (C.k)
#Var.C.k <- dCdbA%*%varb.bA%*%t(dCdbA) + dCdbB%*%varb.bB%*%t(dCdbB)
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
#var.E0 <- t(dEdb.0)%*%varb.b1%*%dEdb.0
var.E0 <- t(dEdb.0)%*%varbA%*%dEdb.0
s.E0   <- sqrt(var.E0)

#var.C0 <- t(dCdb.0A)%*%varb.b1%*%dCdb.0A + t(dCdb.0B)%*%varb.b2%*%dCdb.0B 
var.C0 <- t(dCdb.0A)%*%varbA%*%dCdb.0A + t(dCdb.0B)%*%varbB%*%dCdb.0B 
s.C0   <- sqrt(var.C0)

R <- E + C
v.names <- names(bA)
if (is.null(v.names)) {
  v.names <- paste("Var", 1:length(bA), sep = "") # if names not available
}

if (printit==TRUE) {
  
  cat("TOTAL...........(E + C)", formatC((E+C)*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Raw Contrib.....(E,  C)", formatC(E*scale, flag=" ", dig=4, wid=10, format="f"), 
      formatC(C*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Std. Errors.....(E,  C)", formatC(s.E0*scale, flag=" ", dig=4, wid=10, format="f"),
      formatC(s.C0*scale, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Z-values........(E,  C)", formatC(E/s.E0, flag=" ", dig=4, wid=10, format="f"),
      formatC(C/s.C0, flag=" ", dig=4, wid=10, format="f"), "\n")
  cat("Pct Contrib.....(E,  C)", formatC(100*E/(E+C), flag=" ", dig=4, wid=10, format="f"),
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
