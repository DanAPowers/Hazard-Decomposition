# Hazard-Decomposition
Multivariate Decomposition of Group Differences in Hazard Rates

Routines to decompose differences in hazard rates estimated from continuous-time and discrete-time hazard models as described by Powers and Yun (2009) MULTIVARIATE DECOMPOSITION FOR HAZARD RATE MODELS, Sociological Methodolgy. 

We develop a regression decomposition technique for hazard rate models, where the difference in observed rates is decomposed into components attributable to group differences in characteristics and group differences in effects. The baseline hazard is specified using a piecewise-constant exponential model, which leads to convenient estimation based on a Poisson regression model fit to person-period, or split-episode data. This specification allows for a flexible representation of the baseline hazard and provides a straightforward way to introduce time-varying covariates and time-varying effects. We provide computational details underlying the method and apply the technique to the decomposition of the Black-White difference in first premarital birth rates into components reflecting characteristics and effect contributions of several predictors, as well as the effect contribution attributable to race differences in the baseline hazard.

Three models are supported: continous-time (piecewise-constant exponential), discrete-time logit, and discrete-time complementary log log. A user-specified formula and group design lists are passed to hazard_decomp_functions_svy.R. Although this requires data to be setup as a complex design, simple designs are also supported.

        df <- read.csv(file='hazdata.csv'

        require(survey)

clusters, no weights 

        dfsvy <- svydesign(ids=~iid, weights=~1, cluster=~famid, nest=TRUE, data=df)

no weights, no clusters 

        dfsvy <- svydesign(ids=~iid, weights=~1, nest=TRUE, data=df)

make subsets based on the grouping variable

        Asub <- subset(dfsvy, race == "Black")

        Bsub <- subset(dfsvy, race == "White")

#source main subroutines

        source('hazard_decomp_functions_svy.R')

Hazdecomp_Example.R shows typical calls (and default values)

        call: decomp_model(formula, Asub, Bsub, scale=1, reverse=FALSE, prinitit=FALSE)

where, Asub and Bsub are the design lists for each group, scale is a rate multiplier, reverse=TRUE if groups are swapped, printit to print output.

models: 

        decomp.pwcexp  (continuous-time model piecewise constant exponential via poisson trick (log-exposure as offset))

        decomp.logit   (discrete-time logit)
        
        decomp.cloglog (discrete-time complementary log-log)
        
For example: two decomposions are performed and averaged (see example)

m1a <- decomp.pwcexp(devnt ~ age + pctsmom + nfamtran + medu + 
                      inc1000 + nosibs + magebir + offset(logexp) - 1,
                    Asub, Bsub, scale=1000, reverse=FALSE)
                    
m1b <- decomp.pwcexp(devnt ~ age + pctsmom + nfamtran + medu + 
                       inc1000 + nosibs + magebir + offset(logexp) - 1,
                     Asub, Bsub, scale=1000, reverse=TRUE)




