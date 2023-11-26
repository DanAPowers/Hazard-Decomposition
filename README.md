# Hazard-Decomposition
Multivariate Decomposition of Group Differences in Hazard Rates

Routines to decompose differences in hazard rates estimated from continuous-time and discrete-time hazard models as described by Powers and Yun (2009) MULTIVARIATE DECOMPOSITION FOR HAZARD RATE MODELS, Sociological Methodolgy. 

We develop a regression decomposition technique for hazard rate models, where the difference in observed rates is decomposed into components attributable to group differences in characteristics and group differences in effects. The baseline hazard is specified using a piecewise-constant exponential model, which leads to convenient estimation based on a Poisson regression model fit to person-period, or split-episode data. This specification allows for a flexible representation of the baseline hazard and provides a straightforward way to introduce time-varying covariates and time-varying effects. We provide computational details underlying the method and apply the technique to the decomposition of the Black-White difference in first premarital birth rates into components reflecting characteristics and effect contributions of several predictors, as well as the effect contribution attributable to race differences in the baseline hazard.

Three models are supported: continous-time (piecewise-constant exponential), discrete-time logit, and discrete-time complementary log log. A user-specified formula and group definition (as well as cluster variable if applicable) are passed to glm, components of the glm objects are processed into lists passed to model_setup.R, whose reults are passed to decomp_functions.R.

The Example.R code contains routines for this. A typical call might consist of performing two decompositions (one from each group's perspective) and averaging the results as a solution to the indexing problem.

For example:
# two decomposions: one from each perspective (to be averaged)
C1 <- decomp.cthaz(m1stuff, m0stuff, printit=TRUE, scale=1000)
C2 <- decomp.cthaz(m0stuff, m1stuff, printit=TRUE, scale=1000)
