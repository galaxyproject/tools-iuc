glmnet wrappers
===============

This is a self installing Galaxy tool exposing the glmnet_ R package which has excellent documentation at
glmnet_ Minimal details are provided in this wrapper - please RTM to get the best out of it.

The tool exposes the entire range of penalised maximum likelihood
GLM models ranging from pure lasso (set alpha to 1) to pure ridge-regression (set alpha to 0). 

These models can be k-fold internally cross validated to help select an "optimal" predictive or classification
algorithm. Predictive coefficients for each included independent variable are output for each model. 

Predictors can be forced into models to adjust for known confounders or explanatory factors.

The glmnet_ implementation of the coordinate descent algorithm is fast and efficient even on relatively large problems
with tens of thousands of predictors and thousands of samples - such as normalised microarray intensities and anthropometry
on a very large sample of obese patients. 

The user supplies a tabular file with rows as samples and columns containing observations, then chooses 
as many predictors as required. A separate model will be output for each of potentially multiple dependent
variables. Models are reported as the coefficients for terms in an 'optimal' model.
These optimal predictors are selected by repeatedly setting
aside a random subsample, building a model in the remainder and estimating AUC or deviance 
using  k (default 10) fold internal cross validation. For each of these steps, a random 1/k 
of the samples are set aside and used to estiamte performance of an optimal model estimated 
from the remaining samples. Plots are provided showing the range of these (eg 10) internal validation 
estimates and mean model AUC (binomial) or residual deviance plots at each penalty increment step.

A full range of link functions are available including Gaussian, Poisson, Binomial and
Cox proportional hazard time to failure for censored data in this wrapper.

Note that multinomial and multiresponse gaussian models are NOT yet implemented since I have not yet
had use for them - send code!

.. _glmnet: http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html

Wrapper author: Ross Lazarus
19 october 2014

