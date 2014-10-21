glmnet wrappers
===============

This is a self installing Galaxy tool exposing the glmnet_ R package which has excellent documentation at
glmnet_ Minimal details are provided in this wrapper - please RTM to get the best out of it.

The Galaxy user can select from the entire range of pure and hybrid penalised maximum likelihood regression 
models from pure lasso (set alpha=1) to pure ridge-regression (set alpha=0). 

THe implementation and coordinate descent make ti wickedly fast and efficient even on relatively large problems
with thousands of predictors and thousands of samples.

Predictors can be forced into the final models.

The user supplies a tabular file with rows as samples and columns containing observations, then chooses 
as many predictors as required. Each of potentially multiple dependent (predicted) variables will be 
reported with the coefficients for terms in an 'optimal' model. Any number of dependent (y) variables 
can be selected. A separate model will be run and reported for each one.

The optimal predictors reported are selected based on empirical internal classification performance 
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

