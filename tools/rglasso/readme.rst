glmnet wrappers
===============

This is a self installing galaxy tool wrapping the glmnet_ R package from Trevor Hastie and others.

It exposes the entire range of pure and hybrid penalised maximum likelihood regression models from 
pure lasso (alpha=1) to pure ridge-regression (alpha=0). 

It uses coordinate descent and so is wickedly fast and efficient even on relatively large problems
with thousands of predictors and thousands of samples.

Predictors can be forced into the final models.

Any number of dependent (y) variables can be selected. A separate model will be run and reported for
each one.

It will select the 'best' set of parsimonious predictors based on (default 10 fold) cross validation using 
either AUC for binomial dependent variables or mean square error for others to select and report
coefficients in the model at an 'optimal' lambda.

A full range of link functions are available including Gaussian, Poisson, Binomial and
Cox proportional hazard time to failure for censored data in this wrapper.

Note that multinomial and multiresponse gaussian models are NOT yet implemented since I have not yet
worked out exactly what they are for - send code!

.. _glmnet: http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html

Wrapper author: Ross Lazarus
19 october 2014

