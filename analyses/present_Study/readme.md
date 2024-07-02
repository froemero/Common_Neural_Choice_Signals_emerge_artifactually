# Overview

Here we analyze the primary dataset (Study 1) and some data from a different project (Study 2).

R-scripts: 

For study 1 we analyze behavior (Table S1), generate Figure 2B panels and analyze the average response locked EEG activity at Pz -700 - -200 ms as a function of RT and value difference.


For study 2 we analyze the average response locked EEG activity at Pz -700 - -200 ms as a function of RT and value difference and replicate the findings in study 1.


For EEG analyses:

We preprocess the data (Study 1 and 2) and then submit to 1 of 2 analysis pipelines (Study 1 only):

## 1) MASS-univariate analses in Matlab

1) Preprocessing
2) Epoching
3) massStats
4) getResMass
5) Export (for component analyses in R)

## 2) MASS-univariate and deconvolution analyses in Matlab using unfold for matlab


