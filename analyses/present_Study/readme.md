# Overview

Here we analyze the primary dataset (Study 1) and some data from a different project (Study 2).

R-scripts: 

For study 1 we analyze behavior (Table S1), generate Figure 2B panels and analyze the average response locked EEG activity at Pz -700 - -200 ms as a function of RT and value difference.


For study 2 we analyze the average response locked EEG activity at Pz -700 - -200 ms as a function of RT and value difference and replicate the findings in study 1 (files are includeded in this folder).


For EEG analyses:

We preprocess the data (Study 1 and 2) and then submit to 1 of 2 analysis pipelines (Study 1 only):

## 1) MASS-univariate analses in Matlab

1) Preprocessing
2) Epoching
3) massStats
4) getResMass
5) export (for component analyses in R)

## 2) MASS-univariate and deconvolution analyses in Matlab using unfold for matlab

takes preprocessed files:
1) unfold deconvolution: BASB_unfold.m (Study 1), BASS_unfold.m (Study 2)
2) Permutation tests: Permutation_test_on_ufresults.m (Study 1), Permutation_test_on_ufresults_BASS.m (Study 2)
3) Extract results: GetResMASS_uf.m (Study1), GetResMASS_uf_BASS.m (Study 2)
