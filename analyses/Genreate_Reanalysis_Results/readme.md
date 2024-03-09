# Procesdure

These scripts take in outputs from unfold analyses, perform permutation tests to identify significant clusters, and plot the findings.

3 sets of re-analysis results are generated:

## 1) mass-univariate vs deconvolution with only r-locked activity varying by RT
      this is to see how much of the typical effect remains after controlling for overlap with a constant s-locked component.

## 2) mass-univariate vs deconvolution with r-locked AND s-locked activity varying by RT
      this is to compare results with simulation findings for signatures of evidence accumulation under varying s and r variability regimes

## 3) mass-univariate of regular data vs mass-univariate data with s- and r-components removed
      this is a critical test of whether anything like the typical evidence accumulation signals remain to be found if prior to analyses s- and r-    
      locked activity is removed without accounting for anything related to evidence accumulation, only for s- and r-locked activity and their 
      overlap. (Spoiler, we do not find anything after correcting for overlap.) 
