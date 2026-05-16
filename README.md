# ocean

### Calculation of Slepian Wavelet Variance on irregularly sampled time series

### Original author: [Matthew J. Graham](https://sites.astro.caltech.edu/~mjg/), California Institute of Technology, USA

### Modified, documented and maintained by [Adrien Hélias](https://ahelias-astro.github.io/), Western University, Canada

#### Version 1.0.0 (Last updated: May 16th 2026)

This package allows the user the run Slepian Wavelet Variance analysis on their irregularly sampled time series. Specifically, *ocean* calculates the variance of the time series at multiple timescales, to give more insight on the type of variability seen in the data, and the timescales for which the variability is the strongest. The code is originally intended to classify quasars based on the variance curves of their astronomical light curves, but it is general enough to be used on other kinds of time series.

![example](ocean/example.png)

See **ocean_tutorial.py** to learn how to use *ocean*.
