# ocean

![logo](ocean/ocean logo.png)

### Calculation of Slepian Wavelet Variance on irregularly sampled time series

#### Original author: [Matthew J. Graham](https://sites.astro.caltech.edu/~mjg/), California Institute of Technology, USA

#### Modified, documented and maintained by [Adrien Hélias](https://ahelias-astro.github.io/), Western University, Canada

#### Version 1.0.0 (Last updated: May 16th 2026)

This package allows the user the run Slepian Wavelet Variance analysis on irregularly sampled time series. Specifically, *ocean* calculates the variance of the time series at multiple timescales, to give more insight on the type of variability seen in the data, and the timescales for which the variability is the strongest. The code is originally intended to calculate the variance curves of astronomical time series, but it is general enough to be used on other kinds of time series.

![example](ocean/example.png)

See **ocean_tutorial.py** to learn how to use *ocean*.

### Installation:

Option 1: Clone the repository, place the files and the *ocean* folder containing **ocean_functions.py** in your working folder, and then run:
```
pip install -r requirements.txt
```

Option 2: Open your terminal and run the following command to install directly from Github:
```
pip install git+https://github.com/ahelias-astro/ocean.git
```

### Citation:

If you make use of this package, please cite
```
Hélias et al. (2026) [in review]
```
