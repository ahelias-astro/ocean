# ocean

![logo](ocean/ocean_logo.png)

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
@ARTICLE{2014MNRAS.439..703G,
       author = {{Graham}, Matthew J. and {Djorgovski}, S.~G. and {Drake}, Andrew J. and {Mahabal}, Ashish A. and {Chang}, Melissa and {Stern}, Daniel and {Donalek}, Ciro and {Glikman}, Eilat},
        title = "{A novel variability-based method for quasar selection: evidence for a rest-frame {\ensuremath{\sim}}54 d characteristic time-scale}",
      journal = {\mnras},
     keywords = {methods: data analysis, techniques: photometric, surveys, quasars: general, Astrophysics - Cosmology and Extragalactic Astrophysics},
         year = 2014,
        month = mar,
       volume = {439},
       number = {1},
        pages = {703-718},
          doi = {10.1093/mnras/stt2499},
archivePrefix = {arXiv},
       eprint = {1401.1785},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2014MNRAS.439..703G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

Hélias et al. (2026) [in review]
```
