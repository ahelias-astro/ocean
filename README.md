# ocean

![logo](ocean/ocean_logo.png)

### Calculation of Slepian Wavelet Variance on irregularly sampled time series

#### Original author: [Matthew J. Graham](https://sites.astro.caltech.edu/~mjg/), California Institute of Technology, USA

#### Modified, documented and maintained by [Adrien Hélias](https://ahelias-astro.github.io/), Western University, Canada

#### Version 1.0.1 (Last updated: July 11, 2026)

This package allows the user to run Slepian Wavelet Variance analysis on irregularly sampled time series. Specifically, *ocean* calculates the variance of the time series at multiple timescales, to give more insight on the type of variability seen in the data, and the timescales for which the variability is the strongest. The code is originally intended to calculate the variance curves of astronomical time series, but it is general enough to be used on other kinds of time series.

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

If you make use of this package, please cite [Hélias et al. (2026)](https://iopscience.iop.org/article/10.1088/1538-3873/ae8f75):
```
@ARTICLE{2026PASP..138h4507H,
       author = {{H{\'e}lias}, Adrien and {Barmby}, Pauline and {Gallagher}, Sarah C. and {Abbassi}, Shahram and {Graham}, Matthew J.},
        title = "{Classifying Quasar Types Without a Spectrum}",
      journal = {\pasp},
     keywords = {Light curves, Quasars, Active galactic nuclei, Wavelet analysis, Time series analysis, Time domain astronomy, 918, 1319, 16, 1918, 1916, 2109, Astrophysics of Galaxies},
         year = 2026,
        month = aug,
       volume = {138},
       number = {8},
          eid = {084507},
        pages = {084507},
          doi = {10.1088/1538-3873/ae8f75},
archivePrefix = {arXiv},
       eprint = {2608.11916},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2026PASP..138h4507H},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
