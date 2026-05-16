# -*- coding: utf-8 -*-
"""
Created on Thu May 14 22:43:34 2026

Author: Adrien Hélias
"""

from ocean import run_slepian_wavelet_variance as rswv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Import and prepare the time series
#####################################
data = pd.read_csv("lightcurve(681203100005718).csv")  # This is a ZTF light curve used for the example
data = data[data["catflags"] != 32768]
data = data.sort_values(by=['mjd']).reset_index(drop=True)
t = data["mjd"].to_numpy()
mag = data["mag"].to_numpy()
magerr = data["magerr"].to_numpy()


# Ways to run the SWV analysis
###############################
# Example 1: Simple SWV run
results = rswv(t, mag, magerr)

# Example 2: Apply 3 sigma filter to the data
# results = rswv(t, mag, magerr, three_sigma_filter=True)

# Example 3: Bin the light curve into 2-day bins
# results = rswv(t, mag, magerr, bin_light_curve=True, bin_value=2)


# Plot the variance curve
##########################
# If you have the redshift of the quasar, you can use it to look at the rest-frame variance curve:
# restframe_timescales = np.log2((2**results[0])/(1 + data["z"]))
# Otherwise, you will look at the observed-frame variance curve

plt.rcParams.update({'font.size': 22})
plt.figure(figsize=(10, 7))
plt.errorbar(results[0], results[1], yerr=results[2], fmt="o-", color="#FF1A1A", ecolor="#000000",
             markeredgecolor='black', elinewidth=1.5, linewidth=2, capsize=3, markersize=9)
plt.xlabel(r"log$_2$(timescale) [days]",fontsize=20)
plt.ylabel(r"log$_2$(variance) [mag]",fontsize=20)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.grid(linestyle='--', linewidth=1.5, alpha=0.25)
plt.show()


