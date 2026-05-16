# -*- coding: utf-8 -*-
"""
Created on Thu May 14 22:50:47 2026

Author: Adrien Hélias
"""

from setuptools import setup, find_packages

setup(name="ocean",
      version="1.0.0",
      description="Classification of quasars using Slepian Wavelet Variance on light curves",
      author="Adrien Helias",
      packages=find_packages(),
      install_requires=["numpy", "pandas", "irlb", "scipy", "spectrum"])

