---
layout: default
title: Partial Least Squares regression modeling
description: Matlab Toolbox for the creation of a PLSR model
---

# Matlab Toolbox for the creation of a PLSR model

Please cite :

Jacq, K., Perrette, Y., et al. (2019) High-resolution prediction of organic matter concentration with hyperspectral imaging on a sediment core. Science of the Total Environment 663: 236–244

or

Jacq, K., Giguet-Covex, C., et al. (2019) High-resolution grain size distribution of sediment core with hyperspectral imaging. Sedimentary Geology 393–394:

# Worflow

<img src="CreateModel.jpg" width="600"/>

## Raw data

Two types of data are mandatory to use this toolbox. First, an hyperspectral image (M) and some information, depth (dm) and wavelength (wl). Secondly, the variable(s) to predict (Y), that can be a vector or a matrix, its sampling depth vector (dy) and the label (Yn).

## Reflectance or pseudo-absorbance
Generally, the hyperspectral image unity is reflectance (R). According to Beer-Lambert law, the chemical concentration of a compound is related to the absorbance (A). Thus, the hyperspectral image can be converted in pseudo-absorbance with the formula: A=log(1/R).

## Depth calibration
The two depth vectors are then compared to find the hyperspectral pixels related to a sampling area.

## Spectral preprocessing
![](Preprocessing.jpg)

In this toolbox, four main preprocessing are used, also nine combinations of them. Thus, thirteen preprocessing and the raw data are compared to find an optimal one. They allow to reduce noise and highlights relevant spectral information to predict the interest variable(s).

Detending is used to remove the baseline. Multiplicative Scatter Correction (MSC) and Standard Normal Variate (SNV) are used to correct the spectra from light scaterring. The Savitzky–Golay filter is used to derivate the spectra and to reduce additive effects (baseline offset and slope).

```markdown
Vidal, M., Amigo, J.M. (2012) Pre-processing of hyperspectral images. Essential steps before image analysis. Chemometrics and Intelligent Laboratory Systems 117: 138–148

Rinnan, Å., Berg, F. van den, Engelsen, S.B. (2009) Review of the most common pre-processing techniques for near-infrared spectra. TrAC Trends in Analytical Chemistry 28: 1201–1222

Barnes, R.J., Dhanoa, M.S., Lister, S.J. (1989) Standard Normal Variate Transformation and De-Trending of Near-Infrared Diffuse Reflectance Spectra. Applied Spectroscopy 43: 772–777

Savitzky, A., Golay, M.J.E. (1964) Smoothing and Differentiation of Data by Simplified Least Squares Procedures. Analytical Chemistry 36: 1627–1639
```

## Normalization


## Subsampling

## Optimal preprocessing

## Complete model estimation

## Variable selection

## Reduce model estimation
