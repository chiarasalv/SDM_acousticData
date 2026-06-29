# SDM_acousticData

Code and derived data associated with the manuscript:

“...and the forests will echo with laughter!”:
Remote Sensing and Bioacoustic Integration for Assessing Picus viridis Distribution in Structurally Heterogeneous Riparian Forests

## Authors
Chiara Salvatori et al.

## Overview

This repository contains the R scripts and derived datasets used to model the acoustic habitat suitability of Picus viridis in riparian forests of South Tyrol (Italy) using:

- AudioMoth passive acoustic monitoring BirdNET detections
- LiDAR-derived structural metrics
- Sentinel-2 NDVI descriptors
- Partial Least Squares Regression (PLSR)


## Data availability
Processed datasets are available on Zenodo: DOI: XXXXX

Derived modelling datasets required to reproduce the analyses are provided where permitted.

## Requirements

Main R packages used:

- terra
- caret
- pROC
- pls
- ggplot2
- dplyr
- ape

## Reproducibility

Analyses were performed using R version 4.5.2
