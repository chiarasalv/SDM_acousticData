# SDM_acousticData

Code and processed datasets associated with the manuscript:

> **"...and the forests will echo with laughter!": Remote Sensing and Bioacoustic Integration for Assessing *Picus viridis* Distribution in Structurally Heterogeneous Riparian Forests**

---

# Authors

Chiara Salvatori, Francesco Ceresa, Michele Torresani, Irene Menegaldo, Vincenzo Saponaro, Luca Da Ros, Enrico Tomelleri

---

# Overview

This repository contains the complete R workflow used to reproduce the analyses presented in the manuscript.

The study integrates:

- AudioMoth passive acoustic monitoring
- BirdNET detections
- LiDAR-derived forest structural metrics
- Sentinel-2 vegetation descriptors (NDVI)
- Partial Least Squares Regression (PLSR)

to model the acoustic occurrence and habitat suitability of the European Green Woodpecker (*Picus viridis*) within structurally heterogeneous riparian forests of South Tyrol (Italy).

The repository reproduces:

- data preparation
- model calibration
- leave-one-logger-out validation
- leave-one-area-out transferability assessment
- Moran's I analysis
- habitat suitability mapping
- environmental novelty analysis
- all figures and tables presented in the manuscript

---

# Repository structure

```
SDM_acousticData/

├── scripts/
│      R scripts used in the analyses
│
├── functions/
│      Custom functions (if applicable)
│
├── data/
│      prepared_resampled90_data/
│
├── outputs/
│      Generated figures, tables and rasters
│
└── README.md
```

---

# Data availability

The processed datasets required to reproduce the analyses are archived on Zenodo:

**DOI**

https://doi.org/10.5281/zenodo.21038768

Direct download

https://zenodo.org/records/21038768/files/zenodo_data.zip?download=1

The archive contains:

- PicusViridis_daily_model_data.csv
- logger_coordinates.csv
- raster predictor stacks for the three study areas
- supplementary documentation

The complete AudioMoth recordings and raw multi-species BirdNET outputs are not redistributed. The processed datasets provided through Zenodo are sufficient to reproduce all analyses presented in the manuscript.

---

# Installation

Clone the repository

```bash
git clone https://github.com/chiarasalv/SDM_acousticData.git
```

Download the Zenodo archive.

Extract the archive into

```
data/prepared_resampled90_data/
```

without changing the folder structure.

---

# Requirements

Analyses were performed using

- R 4.5.2

Main packages

- terra
- caret
- pROC
- pls
- ggplot2
- dplyr
- ape

---

# Running the workflow

After downloading the processed datasets, run

```r
source("scripts/my_personal_code.R")
```

The workflow automatically generates all figures, tables and raster outputs reported in the manuscript.

---

# Citation

If you use this repository, please cite:

1. the accompanying manuscript;
2. this GitHub repository;
3. the Zenodo dataset.
