# HAB_EDM
Analysis related to the manuscript: A Prototype Coupled Modelling Approach for Predicting Harmful Algal Blooms: A Case Study in Chile
# HAB forecasting using Empirical Dynamic Modeling

This repository contains the Python code and data used in the paper:

"Analysis related to the manuscript: A Prototype Coupled Modelling Approach for Predicting Harmful Algal Blooms: A Case Study in Chile"

## Files
- `pyEDM_RedTide.py`: pyEDM analysis (CCM, S-map, surrogate tests)
- `pyEDM_Forecast.py`: pyEDM analysis (Multivariate S-map for prediction)
- `BD_FITOSenoReloncaviMetriReady.csv`: phytoplankton time-series data

## Requirements
Python >= 1.13.1
pyEDM, pandas, numpy, matplotlib, scipy, seaborn

## Usage
```bash
python pyEDM_CCM-SMap.py
python pyEDM_Forecast.py
