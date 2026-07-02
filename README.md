# ICBT-Biomarker-Panel-Optimizer

[cite_start]This repository contains a Python implementation of the **ICBT (Iterative Combination of Biomarkers and Thresholds)** algorithm, based on the methodology published by Robin et al. in *Translational Proteomics (2013)*[cite: 1].

## Objective
This tool allows for the combination of multiple continuous biomarkers to create a clinical decision rule ("white box") that classifies patients with optimized precision (AUC), often outperforming individual markers.

## Features
- **Automatic directionality detection**: The model analyzes whether the elevation of a marker is a risk factor or a protective one.
- **Local extrema search**: Applies Youden's logic to identify relevant biological thresholds.
- **Exhaustive combinatorial search**: Evaluates all possible combinations of markers to find the panel with the highest predictive power.
- **Bootstrapping**: Automatic calculation of 95% confidence intervals (CI) to ensure statistical robustness.

## Usage
```python
from icbt_analysis import ICBTAnalyzer

analyzer = ICBTAnalyzer(df, target_col='Pneumonia_CDC_YN')
results = analyzer.run(biomarkers=['CRP_T1', 'SAA_T1', 'MRproADM_T1', 'IL6_T1'])
print(f"AUC achieved: {results['auc']}")
