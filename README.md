# ICBT-Biomarker-Panel-Optimizer

[cite_start]This repository provides a professional-grade implementation of the **ICBT (Iterative Combination of Biomarkers and Thresholds)** algorithm, based on the methodology established by Robin et al. in *Translational Proteomics (2013)*[cite: 1].

## Objective
The tool creates a clinical "white box" classification model by combining multiple continuous biomarkers into a single panel. [cite_start]It generates an interpretable scoring system, which often outperforms the predictive power of individual biomarkers by identifying optimal threshold combinations[cite: 27, 28].

## How it works: Technical Breakdown
The algorithm follows a rigorous, multi-stage process to ensure clinical validity and statistical robustness:

1.  **Biological Directionality Analysis**: 
    The model automatically detects whether a biomarker is a risk factor or a protective factor by comparing the medians between clinical groups. It adjusts the mathematical directionality to ensure all markers are aligned with "higher risk"[cite: 97].

2.  **Threshold Optimization (Local Extrema)**: 
    Instead of testing every possible value, the algorithm identifies "local extrema" on the ROC curve. These are the specific biological points where sensitivity and specificity are optimally balanced, as defined by Robin et al.[cite: 103, 104].

3.  **Exhaustive Combinatorial Search (ICBT)**: 
    The algorithm utilizes a vectorized search engine to evaluate thousands of combinations of markers and thresholds. It calculates a patient score ($S_p$), which represents the count of biomarkers exceeding their respective thresholds[cite: 93, 94].

4.  **Statistical Robustness (Bootstrapping)**: 
    To prevent overfitting—a common issue in biomarker research—the tool performs 1,000-fold bootstrapping to calculate the 95% Confidence Interval (CI) of the Area Under the Curve (AUC)[cite: 63, 64].

5.  **Clinical Decision Rule ($T_s$)**: 
    The algorithm defines a final cutoff ($T_s$) for the panel. [cite_start]This transforms complex continuous protein levels into a simple, actionable clinical rule (e.g., *"If patient has >= 2 markers elevated, classify as high risk"*), which is essential for medical decision-making at the point of care[cite: 98, 99].

## Features
- **Automatic directionality detection**: Ensures all biomarkers are mathematically aligned with the outcome.
- **Robustness metrics**: Returns AUC values accompanied by 95% Confidence Intervals via bootstrapping.
- **Interpretable output**: Generates a clear, additive scoring rule for clinical use.

## Usage
```python
from icbt_analysis import ICBTAnalyzer

# Initialize with your dataframe and the target variable
analyzer = ICBTAnalyzer(df, target_col='Pneumonia_CDC_YN')

# Execute the ICBT search
results = analyzer.run(biomarkers=['CRP_T1', 'SAA_T1', 'MRproADM_T1', 'IL6_T1'])

print(f"AUC achieved: {results['auc']}")

**Citation**
Robin, X. et al. "Panelomix: A threshold-based algorithm to create panels of biomarkers." Translational Proteomics, 2013. 
