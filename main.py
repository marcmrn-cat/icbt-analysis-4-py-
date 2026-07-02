from icbt_analysis import ICBTAnalyzer
import numpy as np

# Assuming 'df' is already loaded in your environment
# Initialize the analyzer
analyzer = ICBTAnalyzer(df, target_col='Pneumonia_CDC_YN')

# Run analysis
results = analyzer.run(biomarkers=['CRP_T1', 'SAA_T1', 'MRproADM_T1', 'IL6_T1'])
model = results['model']

print("OPTIMIZED BIOMARKER PANEL RESULTS")
print("=" * 60)
print(f"Global AUC : {results['auc']:.3f}")
print(f"Sensitivity: {model['Sens']:.1f}%")
print(f"Specificity: {model['Spec']:.1f}%\n")

print(f"CLINICAL DECISION RULE (Score Ts >= {int(model['Ts'])}):")
for bio, signe, llindar in zip(model['panell'], model['signes'], model['llindars']):
    print(f"   [+1 pt] if {bio:<12} {signe} {llindar:>8.2f}")

# Update dataframe
df['Score_Panell_Sang'] = results['scores']
print("\nDone. 'Score_Panell_Sang' has been added to the DataFrame.")
