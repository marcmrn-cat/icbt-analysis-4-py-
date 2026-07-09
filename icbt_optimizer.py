"""
ICBT (Iterative Combination of Biomarkers and Thresholds) Optimizer
-------------------------------------------------------------------
This module implements a highly optimized, vectorized engine to find the 
best combinations and cut-off values for biomarker panels.

Features:
- Fast vectorization using NumPy.
- Configurable optimization strategies: 'balanced' (Youden), 'max_spec', 'max_sens'.
- Top-K threshold pre-selection to prevent combinatorial explosion.
"""

import numpy as np
import pandas as pd
import itertools
from sklearn.metrics import roc_curve, roc_auc_score
import warnings
from typing import List, Dict, Tuple, Any

warnings.filterwarnings("ignore")

class ICBTOptimizer:
    """
    Engine to find optimal biomarker thresholds for clinical decision rules.
    """
    def __init__(self, data: pd.DataFrame, target_col: str):
        self.df = data.dropna(subset=[target_col]).copy()
        self.target = target_col
        self.y_true = self.df[self.target].values

    def _get_fast_thresholds(self, bio: str, max_thresh: int = 15) -> Tuple[List[float], bool]:
        """
        Pre-selects the best cut-offs for a single biomarker based on the Youden Index.
        """
        vals = self.df[bio].values
        # Detect biological directionality
        is_inverse = np.median(vals[self.y_true == 0]) > np.median(vals[self.y_true == 1])
        vals_roc = -vals if is_inverse else vals
        
        fpr, tpr, thresh = roc_curve(self.y_true, vals_roc)
        sp, se = 1 - fpr, tpr
        youden = se + sp - 1
        
        # Keep only the top 'max_thresh' thresholds to optimize speed
        best_indices = np.argsort(youden)[::-1][:max_thresh]
        best_thresh = thresh[best_indices]
        
        final_thresh = [(-t if is_inverse else t) for t in best_thresh if not np.isinf(t)]
        return list(set(final_thresh)), is_inverse

    def run(self, biomarkers: List[str], strategy: str = 'balanced', min_req: float = 0.40) -> Dict[str, Any]:
        """
        Runs the iterative combination across all provided biomarkers.
        
        Parameters:
        -----------
        biomarkers : List[str]
            List of column names corresponding to the biomarkers.
        strategy : str
            'balanced' (max AUC), 'max_spec' (high specificity), 'max_sens' (high sensitivity).
        min_req : float
            Minimum required sensitivity/specificity when forcing max_spec or max_sens.
        """
        rules = {bio: self._get_fast_thresholds(bio) for bio in biomarkers}
        X_dict = {bio: self.df[bio].values for bio in biomarkers}
        
        best_metric_global = -1
        best_auc, best_model, best_scores = 0, None, None
        
        for r in range(1, len(biomarkers) + 1):
            for comb in itertools.combinations(biomarkers, r):
                threshold_list = [rules[b][0] for b in comb]
                is_inv_list = [rules[b][1] for b in comb]
                
                # Fast vectorized evaluation
                for thresh_comb in itertools.product(*threshold_list):
                    score = np.zeros(len(self.y_true))
                    
                    for i, bio in enumerate(comb):
                        val = X_dict[bio]
                        t = thresh_comb[i]
                        if is_inv_list[i]:
                            score += (val <= t)
                        else:
                            score += (val >= t)
                    
                    if len(np.unique(score)) > 1:
                        fpr_p, tpr_p, thr_p = roc_curve(self.y_true, score)
                        sens = tpr_p
                        spec = 1 - fpr_p
                        auc = roc_auc_score(self.y_true, score)
                        
                        # Apply clinical optimization strategy
                        if strategy == 'max_spec':
                            valid_idx = np.where(sens >= min_req)[0]
                            idx = valid_idx[np.argmax(spec[valid_idx])] if len(valid_idx) > 0 else np.argmax(spec)
                            metric_score = spec[idx] * 0.8 + sens[idx] * 0.2
                            
                        elif strategy == 'max_sens':
                            valid_idx = np.where(spec >= min_req)[0]
                            idx = valid_idx[np.argmax(sens[valid_idx])] if len(valid_idx) > 0 else np.argmax(sens)
                            metric_score = sens[idx] * 0.8 + spec[idx] * 0.2
                            
                        else:
                            idx = np.argmax(sens + spec - 1)
                            metric_score = auc
                            
                        # Update best model
                        if metric_score > best_metric_global:
                            best_metric_global = metric_score
                            best_auc = auc
                            best_scores = score.copy()
                            best_model = {
                                'panel': comb,
                                'thresholds': thresh_comb,
                                'signs': ["<=" if is_inv_list[i] else ">=" for i in range(len(comb))],
                                'Ts': thr_p[idx],
                                'Sens': sens[idx] * 100,
                                'Spec': spec[idx] * 100
                            }
                            
        return {'model': best_model, 'auc': best_auc, 'scores': best_scores}


# EXECUTION PIPELINE (Example usage for the study)

if __name__ == "__main__":
    # This block only runs if you execute this file directly
    # (Assumes 'df' is previously loaded with your clinical data)
    
    # --- 1. T1 EVALUATION (Admission) -> Objective: Max Specificity ---
    vars_t1 = ['x', 'x', 'x', 'x']
    df_clean_t1 = df.dropna(subset=vars_t1 + ['Pneumonia_CDC_YN']).copy()
    
    analyzer_t1 = ICBTOptimizer(df_clean_t1, target_col='Pneumonia_CDC_YN')
    print("--- PROCESSING T1 (Strategy: MAX SPECIFICITY) ---")
    res_t1 = analyzer_t1.run(biomarkers=vars_t1, strategy='max_spec')
    mod_t1 = res_t1['model']
    
    print(f"Sensitivity: {mod_t1['Sens']:.1f}%")
    print(f"Specificity: {mod_t1['Spec']:.1f}%")
    print(f"T1 Rule: Positive if Score >= {int(mod_t1['Ts'])}")
    for b, sig, ll in zip(mod_t1['panel'], mod_t1['signs'], mod_t1['thresholds']):
        print(f"  [+1] if {b:<15} {sig} {ll:>6.2f}")
    
    # Save the score back to the main DataFrame
    df.loc[df_clean_t1.index, 'Score_Blood_T1_MaxSpec'] = res_t1['scores']

    # --- 2. T2 EVALUATION (24-48h) -> Objective: Max Sensitivity ---
    vars_t2 = ['PGRP_S_T2_ng_mL', 'IL6_T2', 'MRproADM_T2', 'CRP_T2']
    df_clean_t2 = df.dropna(subset=vars_t2 + ['Pneumonia_CDC_YN']).copy()
    
    analyzer_t2 = ICBTOptimizer(df_clean_t2, target_col='Pneumonia_CDC_YN')
    print("\n--- PROCESSING T2 (Strategy: MAX SENSITIVITY) ---")
    res_t2 = analyzer_t2.run(biomarkers=vars_t2, strategy='max_sens')
    mod_t2 = res_t2['model']
    
    print(f"Sensitivity: {mod_t2['Sens']:.1f}%")
    print(f"Specificity: {mod_t2['Spec']:.1f}%")
    print(f"T2 Rule: Positive if Score >= {int(mod_t2['Ts'])}")
    for b, sig, ll in zip(mod_t2['panel'], mod_t2['signs'], mod_t2['thresholds']):
        print(f"  [+1] if {b:<15} {sig} {ll:>6.2f}")
    
    # Save the score back to the main DataFrame
    df.loc[df_clean_t2.index, 'Score_Blood_T2_MaxSens'] = res_t2['scores']
