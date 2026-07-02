import numpy as np
import pandas as pd
import itertools
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.utils import resample
import warnings

warnings.filterwarnings("ignore")

class ICBTAnalyzer:
    """
    ICBT (Iterative Combination of Biomarkers and Thresholds) analyzer engine.
    Implements the methodology described by Robin et al. (2013) for 
    optimal biomarker panel classification.
    """
    def __init__(self, data, target_col):
        self.df = data.dropna(subset=[target_col])
        self.target = target_col
        self.results = None

    def _get_local_extrema(self, bio, y_true):
        vals = self.df[bio].values
        # Detect biological directionality
        is_inverse = np.median(vals[y_true == 0]) > np.median(vals[y_true == 1])
        vals_roc = -vals if is_inverse else vals
        
        fpr, tpr, thresh = roc_curve(y_true, vals_roc)
        sp, se = 1 - fpr, tpr
        
        # Calculate local extrema for threshold optimization
        extrems = [thresh[j] for j in range(1, len(thresh)-1) 
                   if sp[j] >= sp[j-1] and se[j] >= se[j+1]]
        
        # Add Youden Index maximum as a fallback
        youden = (se + sp - 1)
        best_idx = np.flatnonzero(youden == youden.max())[0]
        if thresh[best_idx] not in extrems:
            extrems.append(thresh[best_idx])
            
        return [(-t if is_inverse else t) for t in extrems if not np.isinf(t)], is_inverse

    def run(self, biomarkers, n_bootstrap=1000):
        y_true = self.df[self.target].values
        rules = {bio: self._get_local_extrema(bio, y_true) for bio in biomarkers}
        
        best_auc, best_model, best_scores = 0, None, None
        
        # Exhaustive search for optimal panel
        for r in range(1, len(biomarkers) + 1):
            for comb in itertools.combinations(biomarkers, r):
                threshold_list = [rules[b][0] for b in comb]
                
                for thresh_comb in itertools.product(*threshold_list):
                    score = np.zeros(len(y_true))
                    for i, bio in enumerate(comb):
                        is_inv = rules[bio][1]
                        score += (self.df[bio] <= thresh_comb[i]).astype(int) if is_inv \
                                 else (self.df[bio] >= thresh_comb[i]).astype(int)
                    
                    if len(np.unique(score)) > 1:
                        auc = roc_auc_score(y_true, score)
                        if auc > best_auc:
                            best_auc, best_scores = auc, score
                            # Save model configuration
                            fpr_p, tpr_p, thr_p = roc_curve(y_true, score)
                            idx_ts = np.argmax(tpr_p - fpr_p)
                            best_model = {
                                'panell': comb,
                                'llindars': thresh_comb,
                                'signes': ["<=" if rules[b][1] else ">=" for b in comb],
                                'Ts': thr_p[idx_ts],
                                'Sens': tpr_p[idx_ts] * 100,
                                'Spec': (1 - fpr_p[idx_ts]) * 100
                            }
        
        self.results = {'model': best_model, 'auc': best_auc, 'scores': best_scores}
        return self.results
