"""
Distribution Structure Analysis (DSA) Algorithm
Python Implementation
License: MIT
Author: Michio Okazaki
Contact: senryaku@si-lab.work
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
import warnings
from datetime import datetime
from typing import Dict, List, Tuple, Optional

warnings.filterwarnings('ignore')


class DSAAlgorithm:
    """
    Distribution Structure Analysis (DSA) Algorithm
    
    A comprehensive framework for rigorous and transparent statistical analysis
    in medical research with automated audit trail generation.
    """
    
    def __init__(self, data: np.ndarray, estimand: str = "Mean difference between groups",
                 alpha: float = 0.05):
        """
        Initialize DSA Algorithm
        
        Parameters:
        -----------
        data : np.ndarray
            Numeric array of observations
        estimand : str
            Description of the estimand
        alpha : float
            Significance level for statistical tests (default: 0.05)
        """
        self.data = np.array(data)
        self.estimand = estimand
        self.alpha = alpha
        self.audit_log = {
            'timestamp': datetime.now(),
            'estimand': estimand,
            'sample_size': len(data),
            'data_summary': self._get_data_summary(),
            'distributions_tested': []
        }
        
    def _get_data_summary(self) -> Dict:
        """Get summary statistics of the data"""
        return {
            'mean': np.mean(self.data),
            'median': np.median(self.data),
            'std': np.std(self.data),
            'min': np.min(self.data),
            'max': np.max(self.data),
            'q25': np.percentile(self.data, 25),
            'q75': np.percentile(self.data, 75)
        }
    
    def run(self, distributions: Optional[List[str]] = None) -> Dict:
        """
        Run the DSA algorithm
        
        Parameters:
        -----------
        distributions : List[str], optional
            List of distributions to test. If None, tests all available distributions.
            
        Returns:
        --------
        Dict containing selected distribution, parameters, GOF metrics, and audit log
        """
        if distributions is None:
            distributions = ['normal', 'lognormal', 'exponential', 'weibull', 
                           'gamma', 'poisson', 'negbinom']
        
        print("=== DSA Algorithm Started ===")
        print(f"Estimand: {self.estimand}")
        print(f"Sample size: {len(self.data)}\n")
        
        # Step 1: Data validation
        if np.any(np.isnan(self.data)):
            raise ValueError("Data contains missing values. Please handle missing data before analysis.")
        
        # Step 2: Distribution identification
        print("Step 2: Testing candidate distributions...")
        results = self._test_distributions(distributions)
        
        # Step 3: Goodness-of-fit assessment
        print("\nStep 3: Goodness-of-fit assessment...")
        gof_results = self._assess_goodness_of_fit(results)
        
        # Step 4: Select best distribution
        print("\nStep 4: Selecting best distribution...")
        best_dist = self._select_best_distribution(gof_results)
        
        # Step 5: Quality control flags
        qc_flags = self._assign_quality_flags(best_dist)
        
        # Update audit log
        self.audit_log['results'] = results
        self.audit_log['gof_results'] = gof_results
        self.audit_log['selected_distribution'] = best_dist
        self.audit_log['qc_flags'] = qc_flags
        self.audit_log['completion_time'] = datetime.now()
        
        # Print summary
        print("\n=== DSA Algorithm Completed ===")
        print(f"Selected distribution: {best_dist['distribution']}")
        print(f"Quality flag: {qc_flags['flag']}")
        print(f"AIC: {best_dist['aic']:.2f}")
        print(f"BIC: {best_dist['bic']:.2f}\n")
        
        return {
            'distribution': best_dist['distribution'],
            'parameters': best_dist['parameters'],
            'gof_metrics': best_dist['gof_metrics'],
            'qc_flags': qc_flags,
            'audit_log': self.audit_log
        }
    
    def _test_distributions(self, distributions: List[str]) -> Dict:
        """Test multiple candidate distributions"""
        results = {}
        
        for dist in distributions:
            try:
                if dist == 'normal':
                    params = stats.norm.fit(self.data)
                    loglik = np.sum(stats.norm.logpdf(self.data, *params))
                    n_params = 2
                    
                elif dist == 'lognormal':
                    if np.any(self.data <= 0):
                        print(f"  Skipping lognormal (data contains non-positive values)")
                        continue
                    params = stats.lognorm.fit(self.data)
                    loglik = np.sum(stats.lognorm.logpdf(self.data, *params))
                    n_params = 3
                    
                elif dist == 'exponential':
                    if np.any(self.data < 0):
                        print(f"  Skipping exponential (data contains negative values)")
                        continue
                    params = stats.expon.fit(self.data)
                    loglik = np.sum(stats.expon.logpdf(self.data, *params))
                    n_params = 2
                    
                elif dist == 'weibull':
                    if np.any(self.data < 0):
                        print(f"  Skipping Weibull (data contains negative values)")
                        continue
                    params = stats.weibull_min.fit(self.data)
                    loglik = np.sum(stats.weibull_min.logpdf(self.data, *params))
                    n_params = 3
                    
                elif dist == 'gamma':
                    if np.any(self.data <= 0):
                        print(f"  Skipping gamma (data contains non-positive values)")
                        continue
                    params = stats.gamma.fit(self.data)
                    loglik = np.sum(stats.gamma.logpdf(self.data, *params))
                    n_params = 3
                    
                elif dist == 'poisson':
                    if np.any(self.data < 0) or not np.all(self.data == np.floor(self.data)):
                        print(f"  Skipping Poisson (data must be non-negative integers)")
                        continue
                    lambda_param = np.mean(self.data)
                    params = (lambda_param,)
                    loglik = np.sum(stats.poisson.logpmf(self.data.astype(int), lambda_param))
                    n_params = 1
                    
                elif dist == 'negbinom':
                    if np.any(self.data < 0) or not np.all(self.data == np.floor(self.data)):
                        print(f"  Skipping negative binomial (data must be non-negative integers)")
                        continue
                    # Fit negative binomial using method of moments
                    mean = np.mean(self.data)
                    var = np.var(self.data)
                    if var <= mean:
                        print(f"  Skipping negative binomial (variance <= mean)")
                        continue
                    n = mean**2 / (var - mean)
                    p = mean / var
                    params = (n, p)
                    loglik = np.sum(stats.nbinom.logpmf(self.data.astype(int), n, p))
                    n_params = 2
                    
                else:
                    print(f"  Unknown distribution: {dist}")
                    continue
                
                # Calculate AIC and BIC
                n = len(self.data)
                aic = 2 * n_params - 2 * loglik
                bic = n_params * np.log(n) - 2 * loglik
                
                results[dist] = {
                    'params': params,
                    'aic': aic,
                    'bic': bic,
                    'loglik': loglik,
                    'n_params': n_params
                }
                
                print(f"  {dist}: AIC = {aic:.2f}, BIC = {bic:.2f}")
                
            except Exception as e:
                print(f"  Error fitting {dist}: {str(e)}")
        
        self.audit_log['distributions_tested'] = list(results.keys())
        return results
    
    def _assess_goodness_of_fit(self, results: Dict) -> Dict:
        """Assess goodness-of-fit for all candidate distributions"""
        gof_results = {}
        
        for dist_name, result in results.items():
            # Kolmogorov-Smirnov test
            if dist_name == 'normal':
                ks_stat, ks_pvalue = stats.kstest(self.data, 'norm', args=result['params'])
            elif dist_name == 'lognormal':
                ks_stat, ks_pvalue = stats.kstest(self.data, 'lognorm', args=result['params'])
            elif dist_name == 'exponential':
                ks_stat, ks_pvalue = stats.kstest(self.data, 'expon', args=result['params'])
            elif dist_name == 'weibull':
                ks_stat, ks_pvalue = stats.kstest(self.data, 'weibull_min', args=result['params'])
            elif dist_name == 'gamma':
                ks_stat, ks_pvalue = stats.kstest(self.data, 'gamma', args=result['params'])
            else:
                ks_stat, ks_pvalue = np.nan, np.nan
            
            gof_results[dist_name] = {
                'aic': result['aic'],
                'bic': result['bic'],
                'ks_statistic': ks_stat,
                'ks_pvalue': ks_pvalue,
                'params': result['params']
            }
        
        return gof_results
    
    def _select_best_distribution(self, gof_results: Dict) -> Dict:
        """Select best distribution based on hierarchical criteria"""
        # Priority 3: Information criteria (BIC)
        bic_values = {name: result['bic'] for name, result in gof_results.items()}
        best_dist_name = min(bic_values, key=bic_values.get)
        
        return {
            'distribution': best_dist_name,
            'parameters': gof_results[best_dist_name]['params'],
            'aic': gof_results[best_dist_name]['aic'],
            'bic': gof_results[best_dist_name]['bic'],
            'gof_metrics': gof_results[best_dist_name]
        }
    
    def _assign_quality_flags(self, best_dist: Dict) -> Dict:
        """Assign quality control flags (Red/Yellow/Green)"""
        ks_pvalue = best_dist['gof_metrics']['ks_pvalue']
        
        if np.isnan(ks_pvalue):
            flag = "YELLOW"
            message = "K-S test not applicable for discrete distributions"
        elif ks_pvalue < 0.01:
            flag = "RED"
            message = "Poor fit: Consider alternative distributions or transformations"
        elif ks_pvalue < 0.05:
            flag = "YELLOW"
            message = "Marginal fit: Proceed with caution and sensitivity analysis"
        else:
            flag = "GREEN"
            message = "Good fit: Distribution is appropriate"
        
        return {'flag': flag, 'message': message, 'ks_pvalue': ks_pvalue}
    
    def generate_audit_report(self, output_file: str = "dsa_audit_report.txt"):
        """Generate audit report"""
        with open(output_file, 'w') as f:
            f.write("=== DSA Algorithm Audit Report ===\n\n")
            f.write(f"Timestamp: {self.audit_log['timestamp']}\n")
            f.write(f"Estimand: {self.audit_log['estimand']}\n")
            f.write(f"Sample size: {self.audit_log['sample_size']}\n\n")
            
            f.write("Data Summary:\n")
            for key, value in self.audit_log['data_summary'].items():
                f.write(f"  {key}: {value:.4f}\n")
            f.write("\n")
            
            f.write(f"Distributions Tested: {', '.join(self.audit_log['distributions_tested'])}\n\n")
            
            if 'selected_distribution' in self.audit_log:
                f.write(f"Selected Distribution: {self.audit_log['selected_distribution']['distribution']}\n")
                f.write(f"Quality Flag: {self.audit_log['qc_flags']['flag']}\n")
                f.write(f"Quality Message: {self.audit_log['qc_flags']['message']}\n\n")
            
            f.write(f"Completion Time: {self.audit_log['completion_time']}\n")
        
        print(f"Audit report saved to: {output_file}")


# Example usage
if __name__ == "__main__":
    # Generate example data
    np.random.seed(123)
    data = np.random.normal(50, 10, 1000)
    
    # Run DSA algorithm
    dsa = DSAAlgorithm(data, estimand="Mean treatment effect")
    result = dsa.run()
    
    # Generate audit report
    dsa.generate_audit_report()
