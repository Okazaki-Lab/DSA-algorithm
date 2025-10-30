"""
Example Usage of DSA Algorithm
Python Implementation
"""

import numpy as np
import pandas as pd
from dsa_algorithm import DSAAlgorithm

# Example 1: Normal distribution (blood pressure data)
print("=" * 60)
print("Example 1: Normal Distribution (Blood Pressure)")
print("=" * 60)

# Load data
df_normal = pd.read_csv('data_normal.csv')
data_normal = df_normal['blood_pressure'].values

# Run DSA algorithm
dsa1 = DSAAlgorithm(data_normal, estimand="Mean blood pressure difference between treatment groups")
result1 = dsa1.run()

# Generate audit report
dsa1.generate_audit_report('audit_report_example1.txt')

print("\n")

# Example 2: Log-normal distribution (length of stay)
print("=" * 60)
print("Example 2: Log-Normal Distribution (Length of Stay)")
print("=" * 60)

# Load data
df_lognormal = pd.read_csv('data_lognormal.csv')
data_lognormal = df_lognormal['length_of_stay_days'].values

# Run DSA algorithm
dsa2 = DSAAlgorithm(data_lognormal, estimand="Median length of stay difference")
result2 = dsa2.run()

# Generate audit report
dsa2.generate_audit_report('audit_report_example2.txt')

print("\n")

# Example 3: Poisson distribution (adverse event counts)
print("=" * 60)
print("Example 3: Poisson Distribution (Adverse Event Counts)")
print("=" * 60)

# Load data
df_poisson = pd.read_csv('data_poisson.csv')
data_poisson = df_poisson['adverse_event_count'].values

# Run DSA algorithm
dsa3 = DSAAlgorithm(data_poisson, estimand="Rate ratio of adverse events")
result3 = dsa3.run(distributions=['poisson', 'negbinom'])

# Generate audit report
dsa3.generate_audit_report('audit_report_example3.txt')

print("\n")

# Example 4: Mixture distribution (heterogeneous patient populations)
print("=" * 60)
print("Example 4: Mixture Distribution (Biomarker Levels)")
print("=" * 60)

# Load data
df_mixture = pd.read_csv('data_mixture.csv')
data_mixture = df_mixture['biomarker_level'].values

# Run DSA algorithm
dsa4 = DSAAlgorithm(data_mixture, estimand="Mean biomarker difference accounting for population heterogeneity")
result4 = dsa4.run()

# Generate audit report
dsa4.generate_audit_report('audit_report_example4.txt')

print("\n")
print("=" * 60)
print("All examples completed successfully!")
print("Audit reports saved to:")
print("  - audit_report_example1.txt")
print("  - audit_report_example2.txt")
print("  - audit_report_example3.txt")
print("  - audit_report_example4.txt")
print("=" * 60)
