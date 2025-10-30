# Distribution Structure Analysis (DSA) Algorithm

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![medRxiv](https://img.shields.io/badge/medRxiv-preprint-blue)](https://doi.org/10.1101/2025.10.29.XXXXXXX)
[![R](https://img.shields.io/badge/R-4.3.0+-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/)

## Overview

The **Distribution Structure Analysis (DSA) Algorithm** is a comprehensive framework for rigorous and transparent statistical analysis in medical research. It systematically identifies distributional structures, ensures statistical rigor through explicit estimand specification and goodness-of-fit testing, and maintains complete audit trails for regulatory compliance.

### Key Features

- **Systematic Distribution Identification**: Evaluates multiple candidate distributions (normal, log-normal, exponential, Weibull, gamma, power-law, and finite mixture models) and selects the most appropriate model based on hierarchical assessment criteria
- **Estimand Framework Integration**: Implements the ICH E9(R1) estimand framework, ensuring alignment between statistical analysis and research questions
- **Comprehensive Goodness-of-Fit Assessment**: Uses multiple criteria including theoretical plausibility, visual diagnostics, information criteria (AIC/BIC), and statistical tests
- **Causal Inference Support**: Integrates Directed Acyclic Graphs (DAG) to identify confounders, mediators, and colliders
- **Audit-Ready Framework**: Automatically logs all analytical decisions with complete justification and reproducibility
- **Three-Tier Quality Control**: Red/yellow/green flagging system to guide decision-making
- **Open-Source Implementation**: Available in both R and Python with comprehensive documentation

### Performance

- **95.2% accuracy** in correctly identifying distribution types across 1,000 simulated datasets
- **100% reproducibility** in independent verification studies
- **12% error prevention rate** in validation studies, preventing methodological errors that would have led to incorrect conclusions

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage Examples](#usage-examples)
- [Algorithm Components](#algorithm-components)
- [Validation Studies](#validation-studies)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Installation

### R Installation

#### Prerequisites

- R version 4.3.0 or higher
- Required packages: `fitdistrplus`, `MASS`, `mixtools`, `dagitty`, `ggplot2`, `dplyr`

#### Installation Steps

```r
# Install required packages
install.packages(c("fitdistrplus", "MASS", "mixtools", "dagitty", "ggplot2", "dplyr"))

# Download DSA algorithm
# Option 1: Clone the repository
system("git clone https://github.com/Okazaki-Lab/DSA-algorithm.git")

# Option 2: Download directly
download.file("https://github.com/Okazaki-Lab/DSA-algorithm/raw/main/R/DSA_algorithm.R", 
              destfile = "DSA_algorithm.R")

# Load the DSA algorithm
source("DSA_algorithm.R")
```

### Python Installation

#### Prerequisites

- Python version 3.11 or higher
- Required packages: `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, `scikit-learn`, `statsmodels`

#### Installation Steps

```bash
# Create a virtual environment (recommended)
python3 -m venv dsa_env
source dsa_env/bin/activate  # On Windows: dsa_env\Scripts\activate

# Install required packages
pip install numpy scipy pandas matplotlib seaborn scikit-learn statsmodels

# Clone the repository
git clone https://github.com/Okazaki-Lab/DSA-algorithm.git
cd DSA-algorithm

# Or install directly
pip install git+https://github.com/Okazaki-Lab/DSA-algorithm.git
```

---

## Quick Start

### R Quick Start

```r
# Load the DSA algorithm
source("DSA_algorithm.R")

# Load example data
data <- read.csv("data/simulated_data.csv")

# Run DSA analysis
result <- dsa_analyze(
  data = data$values,
  estimand = "mean_difference",
  research_design = "randomized_controlled_trial",
  endpoint = "continuous",
  candidate_distributions = c("normal", "lognormal", "exponential", "weibull", "powerlaw", "mixture"),
  audit_log = TRUE
)

# View results
print(result$summary)
print(result$selected_distribution)
print(result$goodness_of_fit)

# Generate diagnostic plots
plot(result)

# View audit log
print(result$audit_log)
```

### Python Quick Start

```python
from dsa_algorithm import DSAAnalyzer
import pandas as pd

# Load example data
data = pd.read_csv("data/simulated_data.csv")

# Initialize DSA analyzer
analyzer = DSAAnalyzer(
    estimand="mean_difference",
    research_design="randomized_controlled_trial",
    endpoint="continuous",
    audit_log=True
)

# Run DSA analysis
result = analyzer.analyze(data['values'])

# View results
print(result.summary())
print(result.selected_distribution)
print(result.goodness_of_fit)

# Generate diagnostic plots
result.plot()

# View audit log
print(result.audit_log)
```

---

## Usage Examples

### Example 1: Clinical Trial Adverse Event Analysis

This example demonstrates how to analyze adverse event frequencies in a clinical trial, where data often follow power-law distributions rather than normal distributions.

#### R Implementation

```r
# Load adverse event data
ae_data <- read.csv("data/adverse_events.csv")

# Specify estimand
estimand <- list(
  research_question = "Compare adverse event rates between treatment and control",
  target_population = "Adults with Type 2 Diabetes",
  endpoint = "Adverse event frequency (count data)",
  intercurrent_events = "Treatment discontinuation handled as censoring",
  summary_measure = "Rate ratio"
)

# Run DSA analysis
result <- dsa_analyze(
  data = ae_data$ae_count,
  estimand = "rate_ratio",
  research_design = "randomized_controlled_trial",
  endpoint = "count",
  candidate_distributions = c("poisson", "negative_binomial", "zero_inflated_poisson", "powerlaw"),
  audit_log = TRUE,
  quality_control = TRUE
)

# Check quality control flags
if (result$quality_control$flag == "red") {
  cat("STOP: Critical issues detected\n")
  print(result$quality_control$issues)
} else if (result$quality_control$flag == "yellow") {
  cat("RE-EXAMINE: Potential issues detected\n")
  print(result$quality_control$issues)
} else {
  cat("PROCEED: No critical issues detected\n")
}

# View selected distribution
cat("Selected distribution:", result$selected_distribution, "\n")
cat("AIC:", result$goodness_of_fit$AIC, "\n")
cat("BIC:", result$goodness_of_fit$BIC, "\n")

# Generate diagnostic plots
plot(result, type = "qq")
plot(result, type = "pp")
plot(result, type = "rootogram")

# Export audit log
write.csv(result$audit_log, "audit_log_ae_analysis.csv")
```

#### Python Implementation

```python
import pandas as pd
from dsa_algorithm import DSAAnalyzer

# Load adverse event data
ae_data = pd.read_csv("data/adverse_events.csv")

# Specify estimand
estimand = {
    "research_question": "Compare adverse event rates between treatment and control",
    "target_population": "Adults with Type 2 Diabetes",
    "endpoint": "Adverse event frequency (count data)",
    "intercurrent_events": "Treatment discontinuation handled as censoring",
    "summary_measure": "Rate ratio"
}

# Initialize DSA analyzer
analyzer = DSAAnalyzer(
    estimand="rate_ratio",
    research_design="randomized_controlled_trial",
    endpoint="count",
    candidate_distributions=["poisson", "negative_binomial", "zero_inflated_poisson", "powerlaw"],
    audit_log=True,
    quality_control=True
)

# Run DSA analysis
result = analyzer.analyze(ae_data['ae_count'])

# Check quality control flags
if result.quality_control['flag'] == 'red':
    print("STOP: Critical issues detected")
    print(result.quality_control['issues'])
elif result.quality_control['flag'] == 'yellow':
    print("RE-EXAMINE: Potential issues detected")
    print(result.quality_control['issues'])
else:
    print("PROCEED: No critical issues detected")

# View selected distribution
print(f"Selected distribution: {result.selected_distribution}")
print(f"AIC: {result.goodness_of_fit['AIC']}")
print(f"BIC: {result.goodness_of_fit['BIC']}")

# Generate diagnostic plots
result.plot(plot_type='qq')
result.plot(plot_type='pp')
result.plot(plot_type='rootogram')

# Export audit log
result.audit_log.to_csv("audit_log_ae_analysis.csv")
```

### Example 2: Epidemiological Study with Mixture Models

This example demonstrates how to identify population heterogeneity using mixture models in an epidemiological study.

#### R Implementation

```r
# Load disease incidence data
incidence_data <- read.csv("data/disease_incidence.csv")

# Run DSA analysis with mixture models
result <- dsa_analyze(
  data = incidence_data$incidence_rate,
  estimand = "incidence_rate",
  research_design = "cohort_study",
  endpoint = "continuous",
  candidate_distributions = c("normal", "lognormal", "mixture_normal", "mixture_lognormal"),
  mixture_components = 2:4,  # Test 2, 3, and 4 component mixtures
  audit_log = TRUE,
  quality_control = TRUE
)

# If mixture model is selected, examine components
if (grepl("mixture", result$selected_distribution)) {
  cat("Mixture model detected\n")
  cat("Number of components:", result$mixture_info$n_components, "\n")
  cat("Component means:", result$mixture_info$means, "\n")
  cat("Component proportions:", result$mixture_info$proportions, "\n")
  
  # Plot mixture components
  plot(result, type = "mixture_components")
}

# Generate report
generate_report(result, output_file = "incidence_analysis_report.pdf")
```

#### Python Implementation

```python
import pandas as pd
from dsa_algorithm import DSAAnalyzer

# Load disease incidence data
incidence_data = pd.read_csv("data/disease_incidence.csv")

# Initialize DSA analyzer with mixture models
analyzer = DSAAnalyzer(
    estimand="incidence_rate",
    research_design="cohort_study",
    endpoint="continuous",
    candidate_distributions=["normal", "lognormal", "mixture_normal", "mixture_lognormal"],
    mixture_components=range(2, 5),  # Test 2, 3, and 4 component mixtures
    audit_log=True,
    quality_control=True
)

# Run DSA analysis
result = analyzer.analyze(incidence_data['incidence_rate'])

# If mixture model is selected, examine components
if 'mixture' in result.selected_distribution:
    print("Mixture model detected")
    print(f"Number of components: {result.mixture_info['n_components']}")
    print(f"Component means: {result.mixture_info['means']}")
    print(f"Component proportions: {result.mixture_info['proportions']}")
    
    # Plot mixture components
    result.plot(plot_type='mixture_components')

# Generate report
result.generate_report(output_file="incidence_analysis_report.pdf")
```

### Example 3: Health Economics Cost Analysis

This example demonstrates how to analyze healthcare costs, which often exhibit heavy-tailed distributions.

#### R Implementation

```r
# Load healthcare cost data
cost_data <- read.csv("data/healthcare_costs.csv")

# Run DSA analysis
result <- dsa_analyze(
  data = cost_data$total_cost,
  estimand = "mean_cost",
  research_design = "observational_study",
  endpoint = "continuous_positive",  # Costs are always positive
  candidate_distributions = c("lognormal", "gamma", "weibull", "powerlaw", "pareto"),
  audit_log = TRUE,
  quality_control = TRUE
)

# Examine heavy-tailed characteristics
if (result$selected_distribution %in% c("powerlaw", "pareto")) {
  cat("Heavy-tailed distribution detected\n")
  cat("Tail index:", result$tail_index, "\n")
  cat("Proportion of costs in top 1%:", result$tail_proportion_1pct, "\n")
  cat("Proportion of costs in top 5%:", result$tail_proportion_5pct, "\n")
}

# Generate cost distribution plot
plot(result, type = "cost_distribution", 
     highlight_tail = TRUE, 
     tail_threshold = 0.95)
```

#### Python Implementation

```python
import pandas as pd
from dsa_algorithm import DSAAnalyzer

# Load healthcare cost data
cost_data = pd.read_csv("data/healthcare_costs.csv")

# Initialize DSA analyzer
analyzer = DSAAnalyzer(
    estimand="mean_cost",
    research_design="observational_study",
    endpoint="continuous_positive",  # Costs are always positive
    candidate_distributions=["lognormal", "gamma", "weibull", "powerlaw", "pareto"],
    audit_log=True,
    quality_control=True
)

# Run DSA analysis
result = analyzer.analyze(cost_data['total_cost'])

# Examine heavy-tailed characteristics
if result.selected_distribution in ['powerlaw', 'pareto']:
    print("Heavy-tailed distribution detected")
    print(f"Tail index: {result.tail_index}")
    print(f"Proportion of costs in top 1%: {result.tail_proportion_1pct}")
    print(f"Proportion of costs in top 5%: {result.tail_proportion_5pct}")

# Generate cost distribution plot
result.plot(plot_type='cost_distribution', 
            highlight_tail=True, 
            tail_threshold=0.95)
```

---

## Algorithm Components

The DSA algorithm consists of five integrated components:

### 1. Estimand Specification Module

Implements the ICH E9(R1) estimand framework, requiring explicit specification of:
- Research question
- Target population
- Endpoint (variable or outcome)
- Intercurrent events
- Population-level summary measure

### 2. Distribution Identification Module

Evaluates multiple candidate distributions:
- **Continuous**: Normal, log-normal, exponential, Weibull, gamma, power-law, Pareto
- **Count**: Poisson, negative binomial, zero-inflated Poisson, zero-inflated negative binomial
- **Mixture**: Gaussian mixture models, log-normal mixture models, finite mixture models

### 3. Goodness-of-Fit Assessment Module

Uses hierarchical assessment criteria:
1. **Theoretical plausibility**: Based on estimand and data-generating mechanism
2. **Visual diagnostics**: Q-Q plots, P-P plots, rootograms
3. **Information criteria**: AIC, BIC
4. **Statistical tests**: Kolmogorov-Smirnov, Anderson-Darling

### 4. Causal Inference Support Module

Integrates Directed Acyclic Graphs (DAG) to:
- Identify confounders, mediators, and colliders
- Determine minimal adjustment sets
- Verify consistency between distributional model and causal structure

### 5. Audit Trail and Quality Control Module

Provides:
- **Automated audit logging**: All analytical decisions documented
- **Three-tier quality control**:
  - **Red (Stop)**: Critical issues (e.g., estimand-design mismatch, domain violations)
  - **Yellow (Re-examine)**: Potential issues (e.g., poor goodness-of-fit, unstable mixtures)
  - **Green (Proceed)**: No critical or potential issues
- **Reproducibility**: 100% reproducibility in independent verification studies

---

## Validation Studies

### Simulation Study

- **Datasets**: 1,000 synthetic datasets with known distributions
- **Overall accuracy**: 95.2% (95% CI: 93.7%-96.5%)
- **Sensitivity**: 94.8% (95% CI: 93.2%-96.2%)
- **Specificity**: 98.7% (95% CI: 98.1%-99.2%)

### Applied Studies

1. **Clinical Trial Data**: Detected heavy-tailed distributions in adverse event frequencies, leading to more accurate safety assessments
2. **Epidemiological Data**: Revealed multimodal patterns in disease incidence, enabling targeted public health interventions
3. **Health Economics Data**: Identified heavy-tailed cost distributions, improving resource allocation efficiency

---

## Citation

If you use the DSA algorithm in your research, please cite:

### Preprint

```
Okazaki M. A Generalizable Distribution Structure Analysis Algorithm with Audit-Ready Framework for Medical Research. medRxiv. 2025. doi: 10.1101/2025.10.29.XXXXXXX
```

### BibTeX

```bibtex
@article{okazaki2025dsa,
  title={A Generalizable Distribution Structure Analysis Algorithm with Audit-Ready Framework for Medical Research},
  author={Okazaki, Michio},
  journal={medRxiv},
  year={2025},
  doi={10.1101/2025.10.29.XXXXXXX}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### MIT License Summary

- ✅ Commercial use
- ✅ Modification
- ✅ Distribution
- ✅ Private use
- ❌ Liability
- ❌ Warranty

---

## Contact

**Michio Okazaki**
- Email: senryaku@si-lab.work
- Affiliation: S.I Lab Inc., Tokyo, Japan; Visiting Researcher, Chiba University Hospital
- ORCID: [0009-0006-8667-8737](https://orcid.org/0009-0006-8667-8737)

### Reporting Issues

If you encounter any issues or have suggestions for improvements, please:
1. Check the [Issues](https://github.com/Okazaki-Lab/DSA-algorithm/issues) page to see if the issue has already been reported
2. If not, create a new issue with:
   - Clear description of the problem
   - Steps to reproduce
   - Expected behavior
   - Actual behavior
   - System information (R/Python version, OS, package versions)

### Contributing

We welcome contributions from the research community! To contribute:
1. Fork the repository
2. Create a new branch for your feature or bug fix
3. Make your changes with clear commit messages
4. Add tests for new features
5. Submit a pull request with a clear description of your changes

---

## Acknowledgments

We thank the following for their contributions and support:
- Chiba University Hospital for providing research facilities
- S.I Lab Inc. for supporting this research
- The open-source community for developing the statistical packages used in this project

---

## Frequently Asked Questions (FAQ)

### Q1: What types of data can the DSA algorithm analyze?

**A**: The DSA algorithm can analyze:
- Continuous data (e.g., biomarker measurements, survival times, costs)
- Count data (e.g., adverse event frequencies, disease incidence)
- Positive continuous data (e.g., healthcare costs, drug concentrations)

### Q2: How long does the analysis take?

**A**: Computation time depends on dataset size:
- Up to 10,000 observations: < 5 minutes on standard desktop
- 10,000-100,000 observations: 15-30 minutes
- > 100,000 observations: Use parallel processing (available in both R and Python)

### Q3: Can I use the DSA algorithm for regulatory submissions?

**A**: Yes. The audit-ready framework is designed to meet regulatory requirements for transparency and reproducibility. The automated audit logging system documents all analytical decisions, which facilitates regulatory review. However, we recommend discussing the methodology with regulatory statisticians before submission.

### Q4: What if the algorithm flags my analysis as "red" or "yellow"?

**A**: 
- **Red (Stop)**: Critical issues detected (e.g., estimand-design mismatch). Review the issues and resolve them before proceeding.
- **Yellow (Re-examine)**: Potential issues detected (e.g., poor goodness-of-fit). Carefully review the issues, consider sensitivity analyses, and document your decisions.
- **Green (Proceed)**: No critical or potential issues detected. Proceed with confidence.

### Q5: Can I add custom distributions to the algorithm?

**A**: Yes. The algorithm is designed to be extensible. You can add custom distributions by:
1. Defining the probability density function (PDF)
2. Implementing parameter estimation (e.g., maximum likelihood)
3. Adding goodness-of-fit assessment
4. Registering the distribution with the DSA framework

See the [Custom Distributions Guide](docs/custom_distributions.md) for details.

### Q6: How do I interpret mixture models?

**A**: Mixture models indicate population heterogeneity. Each component represents a subpopulation with distinct distributional characteristics. Examine:
- Number of components
- Component means and variances
- Component proportions
- Clinical or epidemiological interpretation of subpopulations

### Q7: Is the DSA algorithm suitable for small sample sizes?

**A**: The algorithm's accuracy increases with sample size. For small samples (N < 50), interpretation should be cautious:
- Visual diagnostics become less reliable
- Information criteria may not discriminate well between models
- Consider using bootstrap methods to assess uncertainty

### Q8: Can I use the DSA algorithm for time-to-event data?

**A**: Yes. Use distributions appropriate for time-to-event data:
- Exponential
- Weibull
- Log-normal
- Gamma

Ensure that censoring is handled appropriately in your analysis.

---

## Version History

### Version 1.0.0 (2025-10-29)
- Initial release
- Core DSA algorithm implementation in R and Python
- Five integrated components: Estimand specification, distribution identification, goodness-of-fit assessment, causal inference support, audit trail
- Three-tier quality control system
- Comprehensive documentation and examples
- Validation studies demonstrating 95.2% accuracy

---

## Roadmap

### Planned Features

- **Version 1.1**: Multivariate distribution support
- **Version 1.2**: Automated causal discovery integration
- **Version 1.3**: Time-to-event data extensions
- **Version 1.4**: Ordinal and compositional data support
- **Version 2.0**: Web-based interface for non-programmers

---

## References

1. Altman DG, Bland JM. The normal distribution. BMJ. 1995;310(6975):298.
2. Limpert E, Stahel WA, Abbt M. Log-normal distributions across the sciences: keys and clues. BioScience. 2001;51(5):341-352.
3. Clauset A, Shalizi CR, Newman MEJ. Power-law distributions in empirical data. SIAM Review. 2009;51(4):661-703.
4. ICH E9(R1) Expert Working Group. ICH E9(R1) addendum on estimands and sensitivity analysis in clinical trials to the guideline on statistical principles for clinical trials. 2019.
5. Pearl J. Causality: Models, Reasoning, and Inference. 2nd ed. Cambridge University Press; 2009.

---

**Last updated**: 2025-10-29
