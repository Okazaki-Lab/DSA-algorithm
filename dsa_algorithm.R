# Distribution Structure Analysis (DSA) Algorithm
# R Implementation
# License: MIT
# Author: Michio Okazaki
# Contact: senryaku@si-lab.work

# Required packages
required_packages <- c("fitdistrplus", "MASS", "mixtools", "ggplot2", "dplyr")

# Install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

install_if_missing(required_packages)

#' DSA Algorithm: Main Function
#' 
#' @param data Numeric vector of observations
#' @param estimand Character string describing the estimand
#' @param distributions Vector of distribution names to test (default: all)
#' @param alpha Significance level for statistical tests (default: 0.05)
#' @return List containing selected distribution, parameters, GOF metrics, and audit log
#' @export
dsa_algorithm <- function(data, 
                          estimand = "Mean difference between groups",
                          distributions = c("normal", "lognormal", "exponential", 
                                           "weibull", "gamma", "poisson", "negbin"),
                          alpha = 0.05) {
  
  # Initialize audit log
  audit_log <- list(
    timestamp = Sys.time(),
    estimand = estimand,
    sample_size = length(data),
    data_summary = summary(data),
    distributions_tested = distributions
  )
  
  cat("=== DSA Algorithm Started ===\n")
  cat("Estimand:", estimand, "\n")
  cat("Sample size:", length(data), "\n\n")
  
  # Step 1: Data validation
  if (any(is.na(data))) {
    stop("Data contains missing values. Please handle missing data before analysis.")
  }
  
  # Step 2: Distribution identification
  cat("Step 2: Testing candidate distributions...\n")
  results <- test_distributions(data, distributions, alpha)
  
  # Step 3: Goodness-of-fit assessment
  cat("\nStep 3: Goodness-of-fit assessment...\n")
  gof_results <- assess_goodness_of_fit(data, results)
  
  # Step 4: Select best distribution
  cat("\nStep 4: Selecting best distribution...\n")
  best_dist <- select_best_distribution(gof_results)
  
  # Step 5: Quality control flags
  qc_flags <- assign_quality_flags(best_dist, gof_results)
  
  # Update audit log
  audit_log$results <- results
  audit_log$gof_results <- gof_results
  audit_log$selected_distribution <- best_dist
  audit_log$qc_flags <- qc_flags
  audit_log$completion_time <- Sys.time()
  
  # Print summary
  cat("\n=== DSA Algorithm Completed ===\n")
  cat("Selected distribution:", best_dist$distribution, "\n")
  cat("Quality flag:", qc_flags$flag, "\n")
  cat("AIC:", round(best_dist$aic, 2), "\n")
  cat("BIC:", round(best_dist$bic, 2), "\n\n")
  
  return(list(
    distribution = best_dist$distribution,
    parameters = best_dist$parameters,
    gof_metrics = best_dist$gof_metrics,
    qc_flags = qc_flags,
    audit_log = audit_log
  ))
}

#' Test multiple candidate distributions
#' @keywords internal
test_distributions <- function(data, distributions, alpha) {
  results <- list()
  
  for (dist in distributions) {
    tryCatch({
      if (dist == "normal") {
        fit <- fitdist(data, "norm")
      } else if (dist == "lognormal") {
        if (any(data <= 0)) {
          cat("  Skipping lognormal (data contains non-positive values)\n")
          next
        }
        fit <- fitdist(data, "lnorm")
      } else if (dist == "exponential") {
        if (any(data < 0)) {
          cat("  Skipping exponential (data contains negative values)\n")
          next
        }
        fit <- fitdist(data, "exp")
      } else if (dist == "weibull") {
        if (any(data < 0)) {
          cat("  Skipping Weibull (data contains negative values)\n")
          next
        }
        fit <- fitdist(data, "weibull")
      } else if (dist == "gamma") {
        if (any(data <= 0)) {
          cat("  Skipping gamma (data contains non-positive values)\n")
          next
        }
        fit <- fitdist(data, "gamma")
      } else if (dist == "poisson") {
        if (any(data < 0) || any(data != floor(data))) {
          cat("  Skipping Poisson (data must be non-negative integers)\n")
          next
        }
        fit <- fitdist(data, "pois", discrete = TRUE)
      } else if (dist == "negbin") {
        if (any(data < 0) || any(data != floor(data))) {
          cat("  Skipping negative binomial (data must be non-negative integers)\n")
          next
        }
        fit <- fitdist(data, "nbinom", discrete = TRUE)
      } else {
        cat("  Unknown distribution:", dist, "\n")
        next
      }
      
      results[[dist]] <- list(
        fit = fit,
        aic = fit$aic,
        bic = fit$bic,
        loglik = fit$loglik
      )
      
      cat("  ", dist, ": AIC =", round(fit$aic, 2), ", BIC =", round(fit$bic, 2), "\n")
      
    }, error = function(e) {
      cat("  Error fitting", dist, ":", e$message, "\n")
    })
  }
  
  return(results)
}

#' Assess goodness-of-fit for all candidate distributions
#' @keywords internal
assess_goodness_of_fit <- function(data, results) {
  gof_results <- list()
  
  for (dist_name in names(results)) {
    fit <- results[[dist_name]]$fit
    
    # Kolmogorov-Smirnov test
    ks_test <- gofstat(fit, fitnames = dist_name)
    
    gof_results[[dist_name]] <- list(
      aic = results[[dist_name]]$aic,
      bic = results[[dist_name]]$bic,
      ks_statistic = ks_test$ks,
      ks_pvalue = ks_test$kstest
    )
  }
  
  return(gof_results)
}

#' Select best distribution based on hierarchical criteria
#' @keywords internal
select_best_distribution <- function(gof_results) {
  # Priority 1: Theoretical plausibility (handled by user)
  # Priority 2: Visual diagnostics (Q-Q plot, P-P plot)
  # Priority 3: Information criteria (AIC/BIC)
  # Priority 4: Statistical tests (K-S test)
  
  # Select distribution with lowest BIC
  bic_values <- sapply(gof_results, function(x) x$bic)
  best_dist_name <- names(which.min(bic_values))
  
  return(list(
    distribution = best_dist_name,
    parameters = gof_results[[best_dist_name]],
    aic = gof_results[[best_dist_name]]$aic,
    bic = gof_results[[best_dist_name]]$bic,
    gof_metrics = gof_results[[best_dist_name]]
  ))
}

#' Assign quality control flags (Red/Yellow/Green)
#' @keywords internal
assign_quality_flags <- function(best_dist, gof_results) {
  ks_pvalue <- best_dist$gof_metrics$ks_pvalue
  
  if (ks_pvalue < 0.01) {
    flag <- "RED"
    message <- "Poor fit: Consider alternative distributions or transformations"
  } else if (ks_pvalue < 0.05) {
    flag <- "YELLOW"
    message <- "Marginal fit: Proceed with caution and sensitivity analysis"
  } else {
    flag <- "GREEN"
    message <- "Good fit: Distribution is appropriate"
  }
  
  return(list(flag = flag, message = message, ks_pvalue = ks_pvalue))
}

#' Generate Q-Q plot for visual diagnostics
#' @export
plot_qq <- function(data, distribution, parameters) {
  qqplot_data <- qqplot(data, distribution, parameters)
  print(qqplot_data)
}

#' Generate audit report
#' @export
generate_audit_report <- function(audit_log, output_file = "dsa_audit_report.txt") {
  sink(output_file)
  
  cat("=== DSA Algorithm Audit Report ===\n\n")
  cat("Timestamp:", as.character(audit_log$timestamp), "\n")
  cat("Estimand:", audit_log$estimand, "\n")
  cat("Sample size:", audit_log$sample_size, "\n\n")
  
  cat("Data Summary:\n")
  print(audit_log$data_summary)
  cat("\n")
  
  cat("Distributions Tested:", paste(audit_log$distributions_tested, collapse = ", "), "\n\n")
  
  cat("Selected Distribution:", audit_log$selected_distribution$distribution, "\n")
  cat("Quality Flag:", audit_log$qc_flags$flag, "\n")
  cat("Quality Message:", audit_log$qc_flags$message, "\n\n")
  
  cat("Completion Time:", as.character(audit_log$completion_time), "\n")
  
  sink()
  
  cat("Audit report saved to:", output_file, "\n")
}

# Example usage (commented out)
# set.seed(123)
# data <- rnorm(1000, mean = 50, sd = 10)
# result <- dsa_algorithm(data, estimand = "Mean treatment effect")
# generate_audit_report(result$audit_log)
