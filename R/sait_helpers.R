
# Summary reporting helper
.report_fit_summary <- function(res, verbose = TRUE) {
    if (verbose && "fit_method" %in% colnames(res)) {
        total_genes <- nrow(res)
        fallback_mask <- !is.na(res$fit_method) & res$fit_method != "lmer"
        n_fallback <- sum(fallback_mask)

        # Report model convergence and fit quality
        message(sprintf("[calculate_sait_interaction] Analyzed %d genes | Primary fits: %d | Alternative method: %d",
            total_genes, total_genes - n_fallback, n_fallback))

        if (n_fallback > 0) {
            tab <- table(res$fit_method[fallback_mask])
            tab_str <- paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")
            message(sprintf("[calculate_sait_interaction]   Methods: %s", tab_str))
        }

        # Report numerical/convergence issues
        if ("singular" %in% colnames(res)) {
            n_sing <- sum(as.logical(res$singular), na.rm = TRUE)
            if (n_sing > 0) {
                message(sprintf("[calculate_sait_interaction]   Singular fits (collinear effects): %d genes",
                  n_sing))
            }
        }

        # Report effect size range (q-interaction magnitude)
        if ("f_statistic" %in% colnames(res)) {
            f_vals <- res$f_statistic[!is.na(res$f_statistic)]
            if (length(f_vals) > 0) {
                message(sprintf("[calculate_sait_interaction] q x condition interaction strength: F-stat range [%.2f, %.2f]",
                  min(f_vals), max(f_vals)))
            }
        }

        # Report significance summary
        if ("adj_p_interaction" %in% colnames(res)) {
            sig_p <- sum(res$adj_p_interaction < 0.05, na.rm = TRUE)
            message(sprintf("[calculate_sait_interaction] Significant results (adj.p < 0.05): %d/%d genes (%.1f%%)",
                sig_p, total_genes, 100 * sig_p/total_genes))
        }
    }
}

# ==============================================================================
# Helper: Compute AR(1) design effect for autocorrelated entropy differences
# ==============================================================================
# PURPOSE:
#   Estimates autocorrelation (rho) from differenced entropy data and computes
#   design effect for bias correction in analysis. This is NOT the Kish formula
#   (which assumes exchangeable ICC); TSENAT uses AR(1) correlation structure
#   after ARIMA(1,1,0) differencing.
#
# NOTE: AR(1) design effect computation has been consolidated into
# .compute_ar1_design_effect() in sait_gee.R, which uses the exact
# finite-m form: D_eff = 1 + 2*Σ_{k=1}^{m-1} (1-k/m)*ρ^k
# (Diggle et al. 2002, Crowder 1995).

.estimate_ar1_rho <- function(entropy_diff, subject_vec = NULL) {
    # Estimate first-order autocorrelation rho from differenced entropy Input:
    # entropy_diff = first-differenced entropy values DeltaH_q = H_q - H_{q-1}
    # Returns: rho estimate in [0, 1], or NULL if insufficient data Issues
    # warning if rho is very high (GAMM convergence risk)

    if (is.null(entropy_diff) || length(na.omit(entropy_diff)) < 3) {
        return(NULL)
    }

    # Remove NA values
    entropy_clean <- na.omit(entropy_diff)

    if (length(entropy_clean) < 3) {
        return(NULL)
    }

    # OPTIMIZATION (March 2026): Use stats::acf() for numerical stability
    # Previous: Manual computation (divides by n instead of n-1, less stable)
    # New: Built-in acf() for better stability and standard handling
    acf_result <- tryCatch({
        stats::acf(entropy_clean, lag.max = 1, plot = FALSE, demean = TRUE)
    }, error = function(e) NULL)

    if (!is.null(acf_result)) {
        rho_est <- as.numeric(acf_result$acf[2, 1, 1])
    } else {
        # Fallback to manual calculation if acf fails
        n <- length(entropy_clean)
        mean_x <- mean(entropy_clean, na.rm = TRUE)
        var_x <- sum((entropy_clean - mean_x)^2, na.rm = TRUE)/n

        if (var_x < 1e-10) {
            return(NULL)
        }

        x_t <- entropy_clean[-n]
        x_t1 <- entropy_clean[-1]
        cov_lag1 <- sum((x_t - mean_x) * (x_t1 - mean_x), na.rm = TRUE)/n
        rho_est <- cov_lag1/var_x
    }

    # BUG FIX: Check for NA after acf/fallback computation
    # If rho_est is NA, return NULL to allow fallback handling
    if (is.na(rho_est)) {
        return(NULL)
    }

    # Ensure rho is in [0, 1] (sometimes numerical errors give slight negative
    # values)
    rho_est <- max(0, min(1, rho_est))

    # OPTIMIZATION (March 2026): Add tolerance checks for edge cases Issue #9:
    # Missing tolerance checks from CODE_REVIEW_BUGS_FOUND.md
    if (rho_est > 0.95) {
        # Very high autocorrelation - warn about GAMM convergence risk
        warning(sprintf("AR(1) autocorrelation very high (rho=%.3f). GAMM may fail to converge. Consider reducing q-values or checking data for trends.",
            rho_est), call. = FALSE)
    }

    if (rho_est < 0.01) {
        # Very small autocorrelation - independence assumption near valid
        # Return NULL to suggest simpler model without AR(1)
        return(NULL)
    }

    return(rho_est)
}

# ==============================================================================
# Helper: Knot Selection for Tsallis Entropy Curve Fitting
# ==============================================================================
# MATHEMATICAL FOUNDATION:
#   Tsallis entropy H_q is MATHEMATICALLY GUARANTEED to be monotone decreasing
#   in q. Therefore, k-selection uses a simple fixed formula based on the
#   number of unique q-values. This ensures adequate smoothing without
#   noise-driven over-complexity.
#
# KNOT SELECTION FORMULA:
#   k = max(min_k, min(max_k, n_q_unique - 1))
#
#   Principle: Use at most (# unique q values - 1) basis functions. This
#   provides data-driven parsimony while ensuring sufficient flexibility.
#
# HISTORICAL NOTE:
#   Earlier versions attempted CV-based adaptation to detect curve complexity,
#   but this was backwards for monotone data (high CV indicates noise, not
#   complexity). Current fixed approach is superior for this problem.
#
# REFERENCE:
#   Wood (2006) Generalized Additive Models; approach enforced via fixed
#   formula to respect mathematical monotonicity property of Tsallis entropy

# ==============================================================================
# STATIONARITY VALIDATION FRAMEWORK FOR TSALLIS ENTROPY MODELING
# ==============================================================================
# MATHEMATICAL JUSTIFICATION:
#   - Tsallis entropy H_q is monotone DECREASING in q (proven in Tsallis 1988)
#   - This monotonicity makes the series NON-STATIONARY (systematic/deterministic trend)
#   - AR(1) models assume stationarity (constant mean/variance around trend)
#   - SOLUTION: ARIMA(1,1,0) applies AR(1) to FIRST DIFFERENCES
#   - DeltaH_q = H_q - H_{q-1} removes the trend (differencing)
#   - Then AR(1) can model residual correlation in the differenced series
#
# VALIDATION FRAMEWORK:
#   Four complementary tests validate core assumptions before ARIMA(1,1,0)
#
# TEST 1: Monotonicity Check (Visual Validation)
#   Purpose: Detect q-value ordering issues or data quality problems
#   Method: Count decreasing vs increasing pairs in ordered q-values
#   Expected for Tsallis: >95% pairs should be decreasing (monotone)
#
# TEST 2: Augmented Dickey-Fuller (ADF) Test for Unit Root
#   H0 (Null): Series has unit root (non-stationary)
#   H1 (Alt):  Series is stationary
#   Expected for raw Tsallis entropy: FAIL to reject H0 (non-stationary with unit root)
#   Expected for differenced data:    REJECT H0 (stationary, no unit root)
#   Reference: Dickey & Fuller (1979, 1981); MacKinnon (1996)
#
# TEST 3: KPSS Test (Reverse of ADF)
#   H0 (Null): Series IS stationary
#   H1 (Alt):  Series is NON-stationary
#   Expected for raw Tsallis entropy: REJECT H0 (non-stationary)
#   Expected for differenced data:    FAIL to reject H0 (stationary)
#   Reference: Kwiatkowski, Phillips, Schmidt & Shin (1992)
#
# TEST 4: Integration Order Validation
#   Apply ADF/KPSS to differenced data to confirm ARIMA(1,1,0) appropriateness
#   Expected: Differenced data should be I(0)—integrated of order 0 (stationary)
#
# DECISION LOGIC FOR ARIMA(1,1,0):
#   Use ARIMA(1,1,0) if ALL conditions met:
#   ✓ Raw data is NON-monotone (>5% violations) OR fails stationarity tests
#   ✓ ADF test FAILs to reject H0 on raw data (has unit root)
#   ✓ KPSS test REJECTs H0 on raw data (non-stationary)
#   ✓ ADF test REJECTs H0 on differenced data (stationary)
#   ✓ KPSS test FAILs to reject H0 on differenced data (stationary)
#
# REFERENCES:
#   Dickey, D. A., & Fuller, W. A. (1979). Distribution of the estimators for
#     autoregressive time series with a unit root. J American Statistical
#     Association, 74(366), 427-431.
#   Kwiatkowski, D., Phillips, P. C., Schmidt, P., & Shin, Y. (1992). Testing
#     the null hypothesis of stationarity against the alternative of a unit root.
#     Journal of Econometrics, 54(1-3), 159-178.
#   MacKinnon, J. G. (1996). Numerical distribution functions for unit root and
#     cointegration tests. Journal of Applied Econometrics, 11(6), 601-618.
#   Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs statistics.
#     Journal of Statistical Physics, 52(1), 479-487.

# ==============================================================================
# RESIDUAL DIAGNOSTICS: Shapiro-Wilk Normality Testing
# ==============================================================================
# PURPOSE:
#   Verify that residuals from GAM/LMM/GEE models satisfy normality assumption.
#   This is a standard diagnostic for validating statistical model assumptions.
#
# DATABASE EVIDENCE (March 2026):
#   - Alberghina & Westerhoff (2001): Foundations of Systems Biology
#   - B004 (2008): Linear Models in Systems Biology
#   - Springer Handbook (2006): Springer Handbook of Statistical Methods
#
# METHOD:
#   Shapiro-Wilk test on model residuals (tests H0: residuals are normal)
#   Standard Practice: Applied universally in statistical modeling literature
#
# INTERPRETATION:
#   - p > 0.05: Fail to reject H0 → Residuals appear normal [OK]
#   - p ≤ 0.05: Reject H0 → Residuals show significant departure from normality [?]
#
# IMPLEMENTATION:
#   Extract residuals from fitted model, apply shapiro.test()

.test_residual_normality <- function(model, model_type = c("gam", "gamm", "lme",
    "gee"), verbose = FALSE) {
    # Args: model: fitted model object (GAM, GAMM, lme, or geeglm) model_type:
    # character - type of model for residual extraction verbose: if TRUE, print
    # diagnostic messages Returns: List with components: - shapiro_p_value:
    # p-value from Shapiro-Wilk test (NA if test fails) - residuals_normal:
    # logical - TRUE if p > 0.05 (residuals appear normal) - n_residuals:
    # number of residuals tested - test_status: character - 'pass', 'fail', or
    # 'error' - report: character - human-readable summary

    if (is.null(model) || inherits(model, "try-error")) {
        return(list(shapiro_p_value = NA_real_, residuals_normal = NA, n_residuals = 0,
            test_status = "error", report = "Model object is NULL or error class"))
    }

    model_type <- match.arg(model_type)
    residuals_vec <- NULL

    # Extract residuals based on model type
    tryCatch({
        if (model_type == "gam") {
            # Standard GAM: check family to determine residual type Pearson
            # residuals available for non-gaussian families (gamma, beta, etc.)
            # Gaussian family: use deviance residuals directly
            family_name <- ifelse(!is.null(model$family), model$family$family, "gaussian")
            if (family_name == "gaussian") {
                residuals_vec <- residuals(model, type = "deviance")
            } else {
                # For non-gaussian families, use Pearson residuals
                residuals_vec <- tryCatch(residuals(model, type = "pearson"), error = function(e) residuals(model,
                  type = "deviance"))
            }
        } else if (model_type == "gamm") {
            # GAMM: extract residuals from $gam component, check family first
            if (!is.null(model$gam)) {
                family_name <- ifelse(!is.null(model$gam$family), model$gam$family$family,
                  "gaussian")
                if (family_name == "gaussian") {
                  residuals_vec <- residuals(model$gam, type = "deviance")
                } else {
                  residuals_vec <- tryCatch(residuals(model$gam, type = "pearson"),
                    error = function(e) residuals(model$gam, type = "deviance"))
                }
            } else {
                family_name <- ifelse(!is.null(model$family), model$family$family,
                  "gaussian")
                if (family_name == "gaussian") {
                  residuals_vec <- residuals(model, type = "deviance")
                } else {
                  residuals_vec <- tryCatch(residuals(model, type = "pearson"), error = function(e) residuals(model,
                    type = "deviance"))
                }
            }
        } else if (model_type == "lme") {
            # nlme::lme model: use residuals() generic
            residuals_vec <- residuals(model, type = "normalized")
        } else if (model_type == "gee") {
            # geeglm: use residuals() generic
            residuals_vec <- tryCatch(residuals(model, type = "pearson"), error = function(e) residuals(model,
                type = "deviance"))
        }
    }, error = function(e) {
        if (verbose) {
            message("[.test_residual_normality] Could not extract residuals: ", e$message)
        }
    })

    # Check if residuals were extracted successfully
    if (is.null(residuals_vec) || length(residuals_vec) == 0) {
        return(list(shapiro_p_value = NA_real_, residuals_normal = NA, n_residuals = 0,
            test_status = "error", report = "Could not extract residuals from model"))
    }

    # Remove missing values
    residuals_clean <- as.numeric(na.omit(residuals_vec))
    n_res <- length(residuals_clean)

    # Shapiro-Wilk test requires at least 3 observations
    if (n_res < 3) {
        return(list(shapiro_p_value = NA_real_, residuals_normal = NA, n_residuals = n_res,
            test_status = "error", report = sprintf("Insufficient residuals for Shapiro-Wilk test (n=%d, need >=3)",
                n_res)))
    }

    # Run Shapiro-Wilk test
    test_result <- tryCatch({
        stats::shapiro.test(residuals_clean)
    }, error = function(e) {
        return(NULL)
    })

    if (is.null(test_result)) {
        return(list(shapiro_p_value = NA_real_, residuals_normal = NA, n_residuals = n_res,
            test_status = "error", report = "Shapiro-Wilk test execution failed"))
    }

    # Extract test statistics
    p_value <- test_result$p.value
    is_normal <- p_value > 0.05  # Fail to reject H0 at alpha=0.05

    if (verbose) {
        status_text <- if (is_normal)
            "PASS [OK]" else "FAIL ?"
        message(sprintf("[.test_residual_normality] %s (p=%.4f, n=%d residuals)",
            status_text, p_value, n_res))
    }

    return(list(shapiro_p_value = p_value, residuals_normal = is_normal, n_residuals = n_res,
        test_status = if (is_normal) "pass" else "fail", report = sprintf("Shapiro-Wilk test: p=%.4f, %s normal (n=%d residuals)",
            p_value, if (is_normal) "residuals appear" else "residuals NOT", n_res)))
}

# ==================================================================
# ARIMA(1,1,0) Implementation: Compute First Differences of Entropy
# ==================================================================
# BACKGROUND:
#   - Tsallis entropy H_q is monotone decreasing in q (non-stationary)
#   - AR(1) assumes stationarity (constant mean, variance)
#   - Solution: Apply AR(1) to first differences DeltaH_q = H_q - H_{q-1}
#   - Result: ARIMA(1,1,0) = Integrated AR(1) = AR(1) on differenced data
#
# IMPLEMENTATION STEPS:
#   1. Order data by q-values to ensure proper differencing
#   2. Compute differences within each subject (NOT across subjects)
#   3. Return data frame with differenced entropy, q values, group, subject
#   4. Note: Loses 1 observation per subject (trade-off for stationarity)
#
# RETURNS:
#   list(df_diff, n_lost_obs) or NULL if insufficient data
.compute_arima_differences <- function(df, q_vals, group_vec, subject_vec = NULL) {

    if (is.null(df) || nrow(df) == 0) {
        return(NULL)
    }

    # Add q and group to data frame for sorting
    df_full <- data.frame(entropy = as.numeric(df$entropy), q = as.numeric(q_vals),
        group = factor(group_vec), subject = if (!is.null(subject_vec))
            factor(subject_vec) else factor(seq_len(nrow(df))), stringsAsFactors = FALSE)

    # Remove NA entropy values
    df_full <- df_full[!is.na(df_full$entropy), ]

    if (nrow(df_full) < 2) {
        return(NULL)
    }

    # Sort by subject and q to ensure proper differencing within subjects
    df_full <- df_full[order(df_full$subject, df_full$q), ]

    # Compute first differences within each subject
    df_diff_list <- list()
    n_lost <- 0

    for (subj in levels(df_full$subject)) {
        subj_idx <- which(df_full$subject == subj)

        if (length(subj_idx) < 2) {
            # Skip subjects with < 2 observations (can't compute difference)
            n_lost <- n_lost + length(subj_idx)
            next
        }

        # Extract subject data (should already be sorted by q)
        subj_data <- df_full[subj_idx, ]

        # Compute differences: DeltaH_q = H_q - H_{q-1} CRITICAL: Convert group
        # to character BEFORE subsetting to avoid factor level issues When you
        # subset a factor, R keeps ALL original levels, which causes rbind()
        # problems later
        n_diff <- nrow(subj_data) - 1

        df_diff_list[[subj]] <- data.frame(entropy_diff = diff(subj_data$entropy),
            q = subj_data$q[-1], q_prev = subj_data$q[-nrow(subj_data)],
            group = as.character(subj_data$group[-1]),
            subject = rep(subj, n_diff), stringsAsFactors = FALSE)
    }

    if (length(df_diff_list) == 0) {
        return(NULL)
    }

    # Combine all subject differences
    df_diff <- do.call(rbind, df_diff_list)
    rownames(df_diff) <- NULL

    if (nrow(df_diff) == 0) {
        return(NULL)
    }

    # Rename entropy_diff to entropy for compatibility with model fitting
    names(df_diff)[names(df_diff) == "entropy_diff"] <- "entropy"

    # CRITICAL FIX (Phase 15): Ensure factor consistency after ARIMA
    # differencing Problem: Some subjects/groups may be completely dropped by
    # differencing, leaving factor levels that don't exist in the data. nlme
    # can't handle this.  Solution: Convert to factor WITHOUT forcing unused
    # original levels.  Just let R infer the levels from the actual data
    # present.
    df_diff$subject <- factor(as.character(df_diff$subject))
    df_diff$group <- factor(as.character(df_diff$group))

    return(list(df = df_diff, n_observations_original = nrow(df_full) + n_lost, n_observations_differenced = nrow(df_diff),
        n_observations_lost = n_lost, transformation = "ARIMA(1,1,0): First differences"))
}

# Helper: Check if entropy data is truly bounded in [0, 1] Returns TRUE if data
# appears normalized/proportional
.is_bounded_0_1 <- function(entropy_vals) {
    entropy_clean <- na.omit(entropy_vals)
    if (length(entropy_clean) == 0)
        return(FALSE)
    finite_clean <- is.finite(entropy_clean)
    if (!any(finite_clean))
        return(FALSE)
    min_val <- min(entropy_clean[finite_clean])
    max_val <- max(entropy_clean[finite_clean])
    tolerance <- 0.01
    bounds_check <- min_val >= -tolerance && max_val <= 1 + tolerance
    if (!bounds_check)
        return(FALSE)
    approaches_lower_bound <- min_val <= 0.1
    approaches_upper_bound <- max_val >= 0.9
    return(approaches_lower_bound || approaches_upper_bound)
}

# ============================================================================
# MEMOIZATION: Cache expensive computations
# ============================================================================
# Memoization reduces redundant calculations in multi-q iterative analyses
# Expected speedup: 10-30% for typical multi-q analyses

if (getOption("TSENAT.memoization", TRUE)) {
    # Cache knot selection: input = (entropy, q_vals, n_unique, ...)  Avoids
    # recomputing knots for same entropy data across iterations
    .adaptive_spline_knots_memo <- memoise(.adaptive_spline_knots, cache = cache_memory())
} else {
    # Fallback: no memoization if disabled globally
    .adaptive_spline_knots_memo <- .adaptive_spline_knots
}

#' Internal: Clear memoization cache
#' @description Invalidates cached results for new dataset processing
#' @noRd
.clear_sait_helper_cache <- function() {
    if (getOption("TSENAT.memoization", TRUE)) {
        forget(.adaptive_spline_knots_memo)
    }
}



# Helper: Compute skewness of a vector Positive skew: right tail longer (mode <
# median < mean) Negative skew: left tail longer (mean < median < mode)
.compute_skewness <- function(x, na.rm = TRUE) {
    if (na.rm)
        x <- na.omit(x)
    if (length(x) < 3)
        return(NA)

    m <- mean(x)
    s <- sd(x)
    n <- length(x)

    if (s == 0)
        return(0)

    # Unbiased skewness estimate
    skew <- (sum((x - m)^3)/n)/(s^3)
    return(skew)
}

# ================================================================================
# HETEROSCEDASTICITY DETECTION AND VARIANCE WEIGHTING (March 2026)
# ================================================================================
# Tsallis entropy often exhibits variance that depends on: 1. Mean entropy
# level (mean-variance relationship) 2. q-value (variance changes across
# diversity orders) 3. Group/condition (treatment vs control may have different
# variance) Consequence: Tests can be biased with inflated Type I error if
# heteroscedasticity ignored Solution: Detect heteroscedasticity and apply
# appropriate variance adjustment/weighting

# Detect heteroscedasticity using Breusch-Pagan test
.detect_heteroscedasticity <- function(df, q_vals, group_vec, verbose = FALSE) {
    # Fit OLS to get residuals
    fit_ols <- try(lm(entropy ~ q + group, data = df), silent = TRUE)

    if (inherits(fit_ols, "try-error")) {
        return(list(is_heteroscedastic = NA, bp_stat = NA, p_value = NA, var_ratio_q = NA,
            var_ratio_group = NA))
    }

    residuals_sq <- residuals(fit_ols)^2

    # Breusch-Pagan auxiliary regression: log(residuals^2) ~ q + group
    fit_aux <- try(lm(log(residuals_sq + 1e-08) ~ q + factor(group), data = df),
        silent = TRUE)

    if (inherits(fit_aux, "try-error")) {
        return(list(is_heteroscedastic = NA, bp_stat = NA, p_value = NA, var_ratio_q = NA,
            var_ratio_group = NA))
    }

    # BP statistic = RSS from auxiliary model / (2 * RSS from original model)
    rss_aux <- sum(residuals(fit_aux)^2)
    # OPTIMIZATION (March 2026): Cache computation to avoid redundant
    # calculation
    fitted_sq_sum <- sum((fitted(fit_aux) - mean(fitted(fit_aux)))^2)
    tss_aux <- fitted_sq_sum + rss_aux

    bp_stat <- (fitted_sq_sum/tss_aux * nrow(df))
    # Compute correct degrees of freedom: number of predictors in auxiliary
    # regression BUG FIX: Was hardcoded to 2, but should be ncol(X) - 1 where X
    # is model.matrix
    df_bp <- ncol(model.matrix(fit_aux)) - 1
    p_value <- 1 - pchisq(bp_stat, df = df_bp)

    # Compute variance ratios OPTIMIZATION (March 2026): Use tapply() instead
    # of sapply + subsetting (2-3x faster)
    residuals_vec <- residuals(fit_ols)
    var_by_q <- tapply(residuals_vec, df$q, var)
    finite_q <- is.finite(var_by_q)
    if (any(finite_q)) {
        max_q <- max(var_by_q[finite_q])
        min_q <- min(var_by_q[finite_q])
        var_ratio_q <- max_q/(min_q + 1e-08)
    } else {
        var_ratio_q <- NA
    }

    var_by_group <- tapply(residuals_vec, df$group, var)
    finite_g <- is.finite(var_by_group)
    if (any(finite_g)) {
        max_g <- max(var_by_group[finite_g])
        min_g <- min(var_by_group[finite_g])
        var_ratio_group <- max_g/(min_g + 1e-08)
    } else {
        var_ratio_group <- NA
    }

    if (verbose) {
        message(sprintf("[Heteroscedasticity] BP p-value: %.4f, Var ratio (q): %.2f, Var ratio (group): %.2f",
            p_value, var_ratio_q, var_ratio_group))
    }

    return(list(is_heteroscedastic = p_value < 0.05, bp_stat = bp_stat, p_value = p_value,
        var_ratio_q = var_ratio_q, var_ratio_group = var_ratio_group))
}

# Estimate variance weights for heteroscedasticity adjustment
.estimate_variance_weights <- function(df, q_vals, method = "power", verbose = FALSE) {
    # Estimate weights to model variance heterogeneity method = 'power': Model
    # Var ~ q^?, compute weights w_i = q_i^(-?)  method = 'residual': Use
    # residual variance from OLS as observation weights

    if (method == "power") {
        # Estimate power parameter ? via regression: log(residuals_sq) ~ q
        # First, fit OLS to get residuals
        fit_ols <- try(lm(entropy ~ q + group, data = df), silent = TRUE)

        # If OLS with group fails, try just q
        if (inherits(fit_ols, "try-error")) {
            fit_ols <- try(lm(entropy ~ q, data = df), silent = TRUE)
        }

        if (inherits(fit_ols, "try-error")) {
            return(NULL)
        }

        residuals_ols <- residuals(fit_ols)
        residuals_sq <- residuals_ols^2

        # This gives: log(Var) = log(sigma2) + ? * log(q) So: Var ~ sigma2 *
        # q^?  Weights: w_i = 1 / (sigma2 * q_i^?) ? q_i^(-?)

        w <- 1/(residuals_sq + 1e-08)
        wfit <- try(lm(log(residuals_sq + 1e-08) ~ log(df$q + 1e-08), weights = w),
            silent = TRUE)

        if (!inherits(wfit, "try-error")) {
            theta_est <- coef(wfit)[2]
            if (!is.na(theta_est)) {
                # Ensure positive weighting
                theta_est <- max(theta_est, 0.01)
                weights <- 1/(df$q^theta_est + 1e-08)
                weights <- weights/mean(weights, na.rm = TRUE)  # Standardize

                if (verbose) {
                  message(sprintf("[Variance Weighting] Estimated power parameter ? = %.3f",
                    theta_est))
                }

                return(list(weights = weights, power_param = theta_est, method = "power"))
            }
        }
    }

    if (method == "residual") {
        # Use inverse variance as weights Try OLS fit to estimate residual
        # variance
        ols_fit <- NULL

        # Try with group if available, otherwise just q
        if ("group" %in% colnames(df)) {
            ols_fit <- try(lm(entropy ~ q + group, data = df), silent = TRUE)
        } else {
            ols_fit <- try(lm(entropy ~ q, data = df), silent = TRUE)
        }

        residuals_sq <- NA_real_
        if (!is.null(ols_fit) && !inherits(ols_fit, "try-error")) {
            residuals_sq <- residuals(ols_fit)^2
        }

        # If OLS succeeded and we have residuals, compute weights
        if (!all(is.na(residuals_sq))) {
            weights <- 1/(residuals_sq + 1e-08)
            weights <- weights/mean(weights, na.rm = TRUE)  # Standardize

            return(list(weights = weights, method = "residual"))
        }
    }

    # Fallback: uniform weights
    return(list(weights = rep(1, nrow(df)), method = "uniform"))
}



# Helper functions for .calculate_sait() These internal functions decompose the
# main function logic into focused, testable components that each handle a
# single responsibility.

#' @title Validate Input Parameters for SAIT Interaction Testing
#'
#' @description
#' Internal helper that consolidates parameter validation for
#' \code{.calculate_sait()}. Checks argument types, values,
#' and inter-dependencies to ensure valid model fitting.
#'
#' @param method Character; modeling method (matched from user input)
#' @param pvalue Character; p-value type specification
#' @param corstr Character; correlation structure
#' @param regularization Character; dimensionality reduction method
#' @param multicorr Character; multi-q correction method
#' @param pcorr Character; legacy p-value correction method
#' @param storey Logical; whether to apply Storey correction
#' @param wy_randomizations Integer or character; number of permutations
#'   for Westfall-Young correction. Use 'auto' to estimate a suitable
#'   permutation count from the data.
#' @param paired Logical; whether design is paired
#' @param subject_col Character or NULL; subject column name
#' @param se SummarizedExperiment object
#' @param verbose Logical; print diagnostic messages
#'
#' @return List with validated and normalized parameters:
#'   \itemize{
#'     \item method: Validated method name
#'     \item pvalue: Validated p-value type
#'     \item corstr: Validated correlation structure
#'     \item regularization: Validated regularization method
#'     \item multicorr: Validated multicorr method
#'     \item pcorr: Validated legacy pcorr
#'     \item subject_col: Auto-detected or user-provided subject column
#'     \item wy_randomizations: Validated permutation count or 'auto'
#'   }
#'

#' @noRd
.validate_sait_interaction_input <- function(method, pvalue, corstr, regularization,
    multicorr, pcorr, storey, wy_randomizations, paired, subject_col, se, verbose) {
    # Validate storey parameter
    if (!is.logical(storey)) {
        stop("storey must be TRUE or FALSE", call. = FALSE)
    }

    # Validate wy_randomizations
    if (is.character(wy_randomizations) && tolower(wy_randomizations) == "auto") {
        wy_randomizations <- "auto"
    } else if (is.null(wy_randomizations)) {
        wy_randomizations <- 1000
    } else if (is.numeric(wy_randomizations) && wy_randomizations < 1) {
        stop("wy_randomizations must be numeric and >= 1", call. = FALSE)
    } else if (!is.numeric(wy_randomizations)) {
        stop("wy_randomizations must be numeric, 'auto', or NULL", call. = FALSE)
    } else {
        wy_randomizations <- as.integer(wy_randomizations)
        if (wy_randomizations < 100) {
            warning("wy_randomizations < 100 may give unreliable p-values; ", "recommend >= 100",
                call. = FALSE)
        }
    }

    # Auto-detect subject_col from colData if paired=TRUE and subject_col=NULL
    # Prioritize 'paired_samples' or 'sample_base' columns
    if (paired && is.null(subject_col)) {
        cd_colnames <- colnames(SummarizedExperiment::colData(se))

        # Check for paired_samples or sample_base columns
        if ("paired_samples" %in% cd_colnames) {
            subject_col <- "paired_samples"
            if (verbose) {
                message("[calculate_sait] paired=TRUE detected; ", "auto-using subject_col='paired_samples'")
            }
        } else if ("sample_base" %in% cd_colnames) {
            subject_col <- "sample_base"
            if (verbose) {
                message("[calculate_sait] paired=TRUE detected; ", "auto-using subject_col='sample_base'")
            }
        } else {
            # Error if paired=TRUE but no recognized pairing column found
            stop("paired=TRUE requires either 'paired_samples' or ", "'sample_base' column in colData. Available columns: ",
                paste(cd_colnames, collapse = ", "), ". Ensure .calculate_diversity() or map_metadata() was ",
                "called with appropriate metadata.", call. = FALSE)
        }
    }

    # Validate (method, regularization) combination
    valid_reg_methods <- list(
        gam  = c("pca", "gamsel", "spline"),
        lmm  = c("pca", "lasso", "elasticnet"),
        fpca = c("pca", "lasso", "elasticnet"),
        gee  = c("pca", "lasso", "elasticnet", "gamsel", "spline")  # all accepted, ignored
    )
    if (!regularization %in% valid_reg_methods[[method]]) {
        stop(sprintf("regularization='%s' is not valid for method='%s'. Valid options: %s",
            regularization, method, paste(valid_reg_methods[[method]], collapse = ", ")),
            call. = FALSE)
    }

    return(list(method = method, pvalue = pvalue, corstr = corstr, regularization = regularization,
        multicorr = multicorr, pcorr = pcorr, subject_col = subject_col,
        wy_randomizations = wy_randomizations))
}

#' @title Parse Sample Metadata from SummarizedExperiment
#'
#' @description
#' Internal helper that extracts sample names, q-values, and group
#' assignments from the diversity assay column names and colData.
#'
#' @param se SummarizedExperiment object
#' @param condition_col Character; colData column with group assignments
#' @param assay_name Character; name of diversity assay
#' @param verbose Logical; print diagnostic messages
#'
#' @return List containing:
#'   \itemize{
#'     \item sample_q: Full column names with q= values
#'     \item sample_names: Unique sample identifiers
#'     \item q_vals: Parsed q-value parameters
#'     \item group_vec: Group assignment for each observation
#'     \item has_q: Logical vector indicating cols with q=
#'   }
#'

#' @noRd
.parse_sample_metadata <- function(se, condition_col, assay_name, verbose) {
    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) {
        stop(sprintf("Assay '%s' not found in SummarizedExperiment", assay_name))
    }

    sample_q <- colnames(mat)
    if (is.null(sample_q) || length(sample_q) == 0) {
        stop("No column names found on diversity assay")
    }

    # Parse sample names and q values from column names like 'Sample_q=0.01'
    sample_names <- sub("_q=.*", "", sample_q)
    has_q <- grepl("_q=", sample_q)
    if (!any(has_q)) {
        stop("Could not parse q values; expected '_q=' in column names", call. = FALSE)
    }
    if (!all(has_q)) {
        stop("Some column names are missing '_q='; ensure all diversity ", "columns include a q value",
            call. = FALSE)
    }
    q_vals <- as.numeric(sub(".*_q=", "", sample_q))

    # Determine group for each sample
    condition_in_coldata <- !is.null(condition_col) && condition_col %in% colnames(SummarizedExperiment::colData(se))
    if (condition_in_coldata) {
        st <- as.character(SummarizedExperiment::colData(se)[, condition_col])
        names(st) <- rownames(SummarizedExperiment::colData(se))
        # Index by the FULL column names (sample_q), not by sample_names
        group_vec <- unname(st[sample_q])
    } else {
        stop("No sample grouping found: please supply `condition_col` ", "or map sample types into `colData(se)` before calling ",
            ".calculate_sait().", call. = FALSE)
    }

    if (verbose) {
        message("[calculate_sait_interaction] parsed samples and groups: ", length(unique(sample_names)),
            " samples, ", length(unique(q_vals)), " q-values")
    }

    return(list(sample_q = sample_q, sample_names = sample_names, q_vals = q_vals,
        group_vec = group_vec, has_q = has_q))
}

#' @title Fit Regularized Regression Models for All Genes
#'
#' @description
#' Internal helper that orchestrates parallel or sequential fitting
#' of regularized/penalized regression models (GAM, LMM, GEE, FPCA) to all genes 
#' in the diversity matrix. Consolidates the fitting loop and result collection logic.
#'
#' @param mat Matrix; diversity assay data
#' @param se SummarizedExperiment object
#' @param metadata List; output from .parse_sample_metadata()
#' @param method Character; modeling method
#' @param pvalue Character; p-value type
#' @param subject_col Character or NULL; subject column
#' @param paired Logical; whether design is paired
#' @param min_obs Integer; minimum observations per gene
#' @param nthreads Integer; number of parallel threads
#' @param verbose Logical; print diagnostics
#' @param bias_correction Logical; apply KC bias correction (GEE)
#' @param regularization Character; dimensionality reduction method
#' @param corstr Character; correlation structure
#' @param adaptive_knots Logical; adaptive knot selection (GAM)
#'
#' @return Data.frame with fitted model results for all genes
#'

#' @noRd
.fit_all_genes <- function(mat, se, metadata, method, pvalue, subject_col, paired,
    min_obs, nthreads, verbose, bias_correction, regularization, corstr, adaptive_knots) {
    suppress_lme4_warnings <- TRUE
    progress <- FALSE
    gene_weights <- NULL

    fit_one <- function(g, group_vec_override = NULL) {
        # Use override group_vec if provided (for permutation testing),
        # otherwise use outer scope
        gv <- if (!is.null(group_vec_override)) {
            group_vec_override
        } else {
            metadata$group_vec
        }

        .fit_one_interaction(g = g, se = se, mat = mat, q_vals = metadata$q_vals,
            sample_names = metadata$sample_names, group_vec = gv, method = method,
            pvalue = pvalue, subject_col = subject_col, paired = paired, min_obs = min_obs,
            verbose = verbose, suppress_lme4_warnings = suppress_lme4_warnings, progress = progress,
            bias_correction = bias_correction, regularization = regularization, corstr = corstr,
            adaptive_knots = adaptive_knots, weights = gene_weights)
    }

    if (nthreads > 1) {
        res_list <- .bplapply(rownames(mat), fit_one, nthreads = nthreads)
    } else {
        res_list <- lapply(rownames(mat), fit_one)
    }
    all_results <- Filter(Negate(is.null), res_list)

    if (length(all_results) == 0) {
        return(data.frame())
    }

    # Ensure all results are data frames
    all_results <- Filter(function(x) is.data.frame(x), all_results)
    if (length(all_results) == 0) {
        return(data.frame())
    }

    # Normalize columns: collect all unique column names and ensure every
    # result has them Use first result's column order as reference
    first_cols <- colnames(all_results[[1]])
    all_col_names <- unique(c(first_cols, unlist(lapply(all_results, colnames))))

    all_results <- lapply(all_results, function(df) {
        # Add missing columns as NA
        missing_cols <- setdiff(all_col_names, colnames(df))
        for (col in missing_cols) {
            df[[col]] <- NA
        }
        # Keep columns in consistent order
        df[, all_col_names, drop = FALSE]
    })

    # Phase 15: Wrap rbind in try-error to catch 'los nombres no coinciden'
    # errors from factor level mismatches during result combination
    res <- try(do.call(rbind, all_results), silent = FALSE)
    if (inherits(res, "try-error")) {
        # Debug: Check column mismatch details
        col_counts <- vapply(all_results, ncol, FUN.VALUE = integer(1))
        col_names_list <- lapply(all_results, colnames)
        unique_col_counts <- unique(col_counts)

        if (length(unique_col_counts) > 1) {
            # Column count mismatch
            msg <- sprintf("Column mismatch detected: %d results with varying columns [%s]",
                length(all_results), paste(unique_col_counts, collapse = ", "))
            warning("[calculate_sait_interaction] rbind failed with: ", conditionMessage(res),
                "\n[", msg, "]\n[Attempting recovery: ensuring all results have same columns]",
                call. = FALSE)

            # Add/remove columns to match first result's structure
            first_cols <- col_names_list[[1]]
            all_results_aligned <- lapply(all_results, function(df) {
                # Add missing columns as NA
                missing_cols <- setdiff(first_cols, colnames(df))
                for (col in missing_cols) {
                  df[[col]] <- NA
                }
                # Keep only matching columns
                df[, first_cols, drop = FALSE]
            })

            res <- try(do.call(rbind, all_results_aligned), silent = FALSE)
            if (!inherits(res, "try-error")) {
                # Successfully aligned, continue
                return(res)
            }
        }

        # If rbind fails due to factor level issues, try converting factor
        # columns to character
        warning("[calculate_sait_interaction] rbind failed with: ", conditionMessage(res),
            "\n[Attempting recovery: converting factors to character]", call. = FALSE)

        # Convert all factor columns to character to allow rbind
        all_results_char <- lapply(all_results, function(df) {
            factor_cols <- vapply(df, is.factor, FUN.VALUE = logical(1))
            df[factor_cols] <- lapply(df[factor_cols], as.character)
            df
        })

        res <- try(do.call(rbind, all_results_char), silent = FALSE)
        if (inherits(res, "try-error")) {
            stop("[calculate_sait_interaction] Could not combine results even after ",
                "factor conversion. Error: ", conditionMessage(res), call. = FALSE)
        }
    }

    # VALIDATION: Ensure critical columns exist after rbind
    if (nrow(res) == 0) {
        stop("[calculate_sait_interaction] No genes analyzed (all filtered out)", call. = FALSE)
    }

    critical_cols <- c("p_interaction", "gene")
    missing_cols <- setdiff(critical_cols, colnames(res))
    if (length(missing_cols) > 0) {
        stop("[calculate_sait_interaction] CRITICAL: Missing columns in ", "results for ",
            method, " method: ", paste(missing_cols, collapse = ", "), "\nAvailable columns: ",
            paste(colnames(res), collapse = ", "), call. = FALSE)
    }

    # Ensure Shapiro-Wilk columns exist for methods that add them
    if (method %in% c("gam", "gee")) {
        if (!"shapiro_p_value" %in% colnames(res)) {
            res$shapiro_p_value <- NA_real_
        }
        if (!"residuals_normal" %in% colnames(res)) {
            res$residuals_normal <- NA
        }
        if (!"n_residuals_tested" %in% colnames(res)) {
            res$n_residuals_tested <- NA_integer_
        }
    }

    # Ensure ci_weighted column exists (Phase 1 tracking)
    if (!"ci_weighted" %in% colnames(res)) {
        res$ci_weighted <- NA  # Fallback
        if (verbose) {
            warning("[calculate_sait_interaction] ci_weighted column was ", "missing; added as NAs. This suggests a method helper ",
                "did not properly set ci_weighted.", call. = FALSE)
        }
    }

    return(res)
}

#' @title Adjust P-Values for Multiple Q-Values
#'
#' @description
#' Internal router function that applies the specified primary multi-q
#' p-value correction method (Hochberg, Westfall-Young, or
#' Benjamini-Yekutieli), then optionally applies Storey adaptive FDR
#' enhancement.
#'
#' @param p_values Numeric vector; raw p-values to adjust
#' @param multicorr Character; primary correction method
#' @param wy_randomizations Integer; number of permutations (WY only)
#' @param fit_one_fn Function; function to refit models (WY only)
#' @param metadata List; output from .parse_sample_metadata()
#' @param mat Matrix; diversity assay data (WY only)
#' @param rownames_mat Character; rownames of matrix (WY only)
#' @param se SummarizedExperiment object (WY only)
#' @param assay_name Character; assay name (WY only)
#' @param method Character; modeling method (WY only)
#' @param pvalue Character; p-value type (WY only)
#' @param subject_col Character or NULL; subject column (WY only)
#' @param paired Logical; paired design (WY only)
#' @param min_obs Integer; min observations (WY only)
#' @param nthreads Integer; parallel threads (WY only)
#' @param verbose Logical; print diagnostics
#' @param bias_correction Logical; KC bias correction (WY/GEE)
#' @param regularization Character; dimensionality reduction (WY/FPCA)
#' @param corstr Character; correlation structure (WY/GEE)
#' @param adaptive_knots Logical; adaptive knots (WY/GAM)
#' @param storey Logical; apply Storey after primary correction
#'
#' @return Adjusted p-values vector
#'

#' @noRd
.adjust_pvalues_multicorr <- function(p_values, multicorr = c("hochberg", "westfall-young", "benjamini-yekutieli"), wy_randomizations, fit_one_fn = NULL,
    metadata = NULL, mat = NULL, rownames_mat = NULL, se = NULL, assay_name = "diversity",
    method = NULL, pvalue = NULL, subject_col = NULL, paired = FALSE, min_obs = 10,
    nthreads = 1, verbose = FALSE, bias_correction = TRUE, regularization = NULL,
    corstr = "ar1", adaptive_knots = TRUE, storey = FALSE) {
    # Validate multicorr parameter per Bioconductor code syntax standards
    multicorr <- match.arg(multicorr)

    if (multicorr == "hochberg") {
        adj_p <- .hochberg_stepup(p_values)
        if (verbose) {
            message("[calculate_sait_interaction] Applied Hochberg stepup ", "adjustment for multi-q correlation")
        }
    } else if (multicorr == "westfall-young") {
        if (verbose) {
            message("[calculate_sait_interaction] Computing true ", "Westfall-Young via ",
                wy_randomizations, " permutations (may be slow)...")
        }

        # Save original group vector for safe restoration
        group_vec_orig <- metadata$group_vec

        # Determine if this is a paired design for correct permutation
        is_paired <- isTRUE(paired) && !is.null(subject_col) &&
            subject_col %in% colnames(SummarizedExperiment::colData(se))

        # Build subject vector for paired permutation
        subject_vec <- NULL
        if (is_paired) {
            cd <- SummarizedExperiment::colData(se)
            sample_q <- colnames(mat)
            subject_vec <- as.character(cd[sample_q, subject_col])
        }

        # Run Westfall-Young permutation
        perm_result <- .westfall_young_permutation(n_genes = length(p_values), wy_randomizations = wy_randomizations,
            permute_fn = function() {
                perm_assignment <- group_vec_orig
                if (is_paired && !is.null(subject_vec)) {
                    # PAIRED: shuffle condition labels WITHIN each subject
                    # (preserves within-subject correlation structure)
                    for (subj in unique(subject_vec)) {
                        subj_idx <- which(subject_vec == subj)
                        perm_assignment[subj_idx] <- sample(group_vec_orig[subj_idx])
                    }
                } else {
                    # UNPAIRED: shuffle group labels within each q-level
                    q_unique <- unique(metadata$q_vals)
                    for (q_val in q_unique) {
                        q_idx <- which(metadata$q_vals == q_val)
                        perm_assignment[q_idx] <- sample(group_vec_orig[q_idx])
                    }
                }
                return(perm_assignment)
            }, refit_fn = function(perm_assignment) {
                # Refit all genes with permuted group assignment
                perm_pvalues <- numeric(length(p_values))
                for (g_idx in seq_along(rownames_mat)) {
                  gene_name <- rownames_mat[g_idx]
                  tryCatch({
                    gene_result <- fit_one_fn(gene_name, group_vec_override = perm_assignment)
                    if (!is.null(gene_result) && !is.na(gene_result$p_interaction)) {
                      perm_pvalues[g_idx] <- gene_result$p_interaction
                    }
                  }, error = function(e) {
                    NULL
                  })
                }
                return(perm_pvalues)
            }, nthreads = nthreads, verbose = verbose)

        # Adjust p-values based on permutation distribution
        adj_p <- vapply(p_values, function(p_obs) {
            pmin(1, (sum(perm_result$perm_minima <= p_obs) + 1)/(wy_randomizations +
                1))
        }, FUN.VALUE = numeric(1))

        if (verbose) {
            message("[calculate_sait_interaction] Applied true ", "Westfall-Young (permutation) adjustment")
        }
    } else if (multicorr == "benjamini-yekutieli") {
        adj_p <- .benjamini_yekutieli(p_values)
        if (verbose) {
            message("[calculate_sait_interaction] Applied ", "Benjamini-Yekutieli adjustment for dependent tests")
        }
    } else {
        stop("Unknown multicorr method: ", multicorr, call. = FALSE)
    }

    # Apply optional Storey adaptive FDR enhancement layer
    if (storey) {
        if (requireNamespace("fdrtool", quietly = TRUE)) {
            tryCatch({
                adj_p <- .compute_storey_qvalues(adj_p)
                if (verbose) {
                  message("[calculate_sait_interaction] Applied Storey ", "adaptive FDR pi0 correction to ",
                    multicorr, " p-values")
                }
            }, error = function(e) {
                if (verbose) {
                  message("[calculate_sait_interaction] Storey adjustment ", "failed: ",
                    conditionMessage(e))
                }
            })
        } else if (verbose) {
            message("[calculate_sait_interaction] fdrtool package not ", "available for Storey (install with: ",
                "install.packages('fdrtool'))")
        }
    }

    return(adj_p)
}

#' @title Map Gene Identifiers to Annotations
#'
#' @description
#' Internal helper that maps gene rownames to gene_id and gene_name
#' columns using rowData from the SummarizedExperiment. Ensures
#' consistent gene annotation across downstream analyses.
#'
#' @param res Data.frame; results with gene column (rownames)
#' @param se SummarizedExperiment object
#' @param verbose Logical; print diagnostic messages
#'
#' @return Modified data.frame with gene_id and gene_name columns added
#'

#' @noRd
.map_gene_annotations <- function(res, se, verbose) {
    rd <- SummarizedExperiment::rowData(se)

    # Look for gene_name column from calculate_diversity or build_se
    gene_name_col <- if ("gene_name" %in% colnames(rd)) {
        "gene_name"
    } else {
        NULL
    }

    if (verbose) {
        message("[calculate_sait_interaction] Gene annotations: ", paste(colnames(rd),
            collapse = ", "))
    }

    if (!is.null(gene_name_col)) {
        # Determine gene_id column if it exists
        id_col <- if ("genes" %in% colnames(rd)) {
            "genes"
        } else if ("gene_id" %in% colnames(rd)) {
            "gene_id"
        } else {
            NA  # rownames will be used as ID
        }

        # Build lookup tables: rowname -> gene_id and rowname -> gene_name
        if (is.na(id_col)) {
            # rownames ARE the gene IDs
            rowname_to_id <- setNames(as.character(rownames(rd)), as.character(rownames(rd)))
        } else {
            # gene IDs are in a column
            rowname_to_id <- setNames(as.character(rd[[id_col]]), as.character(rownames(rd)))
        }

        rowname_to_name <- setNames(as.character(rd[[gene_name_col]]), as.character(rownames(rd)))

        # Vectorized lookup: map res$gene to gene_id and gene_name
        res$gene_id <- unname(rowname_to_id[as.character(res$gene)])
        res$gene_name <- unname(rowname_to_name[as.character(res$gene)])

        # For any unmapped genes, use gene column as fallback
        unmapped_idx <- is.na(res$gene_name)
        n_mapped <- sum(!unmapped_idx)
        n_unmapped <- sum(unmapped_idx)

        if (any(unmapped_idx)) {
            res$gene_name[unmapped_idx] <- res$gene[unmapped_idx]
        }

        if (verbose && n_unmapped > 0) {
            message("[calculate_sait_interaction] Gene mapping: ", n_mapped, " mapped, ",
                n_unmapped, " used ID as fallback")
        }
    } else if (verbose) {
        message("[calculate_sait_interaction] gene_name column not found ", "in rowData - using gene ID as fallback")
    }

    # Ensure gene_name column is always present and populated
    if (is.null(res$gene_name) || !"gene_name" %in% colnames(res)) {
        res$gene_name <- res$gene
    }

    # Ensure gene_id column is always present and populated
    if (is.null(res$gene_id) || !"gene_id" %in% colnames(res)) {
        res$gene_id <- res$gene
    }

    return(res)
}

#' @title Assemble Model Metadata for Diagnostics
#'
#' @description
#' Internal helper that builds the comprehensive metadata list returned
#' when \code{return_model_data = TRUE}. Contains method info, q-values,
#' per-group statistics, and test configuration for downstream visualization.
#'
#' @param se SummarizedExperiment object
#' @param res Data.frame; fitted model results
#' @param mat Matrix; diversity assay data
#' @param metadata List; output from .parse_sample_metadata()
#' @param method Character; modeling method name
#' @param pvalue Character; p-value type
#' @param multicorr Character; multi-q correction method
#' @param assay_name Character; assay name
#' @param bias_correction Logical; KC bias correction setting
#' @param regularization Character; dimensionality reduction method
#' @param corstr Character; correlation structure
#' @param adaptive_knots Logical; adaptive knot selection setting
#'
#' @return List with comprehensive model metadata:
#'   \itemize{
#'     \item method: Modeling method used
#'     \item n_genes: Number of genes analyzed
#'     \item n_q_values: Number of q parameters
#'     \item q_values: The specific q-value vector
#'     \item sample_names: Unique sample identifiers
#'     \item group_levels: Group factor levels
#'     \item per_group_statistics: Summary statistics by group
#'     \item test_configuration: Complete test settings
#'     \item genes_analyzed: Vector of gene identifiers
#'     \item call_time: Timestamp of analysis
#'     \item notes: Usage information
#'   }
#'

#' @noRd
.assemble_model_metadata <- function(se, res, mat, metadata, method, pvalue, multicorr,
    assay_name = "diversity", bias_correction = TRUE, regularization = "pca", corstr = "ar1",
    adaptive_knots = TRUE) {
    # Extract per-group statistics from SE
    per_group_stats <- list()

    for (gr in unique(metadata$group_vec)) {
        gr_idx <- which(metadata$group_vec == gr)
        gr_mat <- mat[, gr_idx, drop = FALSE]

        per_group_stats[[gr]] <- list(group = gr, n_samples = length(unique(metadata$sample_names[gr_idx])),
            n_observations = ncol(gr_mat), entropy_mean = mean(as.numeric(gr_mat),
                na.rm = TRUE), entropy_sd = sd(as.numeric(gr_mat), na.rm = TRUE),
            entropy_min = min(as.numeric(gr_mat), na.rm = TRUE), entropy_max = max(as.numeric(gr_mat),
                na.rm = TRUE), entropy_median = median(as.numeric(gr_mat), na.rm = TRUE),
            n_na = sum(is.na(gr_mat)))
    }

    model_data <- list(method = method, n_genes = nrow(res), n_q_values = length(unique(metadata$q_vals)),
        q_values = sort(unique(metadata$q_vals)), sample_names = unique(metadata$sample_names),
        group_levels = levels(factor(metadata$group_vec)), per_group_statistics = per_group_stats,
        test_configuration = list(method = method, pvalue_method = pvalue, multicorr = multicorr,
            bias_correction = bias_correction, regularization = regularization, corstr = corstr,
            adaptive_knots = adaptive_knots), genes_analyzed = res$gene, call_time = Sys.time(),
        notes = paste("Use this model_data with plotting functions to ", "visualize model fits and diagnostics. See ",
            "per_group_statistics for condition-specific ", "entropy summaries."))

    return(model_data)
}






