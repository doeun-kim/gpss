#' Gaussian Process Interrupted Time Series
#'
#' @param y numeric vector or ts object of outcome values
#' @param dates Date vector of observation dates
#' @param date_treat Date of treatment/intervention
#' @param covariates optional matrix/data.frame of covariates
#' @param kernel_type kernel type for GP
#' @param b bandwidth (NULL for automatic)
#' @param s2 noise variance (default 0.3)
#' @param period period for periodic kernels
#' @param auto_period logical; auto-detect period from pre-treatment data
#' @param time_col character or numeric; which column is the time variable
#' @param scale logical; scale covariates
#' @param optimize logical; optimize s2 via MLE
#' @param interval_type "prediction" or "confidence"
#' @param alpha significance level
#' @param mixed_data logical; whether covariates contain categorical variables
#' @param cat_columns character vector of categorical column names
#' @param placebo_check logical; run placebo checks
#' @param placebo_periods number of placebo periods (NULL for automatic)
#' @param verbose logical; print progress messages
#'
#' @return An object of class \code{"gp_its"}, a list containing:
#' \item{y}{the outcome vector}
#' \item{dates}{the date vector}
#' \item{date_treat}{the treatment date}
#' \item{y0_hat}{counterfactual predictions for post-treatment periods}
#' \item{estimates}{data frame with columns \code{tau_t}, \code{tau_cum}, \code{tau_avg}}
#' \item{se}{data frame of standard errors}
#' \item{ci}{data frame of confidence interval bounds}
#' \item{gp_model}{the fitted GP model from \code{\link{gp_train}}}
#' \item{placebo_estimates}{placebo check results (if requested)}
#'
#' @examples
#' \donttest{
#' # Simulated interrupted time series
#' set.seed(42)
#' n <- 50
#' dates <- seq(as.Date("2020-01-01"), by = "month", length.out = n)
#' y <- rnorm(n) + c(rep(0, 30), rep(2, 20))
#' res <- gp_its(y, dates, date_treat = as.Date("2022-07-01"))
#' print(res)
#' }
#' @export
gp_its <- function(y, dates, date_treat,
                   covariates = NULL,
                   kernel_type = "gaussian",
                   b = NULL,
                   s2 = 0.3,
                   period = NULL,
                   auto_period = FALSE,
                   time_col = NULL,
                   scale = TRUE,
                   optimize = TRUE,
                   interval_type = c("prediction", "confidence"),
                   alpha = 0.05,
                   mixed_data = FALSE,
                   cat_columns = NULL,
                   placebo_check = FALSE,
                   placebo_periods = NULL,
                   verbose = FALSE) {

  # Input validation
  if (!is.numeric(y) && !inherits(y, "ts")) stop("`y` must be numeric or ts object.")
  if (any(is.na(y))) stop("`y` contains NA values.")
  if (inherits(y, "ts")) y <- as.numeric(y)

  n <- length(y)
  if (length(dates) != n) stop("`dates` and `y` must have same length.")
  if (!inherits(dates, "Date")) stop("`dates` must be Date vector.")
  if (length(date_treat) != 1 || !inherits(date_treat, "Date")) stop("`date_treat` must be single Date.")
  if (date_treat <= min(dates) || date_treat > max(dates)) stop("`date_treat` must be within range of `dates`.")
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be between 0 and 1.")

  interval_type <- match.arg(interval_type)
  kernel_type <- match.arg(kernel_type,
                           c("gaussian", "gaussian_linear", "gaussian_quadratic",
                             "gaussian_periodic_linear", "gaussian_periodic_quadratic"))

  if (!is.null(placebo_periods) && !placebo_check) {
    placebo_check <- TRUE
    if (verbose) message("Enabling placebo check (placebo_periods specified)")
  }

  # Pre/post split
  pre_idx <- dates < date_treat
  n_pre <- sum(pre_idx)
  n_post <- n - n_pre
  if (n_post < 1) stop("Need at least 1 post-treatment period.")

  # Period handling
  if (grepl("periodic", kernel_type) && is.null(period)) {
    if (auto_period) {
      period <- estimate_period(y[pre_idx])
      if (verbose) message(sprintf("Auto-detected period: %.2f", period))
    } else {
      stop(sprintf("Period required for kernel '%s'. Set `period` or `auto_period = TRUE`.", kernel_type))
    }
  }

  # Covariates
  if (!is.null(covariates)) {
    if (nrow(covariates) != n) stop("`covariates` must have same rows as `y`.")
    if (anyNA(covariates)) stop("`covariates` contains NA values.")
    covariates <- as.matrix(covariates)
    if (is.null(colnames(covariates))) {
      colnames(covariates) <- paste0("cov", seq_len(ncol(covariates)))
    }
  }

  # Design matrix
  time_id <- seq_len(n)
  X_full <- if (is.null(covariates)) {
    matrix(time_id, ncol = 1, dimnames = list(NULL, "time"))
  } else {
    cbind(time = time_id, covariates)
  }

  # Time indices (needed for plot)
  time_id_pre <- which(pre_idx)
  time_id_post <- which(!pre_idx)
  time_id_treat <- time_id_post[1]

  # Fit GP model
  if (verbose) message("Fitting GP-ITS model...")

  gp_model <- gp_train(
    X = X_full[pre_idx, , drop = FALSE],
    Y = y[pre_idx],
    b = b, s2 = s2, optimize = optimize, scale = scale,
    kernel_type = kernel_type, period = period, time_col = time_col,
    mixed_data = mixed_data, cat_columns = cat_columns
  )

  # Predictions
  gp_pre <- gp_predict(gp_model, X_full[pre_idx, , drop = FALSE])
  gp_post <- gp_predict(gp_model, X_full[!pre_idx, , drop = FALSE])

  # Select covariance based on interval type
  pre_cov <- if (interval_type == "prediction") gp_pre$Ys_cov_orig else gp_pre$f_cov_orig
  post_cov <- if (interval_type == "prediction") gp_post$Ys_cov_orig else gp_post$f_cov_orig

  z_crit <- qnorm(1 - alpha / 2)
  y_pre_se <- sqrt(pmax(diag(pre_cov), 0))
  y0_se <- sqrt(pmax(diag(post_cov), 0))

  # Pre-treatment fit
  y_pre_fitted <- gp_pre$Ys_mean_orig
  y_pre_lwr <- y_pre_fitted - z_crit * y_pre_se
  y_pre_upr <- y_pre_fitted + z_crit * y_pre_se

  # Counterfactual predictions
  y0_hat <- gp_post$Ys_mean_orig
  y0_hat_lwr <- y0_hat - z_crit * y0_se
  y0_hat_upr <- y0_hat + z_crit * y0_se

  # Treatment effects
  y_post <- y[!pre_idx]
  tau_t <- as.numeric(y_post - y0_hat)
  tau_cum <- cumsum(tau_t)
  tau_avg <- tau_cum / seq_len(n_post)

  # Standard errors
  tau_t_se <- y0_se
  tau_cum_var <- vapply(seq_len(n_post), function(i) sum(post_cov[1:i, 1:i]), numeric(1))
  tau_cum_se <- sqrt(pmax(tau_cum_var, 0))
  tau_avg_se <- tau_cum_se / seq_len(n_post)

  # Placebo checks
  placebo_estimates <- NULL
  placebo_summary <- NULL

  if (placebo_check) {
    if (verbose) message("Running placebo checks...")
    placebo_res <- .run_placebo_checks(
      y = y, X_full = X_full, n_pre = n_pre,
      kernel_type = kernel_type, period = gp_model$period_original, time_col = time_col,
      scale = scale, interval_type = interval_type, alpha = alpha,
      optimize = optimize, mixed_data = mixed_data, cat_columns = cat_columns,
      placebo_periods = placebo_periods, verbose = verbose
    )
    if (!is.null(placebo_res)) {
      placebo_estimates <- placebo_res$estimates
      placebo_summary <- placebo_res$summary
    }
  }

  # Assemble output
  res <- list(
    # Core data
    y = as.numeric(y),
    dates = dates,
    date_treat = date_treat,
    dates_pre = dates[pre_idx],
    dates_post = dates[!pre_idx],
    # Time indices (required by plot)
    time_id = time_id,
    time_id_pre = time_id_pre,
    time_id_post = time_id_post,
    time_id_treat = time_id_treat,
    time_id_placebo = if (!is.null(placebo_estimates)) placebo_estimates$time_id else NULL,
    # Pre-treatment fit
    y_pre_fitted = y_pre_fitted,
    y_pre_lwr = y_pre_lwr,
    y_pre_upr = y_pre_upr,
    # Counterfactual
    y0_hat = y0_hat,
    y0_hat_lwr = y0_hat_lwr,
    y0_hat_upr = y0_hat_upr,
    # Effects as data.frames
    estimates = data.frame(tau_t = tau_t, tau_cum = tau_cum, tau_avg = tau_avg),
    se = data.frame(tau_t_se = tau_t_se, tau_cum_se = tau_cum_se, tau_avg_se = tau_avg_se),
    ci = data.frame(
      tau_t_lwr = tau_t - z_crit * tau_t_se, tau_t_upr = tau_t + z_crit * tau_t_se,
      tau_cum_lwr = tau_cum - z_crit * tau_cum_se, tau_cum_upr = tau_cum + z_crit * tau_cum_se,
      tau_avg_lwr = tau_avg - z_crit * tau_avg_se, tau_avg_upr = tau_avg + z_crit * tau_avg_se
    ),
    post_cov = post_cov,
    # Placebo
    placebo_estimates = placebo_estimates,
    placebo_summary = placebo_summary,
    # Model info
    kernel_type = kernel_type,
    gp_model = gp_model,
    alpha = alpha,
    interval_type = interval_type,
    n_pre = n_pre,
    n_post = n_post
  )

  class(res) <- "gp_its"
  res
}

#' Internal function for placebo checks
#' @keywords internal
.run_placebo_checks <- function(y, X_full, n_pre,
                                kernel_type, period, time_col,
                                scale, interval_type, alpha,
                                optimize, mixed_data, cat_columns,
                                placebo_periods, verbose) {

  n_post <- length(y) - n_pre
  z_crit <- qnorm(1 - alpha / 2)

  # Determine placebo periods
  if (is.null(placebo_periods)) {
    placebo_periods <- min(n_post, floor(n_pre / 2))
    if (verbose) message(sprintf("Using %d placebo periods", placebo_periods))
  }

  # Check minimum training data
  min_training <- 5
  if (n_pre - placebo_periods < min_training) {
    placebo_periods <- max(1, n_pre - min_training)
    if (placebo_periods < 1) {
      warning("Insufficient pre-treatment data for placebo checks.")
      return(NULL)
    }
    warning(sprintf("Reduced placebo periods to %d due to limited pre-treatment data", placebo_periods))
  }

  fake_treat_idx <- n_pre - placebo_periods + 1
  placebo_test_idx <- fake_treat_idx:n_pre
  n_tests <- length(placebo_test_idx)

  # Pre-allocate results
  results <- data.frame(
    time_id = integer(n_tests),
    relative_time = integer(n_tests),
    tau = numeric(n_tests),
    se = numeric(n_tests),
    ci_lwr = numeric(n_tests),
    ci_upr = numeric(n_tests),
    z_score = numeric(n_tests),
    cover = logical(n_tests),
    y_actual = numeric(n_tests),
    y_predicted = numeric(n_tests)
  )

  current_train_idx <- 1:(fake_treat_idx - 1)
  valid_tests <- rep(TRUE, n_tests)

  for (i in seq_len(n_tests)) {
    test_point <- placebo_test_idx[i]
    y_actual <- y[test_point]

    tryCatch({
      gp_fit <- gp_train(
        X = X_full[current_train_idx, , drop = FALSE],
        Y = y[current_train_idx],
        optimize = optimize, scale = scale,
        kernel_type = kernel_type, period = period, time_col = time_col,
        mixed_data = mixed_data, cat_columns = cat_columns,
        Xtest = X_full[test_point, , drop = FALSE]
      )

      gp_pred <- gp_predict(gp_fit, X_full[test_point, , drop = FALSE])
      y0_hat <- gp_pred$Ys_mean_orig
      y0_se <- sqrt(if (interval_type == "prediction")
        gp_pred$Ys_cov_orig[1, 1] else gp_pred$f_cov_orig[1, 1])

      tau <- y_actual - y0_hat
      z <- tau / y0_se

      results[i, ] <- list(
        time_id = test_point,
        relative_time = i,
        tau = tau,
        se = y0_se,
        ci_lwr = tau - z_crit * y0_se,
        ci_upr = tau + z_crit * y0_se,
        z_score = z,
        cover = abs(z) <= z_crit,
        y_actual = y_actual,
        y_predicted = y0_hat
      )
    }, error = function(e) {
      warning(sprintf("Placebo test failed at time %d: %s", test_point, e$message))
      valid_tests[i] <<- FALSE
    })

    current_train_idx <- c(current_train_idx, test_point)
  }

  results <- results[valid_tests, , drop = FALSE]
  if (nrow(results) == 0) {
    warning("All placebo tests failed.")
    return(NULL)
  }

  # Summary statistics
  summary_stats <- list(
    n_tests = nrow(results),
    avg_effect = mean(results$tau),
    rmse = sqrt(mean(results$tau^2)),
    covered = sum(results$cover),
    expected_coverage = 1 - alpha,
    mean_z_score = mean(results$z_score),
    sd_z_score = sd(results$z_score),
    fake_treat_idx = fake_treat_idx
  )

  if (verbose) {
    cat("\n", strrep("=", 60), "\n", sep = "")
    cat("Placebo Check Summary\n")
    cat(strrep("=", 60), "\n", sep = "")
    cat(sprintf("Tests performed: %d\n", summary_stats$n_tests))
    cat(sprintf("%d CIs out of %d placebo periods contain 0\n",
                summary_stats$covered, summary_stats$n_tests))
    cat(sprintf("Mean standardized effect: %.3f\n", summary_stats$mean_z_score))
    cat(sprintf("SD of standardized effects: %.3f\n", summary_stats$sd_z_score))
    cat(strrep("=", 60), "\n", sep = "")
  }

  list(estimates = results, summary = summary_stats)
}

#' Print method for gp_its objects
#' @param x a gp_its object
#' @param ... additional arguments (not used)
#' @keywords internal
#' @export
print.gp_its <- function(x, ...) {
  cat("GP-ITS Model\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("Observations: %d (pre: %d, post: %d)\n", length(x$y), x$n_pre, x$n_post))
  cat(sprintf("Treatment date: %s\n", x$date_treat))
  cat(sprintf("Kernel: %s\n", x$kernel_type))
  if (!is.null(x$gp_model$period_original)) {
    cat(sprintf("Period: %.2f\n", x$gp_model$period_original))
  }
  cat(sprintf("Interval type: %s (%.0f%%)\n", x$interval_type, (1 - x$alpha) * 100))

  cat("\nAverage treatment effect (final period):\n")
  i <- x$n_post
  cat(sprintf("  tau_avg = %.3f (SE = %.3f)\n", x$estimates$tau_avg[i], x$se$tau_avg_se[i]))
  cat(sprintf("  %.0f%% CI: [%.3f, %.3f]\n",
              (1 - x$alpha) * 100, x$ci$tau_avg_lwr[i], x$ci$tau_avg_upr[i]))

  if (!is.null(x$placebo_summary)) {
    cat(sprintf("\nPlacebo check: %d tests, coverage = %.1f%%\n",
                x$placebo_summary$n_tests,
                x$placebo_summary$covered / x$placebo_summary$n_tests * 100))
  }

  invisible(x)
}

#' Summary method for gp_its objects
#' @param object a gp_its object
#' @param ... additional arguments (not used)
#' @keywords internal
#' @export
summary.gp_its <- function(object, ...) {
  x <- object

  cat("GP-ITS Summary\n")
  cat(strrep("=", 50), "\n\n")

  cat("DATA\n")
  cat(sprintf("  Total observations: %d\n", length(x$y)))
  cat(sprintf("  Pre-treatment: %d\n", x$n_pre))
  cat(sprintf("  Post-treatment: %d\n", x$n_post))
  cat(sprintf("  Treatment date: %s\n", x$date_treat))

  cat("\nMODEL\n")
  cat(sprintf("  Kernel: %s\n", x$kernel_type))
  cat(sprintf("  Bandwidth (b): %.4f\n", x$gp_model$b))
  cat(sprintf("  Noise variance (s2): %.4f\n", x$gp_model$s2))
  if (!is.null(x$gp_model$period_original)) {
    cat(sprintf("  Period: %.2f\n", x$gp_model$period_original))
  }

  cat("\nTREATMENT EFFECTS\n")
  cat(sprintf("  Interval type: %s (%.0f%% CI)\n", x$interval_type, (1 - x$alpha) * 100))

  # Summary table
  cat("\n  Period-specific effects:\n")
  effects_df <- data.frame(
    Date = x$dates_post,
    tau_t = round(x$estimates$tau_t, 3),
    SE = round(x$se$tau_t_se, 3),
    CI_lower = round(x$ci$tau_t_lwr, 3),
    CI_upper = round(x$ci$tau_t_upr, 3)
  )
  print(effects_df, row.names = FALSE)

  i <- x$n_post
  cat(sprintf("\n  Cumulative: %.3f [%.3f, %.3f]\n",
              x$estimates$tau_cum[i], x$ci$tau_cum_lwr[i], x$ci$tau_cum_upr[i]))
  cat(sprintf("  Average: %.3f [%.3f, %.3f]\n",
              x$estimates$tau_avg[i], x$ci$tau_avg_lwr[i], x$ci$tau_avg_upr[i]))

  if (!is.null(x$placebo_summary)) {
    cat("\nPLACEBO CHECK\n")
    cat(sprintf("  Tests: %d\n", x$placebo_summary$n_tests))
    cat(sprintf("  Coverage: %d/%d (%.1f%%)\n",
                x$placebo_summary$covered, x$placebo_summary$n_tests,
                x$placebo_summary$covered / x$placebo_summary$n_tests * 100))
    cat(sprintf("  Mean z-score: %.3f (SD = %.3f)\n",
                x$placebo_summary$mean_z_score, x$placebo_summary$sd_z_score))
  }

  invisible(x)
}
