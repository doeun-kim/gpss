# Changelog

## gpss 1.0.0

- Initial CRAN release.
- Core GP functions:
  [`gp_train()`](https://doeunkim.org/gpss/reference/gp_train.md),
  [`gp_predict()`](https://doeunkim.org/gpss/reference/gp_predict.md),
  and formula interface
  [`gpss()`](https://doeunkim.org/gpss/reference/gpss.md) with
  [`predict()`](https://rdrr.io/r/stats/predict.html) method.
- Regression discontinuity design via
  [`gp_rdd()`](https://doeunkim.org/gpss/reference/gp_rdd.md) and
  [`gp_rdd_plot()`](https://doeunkim.org/gpss/reference/gp_rdd_plot.md).
- Interrupted time-series design via
  [`gp_its()`](https://doeunkim.org/gpss/reference/gp_its.md) with
  optional placebo checks.
- Multiple kernel types: Gaussian, Gaussian-linear, Gaussian-quadratic,
  Gaussian-periodic-linear, and Gaussian-periodic-quadratic.
- Automatic bandwidth selection via kernel variance maximization.
- Optional noise variance optimization via marginal likelihood.
- Support for mixed continuous/categorical covariates.
- Prior mean support for incorporating external predictions.
- Companion paper: Cho, Kim, and Hazlett (2026), “Inference at the
  data’s edge: Gaussian processes for modeling and inference under
  model-dependency, poor overlap, and extrapolation,” *Political
  Analysis*, <doi:10.1017/pan.2026.10032>.
