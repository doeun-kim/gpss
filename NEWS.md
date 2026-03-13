# gpss 1.0.0

* Initial CRAN release.
* Core GP functions: `gp_train()`, `gp_predict()`, and formula interface `gpss()` with `predict()` method.
* Regression discontinuity design via `gp_rdd()` and `gp_rdd_plot()`.
* Interrupted time-series design via `gp_its()` with optional placebo checks.
* Multiple kernel types: Gaussian, Gaussian-linear, Gaussian-quadratic, Gaussian-periodic-linear, and Gaussian-periodic-quadratic.
* Automatic bandwidth selection via kernel variance maximization.
* Optional noise variance optimization via marginal likelihood.
* Support for mixed continuous/categorical covariates.
* Prior mean support for incorporating external predictions.
* Companion paper: Cho, Kim, and Hazlett (2026), "Inference at the data's edge: Gaussian processes for modeling and inference under model-dependency, poor overlap, and extrapolation," *Political Analysis*, doi:10.1017/pan.2026.10032.
