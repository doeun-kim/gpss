# Data from National Supported Work program and Panel Study in Income Dynamics

Dehejia and Wahba (1999) sample of data from Lalonde (1986). This data
set includes 185 treated units from the National Supported Work (NSW)
program, paired with 2490 control units drawn from the Panel Study of
Income Dynamics (PSID-1).

The treatment variable of interest is `nsw`, which indicates that an
individual was in the job training program. The main outcome of interest
is real earnings in 1978 (`re78`). The remaining variables are
characteristics of the individuals, to be used as controls.

## Usage

``` r
lalonde
```

## Format

A data frame with 2675 rows and 14 columns.

- nsw:

  treatment indicator: participation in the National Supported Work
  program.

- re78:

  real earnings in 1978 (outcome)

- u78:

  unemployed in 1978; actually an indicator for zero income in 1978

- age:

  age in years

- black:

  indicator for identifying as black

- hisp:

  indicator for identifying as Hispanic

- race_ethnicity:

  factor for self-identified race/ethnicity; same information as `black`
  and `hisp` in character form.

- married:

  indicator for being married

- re74:

  real income in 1974

- re75:

  real income in 1975

- u74:

  unemployment in 1974; actually an indicator for zero income in 1974

- u75:

  unemployment in 1975; actually an indicator for zero income in 1975

- educ:

  Years of education of the individual

- nodegr:

  indicator for no high school degree; actually an indicator for years
  of education less than 12

## References

Dehejia, Rajeev H., and Sadek Wahba. "Causal effects in non-experimental
studies: Reevaluating the evaluation of training programs." Journal of
the American statistical Association 94.448 (1999): 1053-1062.

LaLonde, Robert J. "Evaluating the econometric evaluations of training
programs with experimental data." The American economic review (1986):
604-620.
