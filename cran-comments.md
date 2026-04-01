## R CMD check results

0 errors | 0 warnings | 1 note

* This is a resubmission. Changes from v1.0.2:
  - Removed incorrect single quotes around non-package names in DESCRIPTION
  - Added () after function names in DESCRIPTION text
  - Added \value documentation to all exported S3 methods
    (predict.gpss, print.gp_its, print.gp_rdd, summary.gp_its, summary.gp_rdd)
* Changes from earlier versions (carried forward):
  - Added src/Makevars to link LAPACK/BLAS on Unix (fixes installation
    failure on Debian reported in initial submission)
  - Added `notes/` to .Rbuildignore to remove non-standard top-level
    directory from the built package

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Notes

This is the companion R package for:

Cho, Kim, and Hazlett (2026). "Inference at the data's edge: Gaussian
processes for modeling and inference under model-dependency, poor overlap,
and extrapolation." *Political Analysis*, doi:10.1017/pan.2026.10032.
