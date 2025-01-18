sanity_formula_data <- function(formula, data) {
  if (!inherits(data, "data.frame")) {
    msg <- "`data` must be a data frame."
    stop(msg, call. = FALSE)
  }

  if (!inherits(formula, "formula")) {
    msg <- "`formula` must be a formula."
    stop(msg, call. = FALSE)
  }

  # Categorical variables: We need to all levels to ensure the Xtrain and Xtest
  # have the same dimensions after 1-hot encoding. This can be difficult to track
  # categorical variables encoded as character vectors, and for variables encoded
  # as factors directly in the formula. For safety, we return an early error to
  # force users to encode their categorical as factors. This allows us to easily
  # track levels, and is informative for the user.

  # no factor() in formula
  vars <- all.vars(formula)
  for (v in vars) {
    if (v %in% colnames(data)) {
    }
  }

  # no character columns
  if (any(grepl("^factor\\(", as.character(formula)))) {
    msg <- "The `gpss()` function does not support `factor()` in the formula. Please convert your variables to factor in the data frame before supplying it to the `data` argument."
    stop(msg, call. = FALSE)
  }
}
