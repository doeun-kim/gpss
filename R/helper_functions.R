### Helper functions

#' one_hot
#'
#' to convert categorical variables into multiple binary variables with 1 and 0.
#'
#' @param data data frame containing only categorical variables
#
#' @return data frame containing expanded form of binary variables
#
#' @importFrom stats model.matrix contrasts
#' @keywords internal
#' @export
one_hot <- function(data) {
  onehot_data <- data.frame(lapply(data.frame(data),as.factor))
  onehot_data <- model.matrix(~ ., onehot_data,
                              contrasts.arg = lapply(onehot_data[,,drop = F],
                                                     contrasts,
                                                     contrasts=FALSE))
  onehot_data <- onehot_data[, -1]

  return(onehot_data)
}


#' mixed_data_processing
#'
#' to convert categorical variables in a data set into multiple binary variables with 1 and 0.
#'
#' @param X data set or matrix
#' @param cat_columns a numerical or character vector that indicates categorical variables
#' @param Xtest (optional) if there is testing data set separate from the training data, include the testing data
#
#' @return data frame containing expanded form of categorical variables
#' @keywords internal
#' @export
mixed_data_processing <- function(X,
                                  cat_columns = NULL,
                                  Xtest=NULL # This is required when Xtest is separately generated
){

  if(is.numeric(cat_columns)){
    cat_num = cat_columns
  }else if(is.character(cat_columns)){
    cat_num = which(colnames(X) %in% cat_columns)
  }

  if(is.null(cat_columns)){
    print("There is no cat_columns defined")
  }else if(is.null(Xtest)){
    allx_cat <- one_hot(X[, cat_num, drop= F])
    allx_cont <- X[, -cat_num, drop = F]
    X <- cbind(allx_cat, allx_cont)
    X <- as.matrix(X)
    cat_num <- 1:ncol(allx_cat)
  }

  if(!is.null(Xtest)){
    Xtest_start <- nrow(X)+1
    Xtmp <- rbind(X, Xtest)
    Xtest_end <- nrow(Xtmp)
    allx_cat <- one_hot(Xtmp[, cat_num, drop= F])
    allx_cont <- Xtmp[, -cat_num, drop = F]
    X <- cbind(allx_cat, allx_cont)
    X <- as.matrix(X)
    X <- X[Xtest_start:Xtest_end, ]

    cat_num <- 1:ncol(allx_cat)
  }

  X = apply(X, 2, as.numeric)

  return(list(X = X,
              cat_num = cat_num))
}


# chol_stable
#
# numerically stable Cholesky decomposition
#
#' @param X data set or matrix
#
#' @return cholesky-decomposed matrix
#' @keywords internal
#' @export
chol_stable <- function(X){
L <- try(chol(X), silent = TRUE)
g <- 1e-6
while ( inherits(L, "try-error") & g<1 ) {
  L <- try(chol(X + diag(g, nrow = nrow(X))), silent = TRUE)
  g <- 10*g
}
if (inherits(L, "try-error")){
  stop("Cholesky decomposition fails as input matrix is not positive semi-definite")
}
return(L)
}

# getb_maxvar
#
#' @param X data set or matrix
#
#' @return cholesky-decomposed matrix
#
#' @importFrom kbal b_maxvarK
#' @keywords internal
#' @export
getb_maxvar <- function(X){
  X = as.matrix(X)
  best_b = kbal::b_maxvarK(data = as.matrix(X), useasbases=rep(1,nrow(X)))
  b_opt = best_b$b_maxvar
  return(b_opt)
}
