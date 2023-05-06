library(tidyverse)
library(data.table)
library(xgboost)
library(shapr)
library(ggplot2)
library(ggbeeswarm)

# index_given: feature combination, determined by partial order
# Return list of size nrows(index_given)
sample_gaussian <- function(index_given, n_samples, mu, cov_mat, m, x_test) {

  # Check input
  stopifnot(is.matrix(x_test))

  # Handles the unconditional and full conditional separately when predicting
  cnms <- colnames(x_test)
  if (length(index_given) %in% c(0, m)) return(data.table::as.data.table(x_test))

  dependent_ind <- (1:length(mu))[-index_given]
  x_test_gaussian <- x_test[index_given]
  tmp <- condMVNorm::condMVN(
    mean = mu,
    sigma = cov_mat,
    dependent.ind = dependent_ind,
    given.ind = index_given,
    X.given = x_test_gaussian
  )

  # Makes the conditional covariance matrix symmetric in the rare case where numerical instability made it unsymmetric
  if (!isSymmetric(tmp[["condVar"]])) {
    tmp[["condVar"]] <- Matrix::symmpart(tmp$condVar)
  }

  ret0 <- mvnfast::rmvn(n = n_samples, mu = tmp$condMean, sigma = tmp$condVar)

  ret <- matrix(NA, ncol = m, nrow = n_samples)
  ret[, index_given] <- rep(x_test_gaussian, each = n_samples)
  ret[, dependent_ind] <- ret0

  colnames(ret) <- cnms
  return(as.data.table(ret))
}

respects_order <- function(index_given, ordering) {
  
  for (i in index_given) {
    
    idx_position <- Position(function(x) i %in% x, ordering, nomatch = 0)
    
    stopifnot(idx_position > 0) # It should always be in the ordering
    
    # check for precedents (only relevant if not root set)
    if (idx_position > 1) {
      
      # get precedents
      precedents <- unlist(ordering[1:(idx_position-1)])
      # all precedents must be in index
      # print(precedents)
      # print(intersect(precedents, index))
      # print(setequal(precedents, intersect(precedents, index)))
      if (!setequal(precedents, intersect(precedents, index_given))) {
        # print("here")
        return (FALSE)
      }
    }
  }
  
  TRUE
}

explainer_x_test <- function(x_test, feature_labels) {

  # Remove variables that were not used for training
  x <- data.table::as.data.table(x_test)
  cnms_remove <- setdiff(colnames(x), feature_labels)
  if (length(cnms_remove) > 0) x[, (cnms_remove) := NULL]
  data.table::setcolorder(x, feature_labels)

  return(as.matrix(x))
}

# First step of gaussian-based explanation
gaussian_inputs <- function(x_test, explainer, prediction_zero, approach = "gaussian") {
    explainer$x_test <- explainer_x_test(x_test, explainer$feature_labels)
    explainer$approach <- approach
    
    explainer$mu <- unname(colMeans(explainer$x_train))
    cov_mat <- stats::cov(explainer$x_train)
    
    eigen_values <- eigen(cov_mat)$values
      if (any(eigen_values <= 1e-06)) {
        explainer$cov_mat <- as.matrix(Matrix::nearPD(cov_mat)$mat)
      } else {
        explainer$cov_mat <- cov_mat
      }
    return(explainer)
}

# Second step of gaussian-based explanation
# Return sample data for every combination                      
prepare_data_gaussian <- function(explainer, seed = 1, n_samples = 1e3, index_features = NULL, asymmetric = FALSE, ordering = NULL) {

  id <- id_combination <- w <- NULL # due to NSE notes in R CMD check

  n_xtest <- nrow(explainer$x_test)

  dt_l <- list()
  if (!is.null(seed)) set.seed(seed)
  if (is.null(index_features)) {
    features <- explainer$X$features
  } else {
    features <- explainer$X$features[index_features]
  }

  if (asymmetric == TRUE) {

    message("Asymmetric flag enabled. Only using permutations consistent with the ordering.")

    # By default, no ordering in specified, meaning all variables are in one component.
    if (is.null(ordering)) {
      message("Using no ordering by default.")
      ordering <- list(1:ncol(x$x_test))
    }

    # Filter out the features that do not agree with the order
    features <- features[sapply(features, respects_order, ordering = ordering)]
  }

  for (i in seq(n_xtest)) {
    l <- lapply(
      X = features,
      FUN = sample_gaussian,
      n_samples = n_samples,
      mu = explainer$mu,
      cov_mat = explainer$cov_mat,
      m = ncol(explainer$x_test),
      x_test = explainer$x_test[i, , drop = FALSE]
    )

    dt_l[[i]] <- data.table::rbindlist(l, idcol = "id_combination")
    dt_l[[i]][, w := 1 / n_samples]
    dt_l[[i]][, id := i]
    if (!is.null(index_features)) dt_l[[i]][, id_combination := index_features[id_combination]]
  }
  combination_table <- data.table::rbindlist(dt_l, use.names = TRUE, fill = TRUE)
  combination_table[id_combination %in% c(1, 2^ncol(explainer$x_test)), w := 1.0]
  return(combination_table)
}
                             
# Third step of explanation
prediction <- function(dt, prediction_zero, explainer) {

  # Checks on input data
  id <- w <- id_combination <- p_hat <- NULL # due to NSE notes in R CMD check
  stopifnot(
    data.table::is.data.table(dt),
    !is.null(dt[["id"]]),
    !is.null(dt[["id_combination"]]),
    !is.null(dt[["w"]])
  )

  # Setup
  cnms <- colnames(explainer$x_test)
  data.table::setkeyv(dt, c("id", "id_combination"))

  # Check that the number of test observations equals max(id)
  stopifnot(nrow(explainer$x_test) == dt[, max(id)])

  # Predictions
  dt[, p_hat := predict_model(explainer$model, newdata = .SD), .SDcols = cnms]
  dt[id_combination == 1, p_hat := prediction_zero]
  p_all <- predict_model(explainer$model, newdata = explainer$x_test)
  dt[id_combination == max(id_combination), p_hat := p_all[id]]

  # Calculate contributions
  dt_res <- dt[, .(k = sum((p_hat * w) / sum(w))), .(id, id_combination)]
  data.table::setkeyv(dt_res, c("id", "id_combination"))
  dt_mat <- data.table::dcast(dt_res, id_combination ~ id, value.var = "k")
  dt_mat[, id_combination := NULL]
  kshap <-  t(explainer$W %*% as.matrix(dt_mat))
  dt_kshap <- data.table::as.data.table(kshap)
  colnames(dt_kshap) <- c("none", cnms)

  r <- list(dt = dt_kshap, model = explainer$model, p = p_all, x_test = explainer$x_test)
  attr(r, "class") <- c("shapr", "list")

  return(r)
}

make_ordering <- function(ordering, column_names){
  ordering <- read.table(
    text=gsub("(?<=[a-z])\\s+", "\n", ordering, perl=TRUE), 
    header=FALSE,
    col.names = c("income", "race", "zip")
  )
  ordering <- as.data.table(ordering)
  colnames(ordering) <- column_names
  ordering[, ordering_id := 1:nrow(ordering)]
  return(ordering)
}


tag_coalitions <- function(ordering, var_i){
  #TODO: use something like for (col in paste0("V", 20:100))dt[ , (col) := sqrt(dt[[col]])]
  ordering[, incl.income := as.integer(income <= get(var_i))]
  ordering[, incl.race := as.integer(race <= get(var_i))]
  ordering[, incl.zip := as.integer(zip <= get(var_i))]

  ordering[, c("excl.income", "excl.race", "excl.zip") := list(incl.income, incl.race, incl.zip)]
  ordering[, paste0("excl.", var_i) := 0]
  return(ordering)
}
                             
add_asymmetric_weights <- function(ordering, weight){
  proximate_weight <- data.table(
      ordering_id = 1:nrow(ordering),
      weight = weight
  )

  ordering <- proximate_weight[ordering, on = c(ordering_id="ordering_id")]
  return(ordering)
}

make_combination_index <- function(explainer, x_var){
  combination_index <- data.table::as.data.table(explainer$S)
  combination_index[, id_combination := 1:2^ncol(explainer$x_test)]
  # Vectors are immutable
  colnames(combination_index) <- c(x_var, "id_combination")
  return(combination_index)
}

# Add p for orderings that include the target variable
add_combination_index <- function(ordering, combination_index, incl){
  prefix <- if (incl) "incl" else "excl"
  column_names <- colnames(combination_index)
  colnames(combination_index) <- paste("excl", column_names, sep = ".")
  colnames(combination_index) <- paste(prefix, column_names, sep = ".")

  ordering <- (
      combination_index[ordering,
          on = c(
            eval(paste(prefix, "income", sep = ".")),
            eval(paste(prefix, "race", sep = ".")),
            eval(paste(prefix, "zip", sep = "."))
          )
      ]
  )
  
  colnames(combination_index) <- column_names # Useless, col change only in function scope
  return(ordering)
}

# One ordering block per test instance
repeat_ordering <- function(ordering, x_test, order_count){
  repeat_ordering <- rep(seq(1, order_count), nrow(x_test))
  return(ordering[repeat_ordering])
}

# Incremental id per test instance block
add_test_ordering_id <- function(ordering, x_test, order_count){
  test_instance_ids <- sort(rep(seq(1,nrow(x_test)), order_count))
  ordering[, test_instance_id := test_instance_ids]
  return(ordering)
}
                             
initiate_ordering <- function(var_i, asymmetric_weight=c(1/3, 0, 1/3, 0, 1/3, 0)){
  ordering <- "1 2 3 1 3 2 3 1 2 3 2 1 2 1 3 2 3 1"
  ordering <- read.table(text=gsub("(?<=[a-z])\\s+", "\n", ordering, perl=TRUE), 
              header=FALSE, col.names = c("income", "race", "zip"))
  ordering <- as.data.table(ordering)
  colnames(ordering) <- x_var
  ordering[, ordering_id := 1:nrow(ordering)]

  #TODO: use something like for (col in paste0("V", 20:100))dt[ , (col) := sqrt(dt[[col]])]
  ordering[, incl.income := as.integer(income <= get(var_i))]
  ordering[, incl.race := as.integer(race <= get(var_i))]
  ordering[, incl.zip := as.integer(zip <= get(var_i))]

  ordering[, c("excl.income", "excl.race", "excl.zip") := list(incl.income, incl.race, incl.zip)]
  ordering[, paste0("excl.", var_i) := 0]

  proximate_weight <- data.table(
      ordering_id = 1:nrow(ordering),
      weight = asymmetric_weight
  )

  ordering <- proximate_weight[ordering, on = c(ordering_id="ordering_id")]
  ordering <- ordering[weight>0, ]
  return(ordering)
}

# Add combination id for incl and excl combinations
add_combination_index <- function(ordering, explainer, x_var){
  combination_index <- copy(explainer$S)
  combination_index <- data.table::as.data.table(combination_index)
  combination_index[, id_combination := 1:2^ncol(explainer$x_test)]
  column_names <- copy(x_var)
  column_names <- append(column_names, "id_combination")
  colnames(combination_index) <- column_names

  # Add p for orderings that include the target variable
  colnames(combination_index) <- paste("incl", column_names, sep = ".")

  ordering <- (
      combination_index[ordering,
          on = c("incl.income", "incl.race", "incl.zip")
      ]
  )

  # Add p for orderings that include the target variable
  colnames(combination_index) <- paste("excl", column_names, sep = ".")

  ordering <- (
      combination_index[ordering,
          on = c("excl.income", "excl.race", "excl.zip")
      ]
  )
  
  return(ordering)
}
                             
add_p <- function(ordering, dt_res, incl){
  prefix <- if (incl) "incl" else "excl"
  prefix_test_instance_id <- paste(prefix, "test_instance_id", sep = ".")
  prefix_id_combination <- paste(prefix, "id_combination", sep = ".")
  colnames(ordering)[colnames(ordering) == "test_instance_id"] <- prefix_test_instance_id
  colnames(dt_res) <- paste(prefix, colnames(dt_res), sep = ".")
  
  ordering <- dt_res[
    ordering,
    on=c(prefix_id_combination, prefix_test_instance_id)
  ]
  
  colnames(ordering)[colnames(ordering)==prefix_test_instance_id] <- "test_instance_id"
  return(ordering)
}

# Repeat ordering for every test case, adding ordering and test id's
add_test_ids <- function(ordering, explainer){
  order_count <- nrow(ordering)

  repeat_orderings <- rep(seq(1, order_count), nrow(explainer$x_test))

  ordering <- ordering[repeat_orderings]

  test_instance_ids <- sort(rep(seq(1, nrow(explainer$x_test)), order_count))
  
  ordering[, test_instance_id := test_instance_ids]
  return(ordering)
}
                       
add_coalition_prediction <- function(ordering, dt_res){
  dt_res_column_names <- c("test_instance_id", "id_combination", "p")
  colnames(dt_res) <- dt_res_column_names

  ordering <- add_p(ordering, dt_res, incl=TRUE)
  ordering <- add_p(ordering, dt_res, incl=FALSE)
  return(ordering)
}
                        
get_ordering <- function(var_i, explainer, x_var, dt_res){
  ordering <- initiate_ordering(var_i)
  order_count <- nrow(ordering)

  ordering <- add_combination_index(ordering, explainer, x_var)
  
  ordering <- add_test_ids(ordering, explainer)

  ordering <- add_coalition_prediction(ordering, dt_res)

  ordering[, v_difference := incl.p - excl.p]

  ordering[, weighted_v_difference := weight*v_difference]

  return(ordering)
}

