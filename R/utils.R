#' Split formula into terms
#' 
#' @param formula Full formula following `fixest` syntax:
#' `y ~ W | W_FE | T ~ Z | Z_FE`.
#' @param parts_as_formula Logical. If `TRUE`, then each part will be a 
#' right-hand side formula. Default is `FALSE`
#' 
#' @return List of expressions/formula for each part of the formula. It will be of type `symbol`/`language` unless `parts_as_formula = TRUE`. Can be used with 
#' `fixest::xpd` and the dot bracket syntax to create formula. Any missing
#' elements will be given a value of `NULL`. The list contains the following:
#' \item{y_fml}{The LHS}
#' \item{W_lin}{The linear part of the exogenous variables}
#' \item{W_FE}{The fixed effects part of the exogenous variables}
#' \item{T_fml}{The endogenous variable}
#' \item{Z_lin}{The linear part of the instruments}
#' \item{Z_FE}{The fixed effects part of the instruments}
#' 
#' @export
get_fml_parts <- function(formula, parts_as_formula = FALSE) {
  has_lhs <- !is_rhs_only(formula)
  fml_split_tilde <- fml_breaker(formula, "~")

  res <- list(
    y_fml = NULL, W_lin = NULL, W_FE = NULL, T_fml = NULL, Z_lin = NULL, Z_FE = NULL
  )

  # LHS
  if (has_lhs) {
    res$y_fml <- fml_split_tilde[[length(fml_split_tilde)]]
    # Drop y and treat everything as RHS only formula
    fml_split_tilde <- fml_split_tilde[-length(fml_split_tilde)]
  }

  # OLS
  if (length(fml_split_tilde) == 1) {
    fml_split <- fml_breaker(fml_split_tilde[[1]], "|")
    if (length(fml_split) == 2) {
      res$W_lin <- fml_split[[2]]
      res$W_FE <- fml_split[[1]]
    } else {
      res$W_lin <- fml_split[[1]]
    }
  }

  # IV
  if (length(fml_split_tilde) == 2) {
    fml_Z_split <- fml_breaker(fml_split_tilde[[1]], "|")
    fml_W_T_split <- fml_breaker(fml_split_tilde[[2]], "|")

    if (length(fml_Z_split) == 2) {
      res$Z_FE <- fml_Z_split[[1]]
      res$Z_lin <- fml_Z_split[[2]]
    } else {
      res$Z_lin <- fml_Z_split[[1]]
    }

    # This works b/c we know there is an IV
    if (length(fml_W_T_split) == 3) {
      res$T_fml <- fml_W_T_split[[1]]
      res$W_FE <- fml_W_T_split[[2]]
      res$W_lin <- fml_W_T_split[[3]]
    } else {
      res$T_fml <- fml_W_T_split[[1]]
      res$W_lin <- fml_W_T_split[[2]]
    }
  }

  if (parts_as_formula) res = lapply(res, rhs_fml_maker)

  return(res)
}

#' Break apart formula (from right to left) based on a symbole (`~` or `|`)
#'
#' @param fml Formula following `fixest` syntax.
#' @param op String. Either `~` or `|`
#' @return list of `symbol` or `language` from right to left that are split at each occurence of `op`.
fml_breaker <- function(fml, op) {
  res <- list()
  k <- 1
  while (is_operator(fml, op) & length(fml) == 3) {
    res[[k]] <- fml[[3]]
    k <- k + 1
    fml <- fml[[2]]
  }

  if (length(fml) == 2) {
    res[[k]] <- fml[[2]]
  } else {
    res[[k]] <- fml
  }
  res
}

# Makes `symbol`/`language` into rhs formula
rhs_fml_maker <- function(rhs) {
  res <- ~.
  res[[2]] <- rhs
  return(res)
}

# Checks if formula is only rhs
is_rhs_only <- function(fml) {
  # e.g. fml = ~ x
  if (length(fml) == 2) {
    return(TRUE)
  }
  # e.g. fml = ~ x | t ~ z
  if (length(fml[[2]]) == 2) {
    return(TRUE)
  }
  return(FALSE)
}

# From fixest
is_operator <- function(x, op) {
  if (length(x) <= 1) {
    return(FALSE)
  } else {
    return(x[[1]] == op)
  }
}

# helper to check arguments from jive
check_args <- function(
    data, formula,
    cluster = NULL, ssc = FALSE) {
  # Check arguments ------------------------------------------------------------
  dreamerr::check_arg(data, "data.frame")
  dreamerr::check_arg(formula, "ts formula var(data)", .data = data)
  if (!is.null(cluster)) {
    dreamerr::check_arg(cluster, "character var(data) | os formula var(data) right(1, 1)", .data = data)
  }
  dreamerr::check_arg(ssc, "logical scalar")
}

#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Tue Sep 22 15:30:09 2020
# ~: Code snippets to share
#----------------------------------------------#


####
#### Sparse model.matrix ####
####

# Here's a quick and dirty (hence limited) code to create a sparse model.matrix from
# a fixest estimation

# What does it do? It returns a list of two elements:
# - mat_RHS: the *sparse* version of the RHS (excluding the FEs)
# - mat_FE: the sparse matrix of the fixed-effects
# NOTE: both are Matrix matrices from library(Matrix)

# Limitations:
# - does not automatically remove references from factors
# - does not handle estimations containing variables with varying slopes
#

# Benefits (yes there are some!):
# - avoids copies of the data in the process of constructing the matrix
# - handles interactions between any number of factors/numeric variables
# - the non-sparse model matrix is never created, hence it's very fast
#

# Be sure to run the *three* functions sparse_model_matrix, mult_wrap and mult_sparse


# x: fixest estimation
sparse_model_matrix = function(x){

    require(Matrix)

    # Linear formula
    fml_lin = formula(x, "lin")

    data = fixest:::fetch_data(x)
    data = as.data.frame(data)
    data = data[obs(x), ]

    #
    # Step 1: Linear matrix
    #

    vars = attr(terms(fml_lin), "term.labels")

    if(length(vars) == 0){
        # Case only FEs
        mat = NULL
    } else {

        # Since we don't want to evaluate the factors,
        # the code is a bit intricate because we have to catch them before
        # any interaction takes place
        #
        # that's why I wrap interactions in a function (mult_sparse())
        #

        # Below, we evaluate all the variables in a "sparse" way

        vars_calls = lapply(vars, mult_wrap)

        n = length(vars)
        variables_list = vector("list", n)
        for(i in 1:n){
            variables_list[[i]] = eval(vars_calls[[i]], data)
        }

        # To create the sparse matrix, we need the indexes

        total_cols = 0
        running_cols = c(0)
        for(i in 1:n){
            xi = variables_list[[i]]
            if(inherits(xi, "sparse_var")){
                total_cols = total_cols + xi$n_cols
            } else {
                total_cols = total_cols + NCOL(xi)
            }
            running_cols[i + 1] = total_cols
        }

        # We just create a sparse matrix and fill it

        # 1) creating the indexes + names

        # NOTA: I use lists to avoid creating copies
        rowid = 1:nrow(data)
        id_all = values_all = names_all = vector("list", n)
        for(i in 1:n){
            xi = variables_list[[i]]
            if(inherits(xi, "sparse_var")){
                id_all[[i]] = cbind(xi$rowid, running_cols[i] + xi$colid)
                values_all[[i]] = xi$values
                names_all[[i]] = paste0(vars[[i]], "::", xi$col_names)
            } else if(NCOL(xi) == 1){
                id_all[[i]] = cbind(rowid, running_cols[i] + 1)
                values_all[[i]] = xi
                names_all[[i]] = vars[[i]]
            } else {
                colid = rep(1:NCOL(xi), each = nrow(data))
                id_all[[i]] = cbind(rep(rowid, NCOL(xi)), running_cols[i] + colid)
                values_all[[i]] = as.vector(xi)
                if(!is.null(colnames(xi))){
                    names_all[[i]] = paste0(vars[[i]], colnames(xi))
                } else {
                    names_all[[i]] = paste0(vars[[i]], 1:NCOL(xi))
                }
            }
        }

        id_mat = do.call(rbind, id_all)
        values_vec = unlist(values_all)
        names_vec = unlist(names_all)

        # 2) filling the matrix: one shot, no copies

        mat = Matrix(0, nrow(data), total_cols, dimnames = list(NULL, names_vec))
        mat[id_mat] = values_vec
    }

    #
    # Step 2: the fixed-effects
    #

    if(length(x$fixef_id) == 0){
        mat_FE = NULL
    } else {
        # Same process, but easier

        rowid = 1:nrow(data)
        total_cols = sum(x$fixef_sizes)
        running_cols = c(0, x$fixef_sizes)
        n_FE = length(x$fixef_sizes)
        id_all = names_all = vector("list", n_FE)
        for(i in 1:n_FE){
            xi = x$fixef_id[[i]]
            id_all[[i]] = cbind(rowid, running_cols[i] + xi)
            names_all[[i]] = paste0(names(x$fixef_id)[i], "::", attr(xi, "fixef_names"))
        }

        id_mat = do.call(rbind, id_all)
        names_vec = unlist(names_all)

        mat_FE = Matrix(0, nrow(data), total_cols, dimnames = list(NULL, names_vec))
        mat_FE[id_mat] = 1
    }

    res = list(mat_RHS = mat, mat_FE = mat_FE)

    res
}

# Internal: modifies the calls so that each variable/interaction is evaluated with mult_sparse
mult_wrap = function(x){
    # x: character string of a variable to be evaluated
    # ex: "x1" => mult_sparse(x1)
    #     "x1:factor(x2):x3" => mult_sparse(x3, factor(x2), x1)
    #
    # We also add the argument sparse to i()
    #     "x1:i(species, TRUE)" => mult_sparse(x1, i(species, TRUE, sparse = TRUE))

    x_call = str2lang(x)

    res = (~ mult_sparse())[[2]]

    if(length(x_call) == 1 || x_call[[1]] != ":"){
        res[[2]] = x_call

    } else {
        res[[2]] = x_call[[3]]
        tmp = x_call[[2]]

        while(length(tmp) == 3 && tmp[[1]] == ":"){
            res[[length(res) + 1]] = tmp[[3]]
            tmp = tmp[[2]]
        }

        res[[length(res) + 1]] = tmp
    }

    # We also add sparse to i() if found
    for(i in 2:length(res)){
        ri = res[[i]]
        if(length(ri) > 1 && ri[[1]] == "i"){
            ri[["sparse"]] = TRUE
            res[[i]] = ri
        }
    }

    if(length(res) > 2){
        # we restore the original order
        res[-1] = rev(res[-1])
    }

    return(res)
}

# Internal function to evaluate the variables (and interactions) in a sparse way
mult_sparse = function(...){
    # Only sparsifies factor variables
    # Takes care of interactions

    dots = list(...)
    n = length(dots)

    num_var = NULL
    factor_list = list()
    info_i = NULL
    is_i = is_factor = FALSE
    # You can't have interactions between i and factors, it's either

    for(i in 1:n){
        xi = dots[[i]]
        if(is.numeric(xi)){
            # We stack the product
            num_var = if(is.null(num_var)) xi else xi * num_var
        } else if(inherits(xi, "i_sparse")){
            is_i = TRUE
            info_i = xi
        } else {
            is_factor = TRUE
            factor_list[[length(factor_list) + 1]] = xi
        }
    }

    if(!is_i && !is_factor){
        return(num_var)
    }

    if(is_factor){
        factor_list$add_items = TRUE
        factor_list$items.list = TRUE

        fact_as_int = do.call(to_integer, factor_list)

        values = if(is.null(num_var)) rep(1, length(fact_as_int$x)) else num_var

        rowid = seq_along(values)
        res = list(rowid = rowid, colid = fact_as_int$x, values = values,
                   col_names = fact_as_int$items, n_cols = length(fact_as_int$items))
    } else {

        values = info_i$values
        if(!is.null(num_var)){
            num_var = num_var[info_i$rowid]
            values = values * num_var
        }

        res = list(rowid = info_i$rowid, colid = info_i$colid, values = values,
                   col_names = info_i$col_names, n_cols = length(info_i$col_names))
    }

    class(res) = "sparse_var"

    res
}
