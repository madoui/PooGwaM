# Formula covariates are encoded at their declared observation level. In
# particular, a mean of an interaction is not an interaction of the means.
.pg_covariate_environment <- function() {
  env <- new.env(parent = baseenv())
  env$poly <- stats::poly
  env
}

.pg_covariate_expression <- function(x, generated = FALSE, inside_transform = FALSE) {
  if (is.null(x) || is.atomic(x) || is.name(x)) return(invisible(TRUE))
  if (is.list(x) && !is.call(x)) {
    for (value in x) .pg_covariate_expression(value, generated, inside_transform)
    return(invisible(TRUE))
  }
  allowed <- c("~", "+", "-", "*", "/", "^", ":", "(", "I",
               "log", "log1p", "sqrt", "exp", "abs", "scale", "poly")
  if (generated) allowed <- c(allowed, "list")
  if (!is.call(x) || !is.name(x[[1L]]) ||
      !as.character(x[[1L]]) %in% allowed)
    stop("Unsupported covariate transformation: use arithmetic, I, log, log1p, sqrt, exp, abs, scale or poly",
         call. = FALSE)
  name <- as.character(x[[1L]])
  if (inside_transform && name %in% c("scale", "poly"))
    stop("scale and poly must be standalone formula terms, not nested inside another transformation",
         call. = FALSE)
  inside_transform <- inside_transform ||
    name %in% c("I", "log", "log1p", "sqrt", "exp", "abs", "scale", "poly")
  for (value in as.list(x)[-1L]) .pg_covariate_expression(value, generated, inside_transform)
  invisible(TRUE)
}

.pg_covariate_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 2L)
    stop("covariates must be a one-sided formula without a response", call. = FALSE)
  .pg_covariate_expression(formula)
  if ("." %in% all.names(formula))
    stop("List covariate variables explicitly; '.' is not supported", call. = FALSE)
  environment(formula) <- .pg_covariate_environment()
  tt <- stats::terms(formula)
  if (attr(tt, "intercept") != 1L)
    stop("Covariate formulas must retain the model's intercept; ~0 and -1 are not supported",
         call. = FALSE)
  formula
}

.pg_covariate_metadata <- function(metadata, ids, variables) {
  .pg_ids(ids, length(ids), "covariate target rows")
  if (is.null(metadata)) {
    if (length(variables)) stop("A covariate formula requires metadata", call. = FALSE)
    metadata <- data.frame(row.names = ids)
  }
  if (!is.data.frame(metadata)) stop("metadata must be a data frame with row IDs", call. = FALSE)
  .pg_ids(rownames(metadata), nrow(metadata), "metadata rows")
  if (!setequal(rownames(metadata), ids)) stop("metadata IDs do not match", call. = FALSE)
  if (anyNA(names(metadata)) || any(!nzchar(names(metadata))) || anyDuplicated(names(metadata)))
    stop("metadata columns must have unique, nonempty names", call. = FALSE)
  if (!all(variables %in% names(metadata)))
    stop("All formula variables must be columns of metadata; missing: ",
         paste(setdiff(variables, names(metadata)), collapse = ", "), call. = FALSE)
  metadata <- metadata[ids, , drop = FALSE]
  for (name in variables) {
    x <- metadata[[name]]
    if (!is.null(dim(x)) || !(is.numeric(x) || is.logical(x) || is.factor(x) || is.character(x)))
      stop("metadata variable '", name, "' must be numeric, logical, character or factor", call. = FALSE)
    if (anyNA(x) || (is.numeric(x) && any(!is.finite(x))))
      stop("metadata variable '", name, "' contains missing or nonfinite values", call. = FALSE)
    if (is.character(x)) metadata[[name]] <- factor(x)
    if (is.logical(x)) metadata[[name]] <- factor(x, levels = c(FALSE, TRUE))
  }
  metadata
}

.pg_covariate_matrix <- function(frame, terms, contrasts = NULL, columns = NULL) {
  matrix <- stats::model.matrix(terms, frame, contrasts.arg = contrasts)
  matrix <- matrix[, colnames(matrix) != "(Intercept)", drop = FALSE]
  actual_columns <- if (ncol(matrix)) colnames(matrix) else character()
  if (!is.null(columns) && !identical(actual_columns, columns))
    stop("Encoded covariate columns do not match the fitted specification", call. = FALSE)
  if (!ncol(matrix)) return(NULL)
  # Strip model.matrix attributes: the specification records this information.
  attributes(matrix) <- list(dim = dim(matrix), dimnames = dimnames(matrix))
  .pg_matrix(matrix, "encoded covariates")
}

.pg_covariates <- function(covariates, metadata, ids, level) {
  level <- match.arg(level, c("pool", "individual"))
  if (!inherits(covariates, "formula")) {
    if (!is.null(metadata)) stop("metadata requires a covariate formula", call. = FALSE)
    matrix <- if (is.null(covariates)) NULL else
      .pg_align(.pg_matrix(covariates, "covariates"), ids, what = "covariates")
    if (!is.null(matrix) && "(Intercept)" %in% colnames(matrix))
      stop("Do not supply an intercept column; the model includes one", call. = FALSE)
    return(list(matrix = matrix, spec = NULL))
  }
  formula <- .pg_covariate_formula(covariates)
  metadata <- .pg_covariate_metadata(metadata, ids, all.vars(formula))
  frame <- stats::model.frame(formula, data = metadata, na.action = stats::na.fail,
                             drop.unused.levels = FALSE)
  terms <- attr(frame, "terms")
  # makepredictcall records the training centre/scale and polynomial basis in
  # predvars; using these terms later avoids fitting transformations to new data.
  factors <- vapply(frame, is.factor, logical(1))
  xlevels <- lapply(frame[factors], levels)
  single <- names(xlevels)[lengths(xlevels) < 2L]
  if (length(single))
    stop("Covariate factor(s) require 2 or more levels: ", paste(single, collapse = ", "), call. = FALSE)
  contrasts <- lapply(frame[factors], stats::contrasts, contrasts = TRUE)
  if (!length(contrasts)) contrasts <- NULL
  matrix <- .pg_covariate_matrix(frame, terms, contrasts)
  spec <- list(version = 1L, input_level = level, level = level,
               formula = formula, terms = terms, predvars = attr(terms, "predvars"),
               xlevels = xlevels, contrasts = contrasts,
               columns = if (is.null(matrix)) character() else colnames(matrix))
  list(matrix = matrix, spec = spec)
}

.pg_check_covariate_spec <- function(spec) {
  if (!is.list(spec) || !identical(spec$version, 1L) ||
      !identical(spec$input_level, spec$level) ||
      !is.character(spec$input_level) || length(spec$input_level) != 1L ||
      !spec$input_level %in% c("pool", "individual") ||
      !inherits(spec$terms, "terms") || !is.list(spec$xlevels) ||
      !is.call(spec$predvars) || !identical(spec$predvars[[1L]], as.name("list")) ||
      !(is.null(spec$contrasts) || is.list(spec$contrasts)) ||
      !is.character(spec$columns) || anyNA(spec$columns) ||
      any(!nzchar(spec$columns)) || anyDuplicated(spec$columns))
    stop("Invalid covariate_spec", call. = FALSE)
  formula <- .pg_covariate_formula(spec$formula)
  base_terms <- stats::terms(formula)
  semantic_attributes <- c("variables", "factors", "term.labels", "order", "intercept", "response")
  if (any(!vapply(semantic_attributes, function(name)
        identical(attr(base_terms, name), attr(spec$terms, name)), logical(1))) ||
      !identical(spec$predvars, attr(spec$terms, "predvars")) ||
      attr(spec$terms, "intercept") != 1L)
    stop("Inconsistent formula terms in covariate_spec", call. = FALSE)
  .pg_covariate_expression(spec$predvars, generated = TRUE)
  invisible(TRUE)
}

.pg_covariate_groups <- function(clusters, ids, pool_sizes = NULL) {
  if (inherits(clusters, "pool_design")) clusters <- clusters$clusters
  .pg_ids(names(clusters), length(clusters), "clusters individuals")
  if (!length(clusters)) stop("Individual covariates require named clusters", call. = FALSE)
  groups <- as.character(clusters)
  names(groups) <- names(clusters)
  if (anyNA(groups) || any(!nzchar(groups)) || !setequal(unique(groups), ids))
    stop("clusters must cover exactly the prediction pool IDs", call. = FALSE)
  sizes <- setNames(as.numeric(table(factor(groups, levels = ids))), ids)
  if (!is.null(pool_sizes)) {
    if (!is.numeric(pool_sizes)) stop("pool_sizes must be numeric", call. = FALSE)
    .pg_ids(names(pool_sizes), length(pool_sizes), "pool_sizes")
    if (!setequal(names(pool_sizes), ids) || any(!is.finite(pool_sizes)) ||
        any(pool_sizes[ids] != sizes))
      stop("clusters membership sizes do not match pool_sizes", call. = FALSE)
  }
  list(groups = groups, sizes = sizes)
}

.pg_predict_covariates <- function(spec, metadata, ids, clusters = NULL, pool_sizes = NULL) {
  .pg_check_covariate_spec(spec)
  .pg_ids(ids, length(ids), "prediction pools")
  variables <- all.vars(spec$formula)
  if (spec$input_level == "individual") {
    grouped <- .pg_covariate_groups(clusters, ids, pool_sizes)
    input_ids <- names(grouped$groups)
  } else {
    if (!is.null(clusters)) stop("Pool-level covariate formulas do not accept individual clusters", call. = FALSE)
    input_ids <- ids
  }
  metadata <- .pg_covariate_metadata(metadata, input_ids, variables)
  terms <- spec$terms
  environment(terms) <- .pg_covariate_environment()
  frame <- stats::model.frame(terms, data = metadata, na.action = stats::na.fail,
                             xlev = spec$xlevels, drop.unused.levels = FALSE)
  matrix <- .pg_covariate_matrix(frame, terms, spec$contrasts, spec$columns)
  if (is.null(matrix) || spec$input_level == "pool") return(matrix)
  rowsum(matrix, grouped$groups, reorder = FALSE)[ids, , drop = FALSE] /
    as.numeric(grouped$sizes)
}

.pg_validate_covariate_spec <- function(data) {
  if (is.null(data$covariate_spec)) {
    if (!is.null(data$covariate_metadata) || !is.null(data$covariate_clusters))
      stop("Stored covariate metadata requires covariate_spec", call. = FALSE)
    return(invisible(TRUE))
  }
  expected <- .pg_predict_covariates(data$covariate_spec, data$covariate_metadata,
    rownames(data$frequencies), data$covariate_clusters, data$pool_sizes)
  actual <- data$covariates
  if (!is.null(actual)) actual <- actual[rownames(data$frequencies), , drop = FALSE]
  if (!isTRUE(all.equal(actual, expected, tolerance = 1e-12, check.attributes = TRUE)))
    stop("Stored covariates are inconsistent with covariate_spec and metadata", call. = FALSE)
  invisible(TRUE)
}
