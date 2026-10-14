# Exhaustive fallback over recursive, axis-aligned cuts. The usual greedy
# design is attempted first, preserving successful designs from V2.0. States
# are canonical partitions, so splitting independent groups in another order
# does not consume the search budget twice. An iterative DFS avoids R's
# recursion-depth limit when many pools are requested.
.pool_balanced_search <- function(X, k, min_size, search_budget) {
  visited <- new.env(hash = TRUE, parent = emptyenv())
  states <- 0L
  make_frame <- function(members) {
    key <- paste(sort(vapply(members, function(ix)
      paste(sort(ix), collapse = ","), "")), collapse = "|")
    if (exists(key, envir = visited, inherits = FALSE)) return(NULL)
    if (states >= search_budget)
      stop("Balanced pooling search exhausted search_budget=", search_budget,
           "; feasibility is unresolved. Increase search_budget or revise k/min_size.",
           call. = FALSE)
    assign(key, TRUE, envir = visited)
    states <<- states + 1L
    if (length(members) == k) return(list(members = members, complete = TRUE))
    cuts <- list()
    for (g in seq_along(members)) {
      ix <- members[[g]]
      if (length(ix) < 2 * min_size) next
      other_capacity <- sum(floor(lengths(members[-g]) / min_size))
      for (j in seq_len(ncol(X))) {
        ord <- ix[order(X[ix, j], rownames(X)[ix])]
        mid <- which(diff(X[ord, j]) > 0)
        mid <- mid[mid >= min_size & mid <= length(ix) - min_size]
        mid <- mid[other_capacity + floor(mid / min_size) +
                     floor((length(ix) - mid) / min_size) >= k]
        if (!length(mid)) next
        mid <- mid[order(abs(mid - length(ix) / 2), mid)]
        cuts[[length(cuts) + 1L]] <- cbind(group = g, trait = j, at = mid)
      }
    }
    list(members = members, complete = FALSE,
         cuts = if (length(cuts)) do.call(rbind, cuts) else matrix(integer(), 0L, 3L),
         next_cut = 1L)
  }
  stack <- list(make_frame(list(seq_len(nrow(X)))))
  while (length(stack)) {
    level <- length(stack)
    frame <- stack[[level]]
    if (frame$complete)
      return(list(members = frame$members, states_visited = states))
    if (frame$next_cut > nrow(frame$cuts)) {
      stack[[level]] <- NULL
      next
    }
    cut <- frame$cuts[frame$next_cut, ]
    stack[[level]]$next_cut <- frame$next_cut + 1L
    members <- frame$members
    ix <- members[[cut[1L]]]
    ord <- ix[order(X[ix, cut[2L]], rownames(X)[ix])]
    members[[cut[1L]]] <- ord[seq_len(cut[3L])]
    members[[length(members) + 1L]] <- ord[seq.int(cut[3L] + 1L, length(ord))]
    child <- make_frame(members)
    if (!is.null(child)) stack[[level + 1L]] <- child
  }
  stop("No feasible recursive balanced partition exists without splitting ties or violating min_size",
       call. = FALSE)
}

design_pools <- function(phenotypes, k = 4, method = c("quantile", "balanced"),
                         min_size = 2L, search_budget = 10000L) {
  method <- match.arg(method)
  Y <- .pg_matrix(phenotypes, "phenotypes")
  .pg_scalar(min_size, "min_size", 1, nrow(Y), integer = TRUE)
  if (!is.numeric(k) || !length(k) || any(!is.finite(k)) || any(k < 1 | k != floor(k)))
    stop("k must contain positive integers", call. = FALSE)
  n <- nrow(Y)
  d <- ncol(Y)
  # Canonical row order makes designs invariant to input ordering, including ties.
  input_ids <- rownames(Y)
  Y <- Y[order(rownames(Y)), , drop = FALSE]
  breaks <- list()
  bin_index <- NULL
  search <- NULL
  if (method == "quantile") {
    if (!length(k) %in% c(1L, d)) stop("For quantiles, k must be scalar or one value per trait", call. = FALSE)
    k <- rep(k, length.out = d)
    if (any(k > n)) stop("Quantiles per trait cannot exceed the number of individuals", call. = FALSE)
    bin_index <- matrix(1L, n, d, dimnames = dimnames(Y))
    for (j in seq_len(d)) {
      v <- Y[, j]
      # Compute order positions directly: quantile(type=1) at rounded fractions
      # can otherwise move an exact integer boundary by one observation.
      thresholds <- if (k[j] == 1L) numeric() else
        unique(sort(v)[ceiling(seq_len(k[j] - 1L) * n / k[j])])
      thresholds <- thresholds[thresholds < max(v)]
      breaks[[colnames(Y)[j]]] <- thresholds
      bin_index[, j] <- findInterval(v, thresholds, left.open = TRUE) + 1L
    }
    labels <- apply(bin_index, 1L, function(x) paste0("P_", paste(x, collapse = "_")))
  } else {
    .pg_scalar(search_budget, "search_budget", 1, .Machine$integer.max, integer = TRUE)
    if (length(k) != 1L) stop("Balanced pooling requires one total number of pools", call. = FALSE)
    if (k < 2 || k > floor(n / min_size)) stop("Requested number/size of balanced pools is infeasible", call. = FALSE)
    sds <- apply(Y, 2L, stats::sd)
    sds[!is.finite(sds) | sds == 0] <- 1
    X <- sweep(sweep(Y, 2L, colMeans(Y), "-"), 2L, sds, "/")
    members <- list(seq_len(n))
    search <- list(method = "greedy", states_visited = 0L, budget = search_budget)
    while (length(members) < k) {
      best <- NULL
      for (g in seq_along(members)) {
        ix <- members[[g]]
        if (length(ix) < 2 * min_size) next
        for (j in seq_len(d)) {
          ord <- ix[order(X[ix, j], rownames(X)[ix])]
          split_at <- which(diff(X[ord, j]) > 0)
          split_at <- split_at[split_at >= min_size & split_at <= length(ix) - min_size]
          # A greedy half split (e.g. 12 -> 6+6 with min_size=4) can
          # destroy a feasible requested design of three pools of four.
          other_capacity <- sum(floor(lengths(members[-g]) / min_size))
          split_at <- split_at[other_capacity + floor(split_at / min_size) +
                                  floor((length(ix) - split_at) / min_size) >= k]
          if (!length(split_at)) next
          # Closest feasible split to half keeps counts balanced and never splits equal values.
          mid <- split_at[which.min(abs(split_at - length(ix) / 2))]
          left <- ord[seq_len(mid)]
          right <- ord[seq.int(mid + 1L, length(ord))]
          gain <- length(left) * length(right) / length(ix) *
            sum((colMeans(X[left, , drop = FALSE]) - colMeans(X[right, , drop = FALSE]))^2)
          if (is.null(best) || gain > best$gain)
            best <- list(group = g, left = left, right = right, gain = gain)
        }
      }
      if (is.null(best)) {
        recovered <- .pool_balanced_search(X, k, min_size, search_budget)
        members <- recovered$members
        search$method <- "backtracking"
        search$states_visited <- recovered$states_visited
        break
      }
      members[[best$group]] <- best$left
      members[[length(members) + 1L]] <- best$right
    }
    labels <- character(n)
    for (g in seq_along(members)) labels[members[[g]]] <- sprintf("P_%0*d", nchar(k), g)
  }
  names(labels) <- rownames(Y)
  sizes <- table(labels)
  if (length(sizes) < 2L) stop("Phenotypes and ties yield fewer than two occupied pools", call. = FALSE)
  if (any(sizes < min_size))
    stop("Quantile grid contains pools smaller than min_size; reduce k or use balanced pooling", call. = FALSE)
  pool_ids <- sort(unique(labels))
  bounds <- do.call(rbind, lapply(pool_ids, function(pool) {
    rows <- which(labels == pool)
    data.frame(pool_id = pool, trait = colnames(Y), n = length(rows),
               lower = apply(Y[rows, , drop = FALSE], 2L, min),
               upper = apply(Y[rows, , drop = FALSE], 2L, max),
               mean = colMeans(Y[rows, , drop = FALSE]), row.names = NULL)
  }))
  clusters <- factor(labels[input_ids], levels = pool_ids)
  names(clusters) <- input_ids
  if (!is.null(bin_index)) bin_index <- bin_index[input_ids, , drop = FALSE]
  structure(list(clusters = clusters, bounds = bounds, pool_sizes = sizes[pool_ids],
                 method = method, k = k, min_size = min_size, breaks = breaks,
                 bin_index = bin_index, traits = colnames(Y),
                 search = search,
                 ties = "identical values are never split along a cut"), class = "pool_design")
}

quantile_clustering <- function(data, n_quantiles = 4, min_size = 1L) {
  design_pools(data, k = n_quantiles, method = "quantile", min_size = min_size)$clusters
}

print.pool_design <- function(x, ...) {
  cat("PooGwaM", x$method, "design:", length(x$pool_sizes), "occupied pools; sizes",
      min(x$pool_sizes), "to", max(x$pool_sizes), "\n")
  invisible(x)
}

# The public entry point returns a fit, never an uncalibrated p-value vector.
poogwam <- function(data, method = c("ridge", "multitrait"), ...) {
  method <- match.arg(method)
  if (!inherits(data, "pool_gwas_data"))
    stop("V2 expects as_pool_gwas_data(...) or pool_from_individuals(...); see migration guide", call. = FALSE)
  if (method == "ridge") fit_pool_ridge(data, ...) else fit_pool_multitrait(data, ...)
}
