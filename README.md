# PooGwaM 2.2.0

An R package for analyzing allele frequencies and multiple traits measured in pools. It provides ridge regression, a multitrait sparse-group model, formula-based fixed covariates, explicit calibration, experimental knockoffs, and Manhattan plots. The package includes code, documentation, and lightweight tests; research simulations and benchmarks are kept separately.

## Installation

```r
# From a local archive of version 2.2.0:
install.packages("PooGwaM_2.2.0.tar.gz", repos = NULL, type = "source")

# For plotting only:
install.packages(c("qqman", "CMplot"))
```

R >= 4.1.0 is required. The core uses the standard `stats`, `tools`, and `graphics` packages; `qqman` and `CMplot` are optional. The package's GitHub repository is [madoui/PooGwaM](https://github.com/madoui/PooGwaM). V1 and V2 have the same package name: use separate libraries to retain multiple versions.

<a id="comparer-ridge-identity-et-sparse"></a>

## Comparing ridge, identity, and sparse

`ridge` and `multitrait` are the two methods in the `poogwam()` API. **`identity` is a penalty option**, available for both solvers; the benchmark label "sparse" denotes a configuration of `fit_pool_multitrait()`. The three configurations compared in the benchmarks are:

| Label | Function and settings | Behavior |
|---|---|---|
| **PooGwaM ridge** | `fit_pool_ridge()`, `lambda=1`, `penalty="within_pool"` | Shrinks coefficients with a SNP-specific quadratic penalty, without imposing sparse selection. |
| **PooGwaM identity** | `fit_pool_ridge()`, `lambda=1`, `penalty="identity"` | Uses the same ridge solver, with the same quadratic penalty weight for every SNP. |
| **PooGwaM sparse** | `fit_pool_multitrait()`, `lambda_ridge=1`, `lambda_group=0.01`, `lambda_l1=0.002`, `penalty="within_pool"` | Adds penalties that can set coefficients exactly to zero. This benchmark label uses one trait at a time. |

Given allele frequencies $F_{kj}$, pool weights $w_k=n_k/\sum_k n_k$, and dosage $X=2F$, ridge regression minimizes the weighted residual sum of squares plus $\lambda\sum_{j,t}D_j B_{jt}^2$. The intercept and any covariates are unpenalized. The raw weight is:

$$
D_j^{\mathrm{raw}} =
\begin{cases}
\sum_k w_k\,2F_{kj}(1-F_{kj}) & \text{for the within-pool penalty}\\
1 & \text{for the identity penalty}
\end{cases}
$$

The code then applies `D_j = max(D_j_raw, penalty_floor)`, with `penalty_floor=1e-8` by default. `within_pool` generally gives rare variants a smaller weight. This can contribute to differences in coefficient and null-score scales across allele frequencies. `identity` removes this particular weighting; it does not by itself guarantee balanced power. Both options use the same deterministic solver on a pools × pools matrix. See the [penalty preparation and ridge solver](R/fit.R).

The sparse model minimizes half the quadratic criterion, with `lambda_ridge`, plus `lambda_group` times the sum of the norms of the effect vectors for each variant and `lambda_l1` times the sum of the absolute coefficient values. With multiple traits, the group penalty can remove a variant across all traits; the L1 penalty also allows trait-specific zeros. **With a single trait, the two sparse penalties add together: the settings above give an L1 penalty of 0.012**, without sharing information across traits.

The sparse solver uses accelerated proximal gradient descent, a line search, and a KKT convergence check. With `lambda_group=0` and `lambda_l1=0`, it calls the exact ridge solver directly with the same remaining settings. Traits are standardized by default for sparse (`standardize_traits=TRUE`), versus `FALSE` for ridge; returned effects are always converted back to their original units. See the [sparse-group solver](R/multitrait.R).

## Examples: ridge, identity, and sparse

This example creates the required data in R without an external simulation file. The first three fits use the same trait; the final fit separately illustrates information sharing across traits.

```r
library(PooGwaM)
stopifnot(packageVersion("PooGwaM") == "2.2.0")

pool_ids <- paste0("pool", 1:4)
F <- matrix(c(.1, .3, .7, .9, .2, .4, .5, .8), 4, 2,
            dimnames = list(pool_ids, c("snp1", "snp2")))
Y <- cbind(height = c(-2, -1, 1, 2), mass = c(3, 1, 4, 2))
rownames(Y) <- pool_ids
sizes <- setNames(c(20, 25, 20, 30), pool_ids)
variants <- data.frame(variant_id = colnames(F),
                       chromosome = c("1", "2"), position = c(100, 200))
dat <- as_pool_gwas_data(F, Y, sizes, variants = variants)
dat_height <- as_pool_gwas_data(F, Y[, "height", drop = FALSE], sizes,
                                variants = variants)

# Ridge with a penalty based on within-pool frequencies:
ridge <- fit_pool_ridge(dat_height, lambda = 1, penalty = "within_pool",
                        standardize_traits = FALSE)

# The same ridge solver with a uniform penalty across SNPs:
ridge_identity <- fit_pool_ridge(dat_height, lambda = 1, penalty = "identity",
                                 standardize_traits = FALSE)

# Sparse with a single trait: total L1 penalty of 0.01 + 0.002.
sparse <- fit_pool_multitrait(dat_height, lambda_ridge = 1,
  lambda_group = .01, lambda_l1 = .002, penalty = "within_pool",
  standardize_traits = TRUE)

# Joint sparse-group fit for height and mass:
joint <- poogwam(dat, method = "multitrait",
  lambda_ridge = 1, lambda_group = .01, lambda_l1 = .002,
  penalty = "within_pool", standardize_traits = TRUE)

stopifnot(all(vapply(list(ridge, ridge_identity, sparse, joint),
  function(fit) isTRUE(fit$diagnostics$converged), logical(1))))
cbind(ridge = ridge$beta[, "height"],
      identity = ridge_identity$beta[, "height"],
      sparse = sparse$beta[, "height"])
sparse$diagnostics$selected_variants
head(association_table(ridge))
```

Effects are expressed per additional alternative allele. Without calibration, results are exploratory scores and p-values are missing. Data are aligned by identifiers; frequencies and any counts must be consistent. Fitting functions revalidate objects, including those modified after construction.

These hyperparameters illustrate the benchmark configurations; they are not universally optimal settings. SNPs with nonzero coefficients in the sparse model are not automatically significant discoveries. The API also accepts `poogwam(dat_height, method="ridge", penalty="identity")`; `method="identity"` and `method="sparse"` are not valid method values.

## Adjusting for sex, age, batch, site, and treatment

V2.2 accepts `covariates = ~ sex + age + batch + site + treatment` and a `metadata` data frame whose row names are individual or pool identifiers. Covariates and the intercept are **unpenalized fixed effects** in all three configurations. They are neither random effects nor an automatic correction for relatedness or genetic structure.

The following example uses 96 fictitious individuals divided into 16 pools. Sex proportions and other covariates vary across pools. There are enough pools for the seven columns of the covariate design, including the intercept, leaving residual dimensions for the genetic fit.

```r
library(PooGwaM)

individual_ids <- sprintf("person%03d", 1:96)
pool_number <- rep(1:16, each = 6)
within_pool <- rep(1:6, 16)
clusters <- setNames(sprintf("pool%02d", pool_number), individual_ids)

metadata <- data.frame(
  sex = factor(ifelse(within_pool <= 1 + (pool_number - 1) %% 4, "M", "F"),
                levels = c("F", "M")),
  age = 25 + 1.3 * (pool_number - 1) + within_pool + 2 * (pool_number %% 3),
  batch = factor(paste0("L", 1 + (pool_number - 1) %% 3),
               levels = c("L1", "L2", "L3")),
  site = factor(ifelse((pool_number - 1) %% 8 < 4, "North", "South"),
                levels = c("North", "South")),
  treatment = factor(ifelse((pool_number - 1) %% 2 == 0, "control", "treated"),
                      levels = c("control", "treated")),
  row.names = individual_ids
)

# Small deterministic dosage and trait data keep the example self-contained.
G <- vapply(1:4, function(j)
  ((1:96) * (2 * j + 1) + (1:96) %/% (j + 1)) %% 3, numeric(96))
dimnames(G) <- list(individual_ids, paste0("snp", 1:4))
Y <- cbind(
  height = 150 + .6 * metadata$age + 4 * (metadata$sex == "M") +
    1.2 * (metadata$treatment == "treated") + .8 * G[, 1] + 2 * sin((1:96) / 7),
  mass = 45 + .4 * metadata$age + 3 * (metadata$site == "South") +
    1.5 * (metadata$batch == "L3") + G[, 2] + cos((1:96) / 5)
)
rownames(Y) <- individual_ids

dat_cofactors <- pool_from_individuals(G, Y, clusters,
  covariates = ~ sex + age + batch + site + treatment,
  metadata = metadata)

ridge_cofactors <- fit_pool_ridge(dat_cofactors, lambda = 1,
                                 penalty = "within_pool")
identity_cofactors <- fit_pool_ridge(dat_cofactors, lambda = 1,
                                    penalty = "identity")
sparse_cofactors <- fit_pool_multitrait(dat_cofactors,
  lambda_ridge = 1, lambda_group = .01, lambda_l1 = .002,
  penalty = "within_pool", standardize_traits = TRUE)

ridge_cofactors$covariate_effects  # coefficients for each covariate and trait
ridge_cofactors$intercept
ridge_cofactors$diagnostics$covariate_design
stopifnot(all(vapply(list(ridge_cofactors, identity_cofactors, sparse_cofactors),
  function(fit) isTRUE(fit$diagnostics$converged), logical(1))))

# Re-encode individual metadata using the learned schema and the same pools.
# newdata contains pool frequencies, not individual dosages.
pred <- predict(ridge_cofactors, newdata = dat_cofactors$frequencies,
                metadata = metadata, clusters = clusters)
stopifnot(isTRUE(all.equal(pred, ridge_cofactors$fitted)))

# A data object also retains the metadata needed for prediction:
predict(ridge_cofactors, newdata = dat_cofactors)
```

Encoding takes place **before** individuals are aggregated into pools. With R's usual treatment contrasts, an indicator column becomes the proportion of individuals in the corresponding category in each pool. The reference levels in this example are `F`, `L1`, `North`, and `control`; setting factor levels explicitly before fitting lets you choose these references. The learned schema and contrasts are retained for prediction, even if a known category is absent from the new observations. An unknown level is rejected.

A formula such as `~ age + I(age^2) + sex:treatment` is also encoded at the individual level before averaging. The squared-age column therefore contains the **mean of the squares**, rather than the square of mean age; the interaction follows the same rule. For data that are already pooled, `as_pool_gwas_data(..., covariates = ~ ..., metadata = metadata_pools)` instead applies the formula directly to pool rows. Means alone cannot recover transformations of individual values.

Direct `scale(age)` and `poly(age, 2)` terms retain their learned parameters for prediction. The usual transformations `I`, `log`, `log1p`, `sqrt`, `exp`, and `abs` are also supported. Custom functions, variables outside `metadata`, and nested expressions such as `I(scale(age))` or `poly(scale(age), 2)` are rejected; list variables explicitly instead of using `~ .`, and retain the formula intercept. Logical variables retain the levels `FALSE` and `TRUE` for prediction. A factor declared with only one level is explicitly rejected.

Existing numeric matrices remain accepted by `covariates`: they are treated as already encoded, without an intercept column. At prediction time, their column names must match those used for training. For a formula learned on individuals, supply individual metadata and `clusters`; for a formula learned on pools, supply the metadata for those pools. All joins use identifiers. Missing values, missing individuals, and unknown categories raise explicit errors; no individual is silently dropped.

Objects constructed with a formula retain the metadata and pool memberships required for validation and prediction, even with `retain_individuals = FALSE`. This option controls retention of individual genotypes and phenotypes, rather than the covariate schema.

Rank must be checked **after pooling**: a batch confounded with site, an absent category, or identical proportions across all pools can make columns redundant. The program then rejects the fit and identifies the affected columns. `diagnostics$covariate_design$residual_df` describes K minus the rank of `cbind(1, C)`, where K is the number of pools and C contains the encoded covariates; it is not a classical degrees-of-freedom count for ridge t-tests. More individuals cannot compensate for too few pools.

**This initial covariate support covers effect estimation and prediction.** The calibration and knockoff functions have not gained an inference model adjusted for these covariates. Scores from `association_table(ridge_cofactors)` remain exploratory, with missing p-values in the absence of compatible calibration; nonzero sparse coefficients are not significant discoveries by themselves.

## Model history

| Date / version | Code changes |
|---|---|
| **24 November 2022** — [ba8e7db](https://github.com/madoui/PooGwaM/commit/ba8e7db) | The former `R/selection.R` used LASSO and stability selection (`glmnet`, `stabs`), separately for two traits. This procedure differs from the current sparse-group model. |
| **7 September 2023** — [06a68ae](https://github.com/madoui/PooGwaM/commit/06a68ae) | Introduced `adaptRidge()` with stochastic gradient updates and a frequency-dependent penalty; removed the former selection module. |
| **14 September 2023** — [43f69a9](https://github.com/madoui/PooGwaM/commit/43f69a9) | Introduced the `poogwam()` wrapper, which called `adaptRidge()` and then transformed coefficients into p-values. |
| **16 October 2023** — [a92744c](https://github.com/madoui/PooGwaM/commit/a92744c) | Added `singValRidge()`, based on an SVD and a common penalty `mean(alpha)`. This average penalty is not the current explicit `identity` option. |
| **V2, local delivery documented on 25 September 2026; published on 13 October 2026** — [89dbcef](https://github.com/madoui/PooGwaM/commit/89dbcef) | Introduced exact ridge regression on pools, the `identity` option, and the new sparse-group solver. Uses dosage `2F`, unpenalized intercept/covariates, and separate fitting and p-value calibration. |
| **V2.1, 14 October 2026** — [0c5c23a](https://github.com/madoui/PooGwaM/commit/0c5c23a) | Retained the ridge and sparse-group formulations. Strengthened input validation and added studentized maxT, experimental knockoffs, and Manhattan plots. |
| **V2.2.0, 15 October 2026** | Added formula-based fixed covariates and metadata matched by identifiers, with encoding before pooling and reuse of the learned schema for predictions. This extension covers fitting and prediction; the previous inference scopes are unchanged. |

In archived V1 (0.1.51), `adaptRidge()` used `X=ZF` and weights `2*colMeans(X)*(1-colMeans(X))`, calculated from mean frequencies. V2 uses `2F` and the weighted mean of **per-pool** heterozygosities for `within_pool`. V1/V2 coefficients are therefore not directly interchangeable. V2 also removed the automatic transformation of coefficients into p-values: fitting provides effects and scores, after which an explicit inference procedure is required.

## Raw or studentized maxT

A self-contained example using small illustrative data:

```r
set.seed(42)
ids <- paste0("i", 1:80)
snps <- paste0("v", 1:40)
f <- setNames(rep(c(.04, .3), each = 20), snps)
G <- matrix(rbinom(80 * 40, 2, rep(f, each = 80)), 80,
            dimnames = list(ids, snps))
Y <- cbind(height = G[, 1] + rnorm(80), mass = G[, 21] + rnorm(80))
rownames(Y) <- ids
builder <- function(y) design_pools(y[, "height", drop = FALSE],
                                   k = 4, min_size = 2)
cal <- calibrate_pool_permutation(G, Y, builder,
  B = 99, B_pilot = 99, maxT = "studentized", seed = 123,
  fit_args = list(lambda = 1))
tab <- association_table(cal$fit, cal)
head(tab[tab$test == "omnibus", ])
cal$studentization$floored_variants
cal$rng$kind
```

Studentization estimates a SNP-specific omnibus mean and scale from an independent pilot, then freezes this transformation for the main permutations. It aims to reduce differences in null scales across SNPs. The pilot size above is illustrative, rather than a universally validated setting.

**`maxT="raw"` remains the default in V2.2.** Making studentized maxT the default awaits the planned confirmatory study of error rates and power. Both modes retain the **`global_null_only`** scope: studentization establishes neither strong FWER control under partial alternatives nor variant-level FDR control. Marginal p-values and raw maxima remain available. B permutations give a minimum resolution of `1/(B+1)`.

The pilot and main series use sampling with replacement, without removing duplicates. A degenerate pilot produces a diagnostic, never proof that the null variance is zero. The RNG state is restored after the call; algorithms with hidden state that cannot be restored (Box-Muller or user-supplied generators) are rejected. To replay in the same environment, reuse the inputs, settings, seed, and `rng_kind = cal$rng$kind`; the realized plan is also recorded.

`calibrate_pool_counts()` remains an experimental count bootstrap under homogeneous frequencies, with raw maxT. It does not automatically account for LD, relatedness, ancestry, overdispersion, or sequencing errors.

## Experimental pooled knockoffs

```r
# f is known here because it was used to generate G above.
# Do not replace it with the frequency estimated from G without justification.
G_tilde <- generate_independent_knockoffs(G, frequencies = f, seed = 456)
selection <- select_pool_knockoffs(G, G_tilde, Y, builder,
  q = .1, method = "ridge", fit_args = list(lambda = 1), seed = 123)
selection$selected
selection$W
selection$inference_status

# The multitrait solver is also available:
# selection <- select_pool_knockoffs(G, G_tilde, Y, builder,
#   method = "multitrait", fit_args = list(lambda_group = .01))
```

Originals and copies are pooled, then fitted together. The W contrast compares the norms of standardized coefficients; the knockoff+ threshold uses the `+1` correction. A discovery concerns a variant associated with the trait vector, rather than each variant–trait pair. No p-value is assigned to W.

**Status: `experimental_knockoff`.** The supplied generator is a reference implementation under a known independent binomial model. It does not generate knockoffs suitable for arbitrary LD. An external pair can be supplied, with provenance and allele metadata, but matrix validation does not certify its exchangeability. The theoretical guarantee requires the Model-X assumptions and symmetry of the entire procedure. Estimated HMMs, covariates, relatedness, and read-based frequencies require separate extensions.

## Manhattan plots with qqman or CMplot

```r
coords <- data.frame(variant_id = snps,
  chromosome = rep(c("1", "2"), each = 20),
  position = rep(seq(1000, 20000, by = 1000), 2))

# Omnibus adjusted maxT p-values, with an explicit threshold:
plot_manhattan(cal, variants = coords, value = "p_maxT",
                backend = "qqman", threshold = .05)

# Marginal p-values for one trait:
plot_manhattan(cal, variants = coords, value = "p_raw",
                test = "trait", trait = "height", backend = "CMplot")

# Signed W: no transformation into a p-value.
plot_manhattan(selection, variants = coords, value = "W")

# Prepare the data without any plotting dependency:
prepared <- manhattan_data(cal, variants = coords)
head(prepared$data)

# Export to PDF; png() can also be used:
pdf("manhattan-maxT.pdf", width = 10, height = 5)
plot_manhattan(cal, variants = coords, threshold = .05)
dev.off()
```

Coordinates are joined by `variant_id`. Column names can be customized through `chromosome` and `position`, and chromosome order through `chromosome_order`. For a single chromosome, use `qqman`: CMplot's single-chromosome axis has limitations with small ranges and labels. The plot adds no conventional human GWAS threshold automatically. Zero p-values are rejected, without silent flooring. Missing values are rejected by default; `na_action="omit"` explicitly excludes them with a warning and a list of identifiers.

## Fixes and limitations

V2.1 fixes integer overflow in AUROC, greedy pooling dead ends when another recursive partition is feasible, and validation of objects modified before fitting. Pool search is bounded by `search_budget`; exhausting that budget does not prove infeasibility.

The pools × variants matrix remains in memory. Exploratory multitrait scores, permutation calibration, and knockoffs retain their respective scopes. Alleles must retain their orientation; no method automatically corrects genetic confounding or a misspecified read model.

See [NEWS](NEWS.md), the [methodological note](inst/doc/inference-v2.1.md), and the R help pages for details. The test suite compares solvers with independent references and covers invariances, edge cases, p-values, knockoff swaps, and both plotting backends. Software checks do not replace biological or confirmatory statistical validation.