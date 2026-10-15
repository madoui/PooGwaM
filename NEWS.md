# PooGwaM 2.2.0

## Fixed cofactors from formulas and identified metadata

- Accept a one-sided formula in `covariates`, with a separate `metadata` data
  frame identified by pool or individual row names. Existing numeric covariate
  matrices remain supported without a formula.
- Encode categorical variables, interactions and supported transformations
  before aggregating individual covariates into pools. Pool-level metadata are
  encoded directly. These are unpenalized fixed effects; genetic ridge and
  sparse-group penalties retain their existing definitions.
- Store the training formula, encoding, factor levels and contrasts, and reuse
  them for prediction from `metadata` and, when applicable, `clusters`. New
  formula-prediction arguments are named; previously encoded numeric prediction
  matrices remain usable.
- Preserve learned parameters for direct `scale()` and `poly()` terms. Support
  arithmetic, `I`, `log`, `log1p`, `sqrt`, `exp` and `abs`; reject unsupported
  functions and nested stateful transforms such as `I(scale(age))` or
  `poly(scale(age), 2)`. Logical variables retain both `FALSE`/`TRUE` levels;
  factors declared with fewer than two levels are rejected explicitly.
- Report covariate design columns, rank, number of pools and residual pool
  dimension in `diagnostics$covariate_design`. Reject missing metadata, unknown
  factor levels, unidentified rows and rank-deficient pooled covariate designs
  explicitly instead of silently dropping terms or individuals.
- This release adds cofactor adjustment for effect estimation and prediction.
  Permutation calibration, count bootstrap and pooled knockoff selection do not
  gain covariate-adjusted inference. Fitted effects and sparse selections remain
  exploratory unless a separately justified inferential procedure is used.

# PooGwaM 2.1.0

## Calibration

- Add `maxT="studentized"` to individual permutation calibration. An independent
  pilot estimates SNP-specific omnibus means and scales before the main series.
  Record pilot moments, positive scale floors, degenerate-pilot diagnostics,
  transformed maxima, and `maxT_method` in association tables.
- Keep `maxT="raw"` as the default. Promotion of the studentized default requires
  the confirmatory error-rate and rare/common power study in the development plan;
  exploratory validation does not satisfy that requirement.
- Preserve raw null summaries and marginal p-values. Both maxT modes retain
  `global_null_only`; count bootstrap remains experimental and uses raw maxT.
- Record RNG algorithm, boundary states, R/package versions, realized pooling
  design, and input/settings checksums. `rng_kind` allows explicit replay.
- Reject RNG configurations whose hidden state cannot be restored through R's
  public API (Box-Muller and user-supplied generators), before changing RNG state.

## Experimental knockoffs

- Add `select_pool_knockoffs()` for original/knockoff individual genotype pairs,
  common phenotype-defined pooling, joint ridge or sparse-group fitting, and
  knockoff+ selection. Canonical unordered pairs preserve numerical sign flips.
- Add `generate_independent_knockoffs()` for a known product-binomial reference
  model with externally supplied allele frequencies. It does not estimate LD or
  certify that supplied frequencies describe the original data.
- Return W, threshold, selected variants, model assumptions and provenance,
  separately from p-value calibration. Selection targets a variant's association
  with the full phenotype vector. Real LD, ancestry, kinship and read noise require
  additional models; this API remains explicitly experimental.

## Manhattan plots

- Add `manhattan_data()` and `plot_manhattan()` using optional `qqman` or `CMplot`
  backends. Align coordinates by variant ID and preserve chromosome labels.
- Explicitly choose omnibus or one trait, p-value type, exploratory score or
  signed knockoff W. No p-values or conventional GWAS thresholds are invented.
- Reject zero p-values and missing coordinates by default; explicit omission is
  reported. Draw on the current device, allowing standard PDF/PNG export.

## Corrections

- Prevent integer overflow in AUROC for large positive/negative classes.
- Recover feasible balanced pooling designs after a greedy dead end, using a
  deterministic bounded search. Distinguish exhausted budget from proven
  infeasibility within the recursive axis-aligned partition family.
- Validate the complete data contract before either solver and retain aligned
  data and previous QC exclusions in fitted objects.
- Avoid duplicate preparation in the unpenalized sparse-group shortcut to ridge.

The source package contains code, help, small tests and documentation. Large
simulations, datasets and research benchmark results are not part of the package.
