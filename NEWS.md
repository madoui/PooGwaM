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
