# Inference in PooGwaM 2.1

## Studentized maxT

Let s_j be the norm of the standardized phenotype-effect vector for variant j.
The standardization across phenotype units is distinct from the null
normalization across variants introduced in this release.

The independent pilot draws B0 uniform permutations with replacement from the
allowed group, including its identity. With pilot means m_j and positive scales
a_j, define t_j(s) = (s - m_j) / a_j. The main series contains B new draws from
the same group; it does not remove repeated permutations. Its maxima are
M_b = max_j t_j(s_bj). The adjusted probabilities are
(1 + sum_b I[M_b >= t_j(s_observed,j)]) / (B + 1), with conservative numerical
handling of ties.

**Argument for the proposed construction.** Under global-null exchangeability,
conditional on the permutation orbit and fixed phenotypes, the observation and
independent uniformly permuted datasets have the same conditional distribution.
The independent pilot is conditionally independent of the observation and the
main series. Conditioning further on the pilot fixes all transformations. The
maximum of the transformed observed scores can thus be ranked with the main
maxima. This bounds the probability of any rejection under the global null.
This is an application of the conditional permutation argument, rather than a
new guarantee for individual conditional null hypotheses.
[Hemerik and Goeman (2018)](https://arxiv.org/html/1411.7565).

A pilot SD of zero does not establish constancy over the orbit. The implementation
uses a positive floor fixed from the pilot and records affected variants. If all
pilot scores are zero, a documented fixed unit-scale fallback is used with a
warning. No permutation failure is discarded. The main-series null summaries
retain their original meaning; pilot summaries are stored separately.

The legacy main-series SD summaries can underflow or overflow at extreme
floating-point score scales. They are diagnostic only and are not used to
studentize; the pilot has separate scale protection. Finite coefficients alone
do not guarantee useful numerical accuracy for arbitrarily ill-conditioned fits.

Estimating moments from the main null scores alone and then comparing the
observed score with those same normalized null scores would treat observation
and simulated datasets asymmetrically. That construction is not used here.
The alternative symmetric normalization over observation plus all replicates,
minP, and analytic normalizers are deferred.

`raw` remains the default in 2.1. Both modes retain `global_null_only`; neither
establishes strong family-wise control under partial alternatives. Requirements
such as subset pivotality are separate from studentization itself.
[Ge, Dudoit and Speed (2003)](https://genomics.princeton.edu/storeylab/papers/Ge_etal_Test_2003.pdf).
The plug-in count bootstrap retains its experimental status.

## Knockoff statistics after pooling

The target is H0_j: Y independent of G_j conditional on G_-j, where Y may contain
several phenotypes. This is not a causal claim. Model-X assumes pairwise swap
exchangeability of (G, Gtilde) and Gtilde independent of Y conditional on G. Under
these conditions, null-pair swaps also preserve the joint distribution with Y.
[Candes et al. (2018)](https://doi.org/10.1111/rssb.12265).

**Composition argument for this implementation.** Conditional on any independent
pool-builder randomness, the same phenotype-only pooling map is applied to both
members of every pair. This map commutes with pair swaps. The code orders the
unordered dosage pairs canonically, using individual and variant IDs to fix all
orders. Therefore swapping any subset gives the identical augmented matrix to
the numerical solver; only the mapping back to original and knockoff changes.
The difference of the two coefficient-vector norms changes sign on swapped
pairs and leaves all other contrasts unchanged. Identical copies receive W=0.
This supplies the sign-flip property even when solver summation order matters.

This argument applies to the statistic as a function of the original individual
data. It does not assert that aggregated Ftilde automatically satisfies a new
Model-X conditional-independence property given aggregated F.

The knockoff+ threshold is the smallest positive observed |W| such that
(1 + number(W <= -t)) / max(1, number(W >= t)) <= q. If there is no threshold,
there are no selections. The result contains W and selections, not p-values.
With multivariate Y, the discovery unit is a variant associated with the trait
vector, not each variant-trait pair. A sparse-group penalty across traits does
not itself create group knockoffs across loci.

The implemented reference generator samples independent binomial copies using
externally supplied probabilities. It is appropriate only if the original data
follow that known product-binomial model. Structural matrix checks and declared
provenance cannot establish this assumption. The complete API remains labeled
`experimental_knockoff`.

An exact HMM sampler is exact for the specified HMM, not necessarily the real LD
distribution. Model estimation errors can affect FDR; empirical swap diagnostics
are not proof of exchangeability. See [Sesia et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6373422/)
and [Barber et al. (2020)](https://arxiv.org/abs/1801.03896).
Optional HMM adapters, covariates, related individuals, group-locus hypotheses
and a symmetric sequencing-read observation model are outside this release.

## Display and interpretation

`plot_manhattan()` uses [qqman](https://search.r-project.org/CRAN/refmans/qqman/html/manhattan.html)
or [CMplot](https://search.r-project.org/CRAN/refmans/CMplot/html/CMplot-package.html).
The selected p-value column is transformed to -log10; scores and W retain their
original scale. Coordinates are joined by ID. No plotting transformation
changes the inferential scope of the result, and no default significance line
assumes a species or a testing family.
