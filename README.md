# PooGwaM 2.0.0

Package R consacré à l'analyse de fréquences alléliques et de plusieurs traits au niveau des pools. Ce dépôt contient uniquement le package, son aide et ses tests. Les simulations, benchmarks, données et runtimes du projet de recherche ne sont pas inclus. La V1 reste accessible dans l'historique Git.

## Installation

Depuis R (version 4.1.0 ou ultérieure) :

```r
install.packages("remotes") # si nécessaire
remotes::install_github("madoui/PooGwaM", upgrade = "never")
```

La V2 remplace la V1 dans la bibliothèque choisie : utiliser des bibliothèques distinctes pour conserver les deux installations. Les deux versions portent le même nom de package R.

## Exemple minimal

Cet exemple utilise de petites matrices illustratives définies directement dans R ; aucun fichier de données externe n'est nécessaire.

```r
library(PooGwaM)
stopifnot(packageVersion("PooGwaM") == "2.0.0")

pool_ids <- paste0("pool", 1:4)
F <- matrix(c(.1, .3, .7, .9, .2, .4, .5, .8), 4, 2,
            dimnames = list(pool_ids, c("snp1", "snp2")))
Y <- cbind(height = c(-2, -1, 1, 2), mass = c(3, 1, 4, 2))
rownames(Y) <- pool_ids
sizes <- setNames(c(20, 25, 20, 30), pool_ids)
dat <- as_pool_gwas_data(F, Y, sizes)

ridge <- poogwam(dat, method = "ridge", lambda = 1)
joint <- poogwam(dat, method = "multitrait",
                lambda_ridge = 1, lambda_group = .01, lambda_l1 = .002)
dim(ridge$beta)                    # 2 SNPs x 2 traits
ridge$diagnostics
head(association_table(ridge))    # Scores exploratoires ; p-valeurs NA
predict(ridge)                    # Moyennes ajustées des pools
```

Les paramètres ci-dessus illustrent l'API ; ils ne sont pas annoncés comme optimaux. Pour évaluer un classement avec une vérité connue par trait, consulter `?evaluate_associations`.

## Capacités

- Jointures explicites par identifiants d'individus, pools, SNPs et traits ; erreurs informatives en cas de doublons, manquants ou valeurs incompatibles.
- Génotypes diploïdes et dosage `2F`, intercept et covariables non pénalisés, effets restitués dans les unités originales des traits.
- Ridge déterministe par système K×K, assemblé par blocs de SNPs. La matrice d'entrée K×P reste en mémoire ; le moteur ne fournit pas encore un stockage entièrement sur disque.
- Estimation sparse-group multi-traits : sélection commune et composante spécifique aux traits, avec diagnostics de convergence/KKT. Ce sont des sélections exploratoires.
- Pools par quantiles ou partition récursive à nombre total fixé, tailles minimales explicites, traitement déterministe des ex æquo.
- Comptages explicitement typés comme chromosomes biologiques ou reads ; les fréquences et profondeurs sont vérifiées.
- Évaluation par trait, vérité correctement alignée, AUROC et average precision sans dépendances lourdes. Le nul complet et l'absence de sélection sont définis. Les corrections multiples ne sont jamais appliquées une seconde fois à des p-valeurs déjà corrigées.

## Calibration et domaine de validité

`calibrate_pool_permutation()` s'applique aux **génotypes individuels observés**, sous échangeabilité des individus ou à l'intérieur de strates déclarées. Les vecteurs multivariés de traits restent ensemble. Le plan phénotypique est fixé et les génotypes sont réassignés aux individus avant de recalculer les fréquences des pools. Le résultat contient p-valeurs Monte Carlo, résolution `1/(B+1)`, score omnibus et maxT sous nul global.

```r
# Exemple d'interface, potentiellement coûteux :
# cal <- calibrate_pool_permutation(G, Y,
#   pool_builder = function(y) design_pools(y, k = c(4,4), min_size = 2),
#   B = 999, seed = 123, fit_args = list(lambda = 1))
# table <- association_table(cal$fit, cal)
```

Il s'agit d'une calibration du **nul global d'indépendance**, pas d'une garantie d'identifier causalement chaque SNP en présence d'effets ailleurs dans le génome. Le maxT n'est pas annoncé comme contrôle fort de la FWER et BH n'est pas annoncé comme garantissant la FDR sous toutes les alternatives. La règle de pooling doit utiliser les phénotypes seuls et les paramètres doivent être fixés indépendamment du résultat.

`calibrate_pool_counts()` est un **bootstrap expérimental** sous fréquences homogènes, individus diploïdes non apparentés et contributions ADN égales. Il utilise un étage pour les chromosomes et deux pour les reads. Il refuse les covariables ; il ne modélise pas LD, parenté, ascendance, erreurs de séquençage ou surdispersion. Ses résultats portent ce statut, sans garantie de calibration finie. Il ne doit pas être utilisé comme solution générale à la structure de population.

La perte quadratique utilise les tailles des pools comme poids d'agrégation. Ce choix n'est pas une modélisation de toutes les incertitudes de mesure Pool-seq. Les variances d'imputation dans la pénalité sont des choix de régularisation ; elles ne prouvent pas Hardy-Weinberg au sein de pools sélectionnés sur les traits.

## Données manquantes et invariants

Les données manquantes sont rejetées par défaut. `missing="drop_variants"` permet d'exclure tout SNP incomplet, avec motifs conservés. Aucune imputation silencieuse ne modifie les dénominateurs. Les SNPs sans variation entre pools et les traits constants ont des coefficients nuls et des diagnostics explicites.

Les allèles doivent garder la même orientation dans tous les pools. La MAF est un filtre global ; jamais un retournement indépendant de chaque fréquence locale. Les covariables doivent être numériques et de plein rang avec l'intercept. Leur présence dans une régression n'établit pas à elle seule l'absence de confusion génétique.

## Migration depuis la V1

| V1 | V2 |
|---|---|
| `poogwam(Method, Freq, phenotypes, clusters)` retourne des p-valeurs empiriques | `poogwam(data, method="ridge")` retourne un modèle ; calibration explicite et séparée. |
| `adaptRidge`, `singValRidge` | `fit_pool_ridge` ; un seul objectif documenté et solveur exact. |
| Réponse multi-traits aplatie | Matrices SNP×trait conservées dans tout le calcul. |
| `pvalues.from.test.stat` | Retirée ; aucune loi normale ajustée automatiquement entre SNPs. |
| `PhenoSim` et `h2.estimation` historique | Retirés du cœur ; import de données et de leur provenance. `sim1` reste une référence historique, pas une simulation corrigée artificiellement. |
| `benchPooGwaM` | `evaluate_associations`, avec type d'entrée obligatoire. |
| `GWalpha`, `Chisquare`, `statgenGWASwrap`, `henshall`, `popGWAS` | Comparateurs hors du package ; la copie GWAlpha R n'est pas réutilisée. |
| Pooling implicite | `design_pools`, `pool_from_individuals`, données identifiées et contrôlées. |

Les helpers `quantile_clustering` et `compute_group_MAFs` restent disponibles, avec identifiants explicites. Cette V2 assume une rupture d'API pour éviter des sorties statistiquement ambiguës. La V1 reste disponible pour reproduire les anciennes analyses.

## Vérification du package

Depuis un terminal à la racine d'une copie de ce dépôt, avec R dans le PATH :

```sh
R CMD build .
R CMD check --no-manual PooGwaM_2.0.0.tar.gz
```

Le package de calcul dépend seulement de R et `stats`. Les tests nécessitent `testthat` (version 3 ou ultérieure). Ils couvrent des références primales indépendantes, invariances, solutions analytiques sparse-group, alignements, cas limites, calibration et erreurs d'API. Le contrôle d'archive ne remplace pas une validation biologique ou une étude de puissance multi-scénarios. Les comparateurs, dont GWAlpha, sont externes au package.

## Points connus de la version 2.0.0

L'audit du 13 octobre 2026 a confirmé deux cas non couverts par les 226 assertions existantes : l'AUROC peut devenir `NA` lorsque le produit des effectifs positifs et négatifs dépasse la capacité d'un entier R, et le pooling équilibré peut rejeter certains plans réalisables avec ex æquo. Ces corrections ne sont pas encore intégrées à cette version.

Après modification manuelle d'un objet de données, appeler `validate_pool_data()` avant l'ajustement pour revérifier notamment les comptes et fréquences. Pour rejouer exactement une calibration, conserver `RNGkind()` et la version de R en plus de la graine. Les garanties et limites statistiques décrites ci-dessus restent applicables.
