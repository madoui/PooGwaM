# PooGwaM 2.1.0

Package R pour l'analyse de fréquences alléliques et de plusieurs traits mesurés dans des pools. Il fournit une ridge, un modèle sparse-group multitraits, une calibration explicite, des knockoffs expérimentaux et des Manhattan plots. Le package contient le code, l'aide et les tests légers ; les simulations et benchmarks de recherche restent séparés.

## Installation

```r
# Depuis une archive locale de la version 2.1.0 :
install.packages("PooGwaM_2.1.0.tar.gz", repos = NULL, type = "source")

# Pour les graphiques uniquement :
install.packages(c("qqman", "CMplot"))
```

R >= 4.1.0 est requis. Le noyau utilise les packages standards `stats`, `tools` et `graphics` ; `qqman` et `CMplot` sont optionnels. Le dépôt GitHub du package est [madoui/PooGwaM](https://github.com/madoui/PooGwaM). La V1 et la V2 portent le même nom : utiliser des bibliothèques distinctes pour conserver plusieurs versions.

## Ajuster les modèles

```r
library(PooGwaM)
stopifnot(packageVersion("PooGwaM") == "2.1.0")

pool_ids <- paste0("pool", 1:4)
F <- matrix(c(.1, .3, .7, .9, .2, .4, .5, .8), 4, 2,
            dimnames = list(pool_ids, c("snp1", "snp2")))
Y <- cbind(height = c(-2, -1, 1, 2), mass = c(3, 1, 4, 2))
rownames(Y) <- pool_ids
sizes <- setNames(c(20, 25, 20, 30), pool_ids)
variants <- data.frame(variant_id = colnames(F),
                       chromosome = c("1", "2"), position = c(100, 200))
dat <- as_pool_gwas_data(F, Y, sizes, variants = variants)
ridge <- poogwam(dat, method = "ridge", lambda = 1)
joint <- poogwam(dat, method = "multitrait",
                 lambda_ridge = 1, lambda_group = .01, lambda_l1 = .002)
head(association_table(ridge))
```

Les effets sont exprimés par allèle alternatif supplémentaire. Sans calibration, les résultats sont des scores exploratoires et les p-valeurs sont manquantes. Les données sont alignées par identifiants ; les fréquences et les comptes éventuels doivent être cohérents. Les ajustements revalident les objets, même s'ils ont été modifiés après construction.

## maxT brut ou studentisé

Exemple autonome sur de petites données illustratives :

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

La studentisation estime une moyenne et une échelle omnibus propres à chaque SNP sur un pilote indépendant, puis fige cette transformation pour les permutations principales. Elle vise à réduire les différences d'échelle nulle entre SNPs. La taille du pilote ci-dessus est illustrative ; ce n'est pas un réglage validé universellement.

**`maxT="raw"` reste le défaut de V2.1.** La promotion du mode studentisé attend l'étude confirmatoire de niveau et de puissance prévue au plan. Les deux modes gardent la portée **`global_null_only`** : la studentisation ne prouve ni contrôle fort de FWER sous alternatives partielles, ni FDR par variant. Les p-valeurs marginales et les maxima bruts restent disponibles. B permutations donnent une résolution minimale de `1/(B+1)`.

Le pilote et la série principale utilisent des tirages avec remise, sans éliminer les doublons. Un pilote dégénéré produit un diagnostic, jamais une preuve que la variance nulle est nulle. L'état RNG est restauré après l'appel ; les algorithmes à état caché non restaurable (Box-Muller ou générateurs utilisateur) sont refusés. Pour rejouer dans le même environnement, réutiliser les entrées, les paramètres, la graine et `rng_kind = cal$rng$kind` ; le plan réalisé est aussi enregistré.

`calibrate_pool_counts()` reste un bootstrap expérimental de comptes sous fréquences homogènes, avec maxT brut. Il ne traite pas automatiquement LD, parenté, ascendance, surdispersion ou erreurs de séquençage.

## Knockoffs poolés expérimentaux

```r
# f est ici connu, car il a servi à générer G ci-dessus.
# Il ne doit pas être remplacé sans justification par la fréquence estimée sur G.
G_tilde <- generate_independent_knockoffs(G, frequencies = f, seed = 456)
selection <- select_pool_knockoffs(G, G_tilde, Y, builder,
  q = .1, method = "ridge", fit_args = list(lambda = 1), seed = 123)
selection$selected
selection$W
selection$inference_status

# Le solveur multitraits est également disponible :
# selection <- select_pool_knockoffs(G, G_tilde, Y, builder,
#   method = "multitrait", fit_args = list(lambda_group = .01))
```

Les originaux et copies sont poolés puis ajustés ensemble. Le contraste W compare les normes de coefficients standardisés ; le seuil knockoff+ utilise la correction `+1`. La découverte porte sur un variant associé au vecteur des traits, pas sur chaque couple variant-trait. Aucune p-valeur n'est attribuée à W.

**Statut `experimental_knockoff`.** Le générateur fourni est une référence sous modèle binomial indépendant connu. Il ne génère pas des knockoffs adaptés à un LD arbitraire. Une paire externe peut être fournie, avec provenance et métadonnées alléliques, mais la validation des matrices ne certifie pas son échangeabilité. Le contrôle théorique exige les hypothèses Model-X et la symétrie de toute la procédure. Les HMM estimés, covariables, parenté et fréquences issues de reads demandent des extensions distinctes.

## Manhattan plots avec qqman ou CMplot

```r
coords <- data.frame(variant_id = snps,
  chromosome = rep(c("1", "2"), each = 20),
  position = rep(seq(1000, 20000, by = 1000), 2))

# Omnibus, p-valeurs maxT ajustées, seuil explicite :
plot_manhattan(cal, variants = coords, value = "p_maxT",
                backend = "qqman", threshold = .05)

# Une p-valeur marginale par trait :
plot_manhattan(cal, variants = coords, value = "p_raw",
                test = "trait", trait = "height", backend = "CMplot")

# W signé : pas de transformation en p-valeur.
plot_manhattan(selection, variants = coords, value = "W")

# Préparer les données sans aucune dépendance graphique :
prepared <- manhattan_data(cal, variants = coords)
head(prepared$data)

# Export PDF ; png() est également utilisable :
pdf("manhattan-maxT.pdf", width = 10, height = 5)
plot_manhattan(cal, variants = coords, threshold = .05)
dev.off()
```

Les coordonnées sont jointes par `variant_id`. Les noms de colonnes sont personnalisables via `chromosome` et `position`, et l'ordre via `chromosome_order`. Pour un seul chromosome, utiliser `qqman` : l'axe mono-chromosome de CMplot présente des limites sur les petits intervalles et les étiquettes. Le graphique n'ajoute aucun seuil humain conventionnel. Les p-valeurs nulles sont refusées, sans plancher silencieux. Les valeurs manquantes sont refusées par défaut ; `na_action="omit"` les exclut explicitement avec avertissement et liste des identifiants.

## Corrections et limites

La V2.1 corrige le débordement entier de l'AUROC, les impasses du pooling glouton lorsqu'une autre partition récursive est faisable, et la validation d'objets modifiés avant ajustement. La recherche de pools est bornée par `search_budget` ; un budget épuisé ne signifie pas impossibilité démontrée.

La matrice pools × variants reste en mémoire. Les scores multitraits exploratoires, la calibration par permutation et les knockoffs gardent leurs portées propres. Les allèles doivent conserver leur orientation ; aucune méthode ne répare automatiquement une confusion génétique ou un modèle de reads mal spécifié.

Voir [NEWS](NEWS.md), la [note méthodologique](inst/doc/inference-v2.1.md) et les pages d'aide R pour les détails. La suite de tests compare les solveurs à des références indépendantes et couvre les invariances, les cas limites, les p-valeurs, les échanges knockoffs et les deux moteurs graphiques. Un contrôle logiciel ne remplace pas la validation biologique ou statistique confirmatoire.
