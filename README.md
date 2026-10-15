# PooGwaM 2.2.0

Package R pour l'analyse de fréquences alléliques et de plusieurs traits mesurés dans des pools. Il fournit une ridge, un modèle sparse-group multitraits, des cofacteurs fixes décrits par formule, une calibration explicite, des knockoffs expérimentaux et des Manhattan plots. Le package contient le code, l'aide et les tests légers ; les simulations et benchmarks de recherche restent séparés.

## Installation

```r
# Depuis une archive locale de la version 2.2.0 :
install.packages("PooGwaM_2.2.0.tar.gz", repos = NULL, type = "source")

# Pour les graphiques uniquement :
install.packages(c("qqman", "CMplot"))
```

R >= 4.1.0 est requis. Le noyau utilise les packages standards `stats`, `tools` et `graphics` ; `qqman` et `CMplot` sont optionnels. Le dépôt GitHub du package est [madoui/PooGwaM](https://github.com/madoui/PooGwaM). La V1 et la V2 portent le même nom : utiliser des bibliothèques distinctes pour conserver plusieurs versions.

## Comparer ridge, identity et sparse

`ridge` et `multitrait` sont les deux méthodes de l'API `poogwam()`. **`identity` est une option de pénalité**, disponible pour les deux solveurs ; le libellé de benchmark « sparse » désigne un réglage de `fit_pool_multitrait()`. Les trois configurations comparées dans les benchmarks sont :

| Libellé | Appel et réglages | Comportement |
|---|---|---|
| **PooGwaM ridge** | `fit_pool_ridge()`, `lambda=1`, `penalty="within_pool"` | Réduit les coefficients avec une pénalité quadratique propre à chaque SNP, sans imposer de sélection sparse. |
| **PooGwaM identity** | `fit_pool_ridge()`, `lambda=1`, `penalty="identity"` | Même solveur ridge, avec un poids de pénalité quadratique identique pour tous les SNPs. |
| **PooGwaM sparse** | `fit_pool_multitrait()`, `lambda_ridge=1`, `lambda_group=0.01`, `lambda_l1=0.002`, `penalty="within_pool"` | Ajoute des pénalités pouvant annuler exactement des coefficients. Ce libellé de benchmark utilise un seul trait à la fois. |

Avec les fréquences alléliques $F_{kj}$, les poids de pools $w_k=n_k/\sum_k n_k$ et le dosage $X=2F$, la ridge minimise la somme pondérée des carrés des résidus plus $\lambda\sum_{j,t}D_j B_{jt}^2$. L'intercept et les covariables éventuelles ne sont pas pénalisés. Le poids brut est :

$$
D_j^{\mathrm{brut}} =
\begin{cases}
\sum_k w_k\,2F_{kj}(1-F_{kj}) & \text{pour within\_pool},\\
1 & \text{pour identity}.
\end{cases}
$$

Le code applique ensuite `D_j = max(D_j_brut, penalty_floor)`, avec `penalty_floor=1e-8` par défaut. `within_pool` donne généralement un poids plus faible aux variants rares. Cela peut contribuer aux différences d'échelle des coefficients et des scores nuls entre fréquences alléliques. `identity` supprime cette pondération particulière ; il ne garantit pas à lui seul une puissance équilibrée. Les deux options utilisent le même solveur déterministe sur une matrice pools × pools. Voir la [préparation de la pénalité et le solveur ridge](R/fit.R).

Le modèle sparse minimise la moitié du critère quadratique, avec `lambda_ridge`, plus `lambda_group` fois la somme des normes des vecteurs d'effets par variant et `lambda_l1` fois la somme des valeurs absolues des coefficients. Avec plusieurs traits, la pénalité de groupe peut retirer un variant pour tous les traits ; la pénalité L1 permet aussi des zéros spécifiques à un trait. **Avec un seul trait, les deux pénalités sparse s'additionnent : les réglages ci-dessus donnent une pénalité L1 de 0,012**, sans partage d'information entre traits.

Le solveur sparse utilise un gradient proximal accéléré, une recherche de pas et un contrôle de convergence KKT. Avec `lambda_group=0` et `lambda_l1=0`, il appelle directement la ridge exacte avec les autres réglages identiques. Les traits sont standardisés par défaut pour sparse (`standardize_traits=TRUE`), contre `FALSE` pour ridge ; les effets retournés sont toujours reconvertis en unités originales. Voir le [solveur sparse-group](R/multitrait.R).

## Exemples : ridge, identity et sparse

Cet exemple crée les données nécessaires dans R, sans fichier de simulation externe. Les trois premiers ajustements utilisent le même trait ; le dernier illustre séparément le partage entre traits.

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

# Ridge avec pénalité dépendant des fréquences intra-pool :
ridge <- fit_pool_ridge(dat_height, lambda = 1, penalty = "within_pool",
                        standardize_traits = FALSE)

# Même ridge avec pénalité uniforme entre SNPs :
ridge_identity <- fit_pool_ridge(dat_height, lambda = 1, penalty = "identity",
                                 standardize_traits = FALSE)

# Sparse sur un seul trait : pénalité L1 totale de 0.01 + 0.002.
sparse <- fit_pool_multitrait(dat_height, lambda_ridge = 1,
  lambda_group = .01, lambda_l1 = .002, penalty = "within_pool",
  standardize_traits = TRUE)

# Sparse-group conjoint sur height et mass :
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

Les effets sont exprimés par allèle alternatif supplémentaire. Sans calibration, les résultats sont des scores exploratoires et les p-valeurs sont manquantes. Les données sont alignées par identifiants ; les fréquences et les comptes éventuels doivent être cohérents. Les ajustements revalident les objets, même s'ils ont été modifiés après construction.

Ces hyperparamètres illustrent les configurations des benchmarks ; ils ne constituent pas des réglages optimaux universels. Les SNPs à coefficient non nul du modèle sparse ne sont pas automatiquement des découvertes significatives. L'API `poogwam(dat_height, method="ridge", penalty="identity")` est également utilisable ; `method="identity"` et `method="sparse"` ne sont pas des valeurs de méthode valides.

## Ajuster sur le sexe, l'âge, le lot, le site et le traitement

La V2.2 accepte `covariates = ~ sexe + age + lot + site + traitement` et un tableau `metadata` dont les noms de lignes sont les identifiants des individus ou des pools. Les cofacteurs et l'intercept sont des **effets fixes non pénalisés** dans les trois configurations. Ils ne sont ni des effets aléatoires, ni une correction automatique de la parenté ou de la structure génétique.

L'exemple suivant utilise 96 individus fictifs répartis en 16 pools. Les proportions de sexe et les autres cofacteurs varient entre les pools. Il fournit suffisamment de pools pour les sept colonnes du modèle de cofacteurs, intercept compris, et conserve une dimension résiduelle pour l'ajustement génétique.

```r
library(PooGwaM)

individual_ids <- sprintf("person%03d", 1:96)
pool_number <- rep(1:16, each = 6)
within_pool <- rep(1:6, 16)
clusters <- setNames(sprintf("pool%02d", pool_number), individual_ids)

metadata <- data.frame(
  sexe = factor(ifelse(within_pool <= 1 + (pool_number - 1) %% 4, "M", "F"),
                levels = c("F", "M")),
  age = 25 + 1.3 * (pool_number - 1) + within_pool + 2 * (pool_number %% 3),
  lot = factor(paste0("L", 1 + (pool_number - 1) %% 3),
               levels = c("L1", "L2", "L3")),
  site = factor(ifelse((pool_number - 1) %% 8 < 4, "Nord", "Sud"),
                levels = c("Nord", "Sud")),
  traitement = factor(ifelse((pool_number - 1) %% 2 == 0, "temoin", "traite"),
                      levels = c("temoin", "traite")),
  row.names = individual_ids
)

# Petits dosages et traits déterministes pour rendre l'exemple autonome.
G <- vapply(1:4, function(j)
  ((1:96) * (2 * j + 1) + (1:96) %/% (j + 1)) %% 3, numeric(96))
dimnames(G) <- list(individual_ids, paste0("snp", 1:4))
Y <- cbind(
  height = 150 + .6 * metadata$age + 4 * (metadata$sexe == "M") +
    1.2 * (metadata$traitement == "traite") + .8 * G[, 1] + 2 * sin((1:96) / 7),
  mass = 45 + .4 * metadata$age + 3 * (metadata$site == "Sud") +
    1.5 * (metadata$lot == "L3") + G[, 2] + cos((1:96) / 5)
)
rownames(Y) <- individual_ids

dat_cofactors <- pool_from_individuals(G, Y, clusters,
  covariates = ~ sexe + age + lot + site + traitement,
  metadata = metadata)

ridge_cofactors <- fit_pool_ridge(dat_cofactors, lambda = 1,
                                 penalty = "within_pool")
identity_cofactors <- fit_pool_ridge(dat_cofactors, lambda = 1,
                                    penalty = "identity")
sparse_cofactors <- fit_pool_multitrait(dat_cofactors,
  lambda_ridge = 1, lambda_group = .01, lambda_l1 = .002,
  penalty = "within_pool", standardize_traits = TRUE)

ridge_cofactors$covariate_effects  # coefficients par cofacteur et par trait
ridge_cofactors$intercept
ridge_cofactors$diagnostics$covariate_design
stopifnot(all(vapply(list(ridge_cofactors, identity_cofactors, sparse_cofactors),
  function(fit) isTRUE(fit$diagnostics$converged), logical(1))))

# Réencoder la metadata individuelle avec le schéma appris et les mêmes pools.
# newdata contient des fréquences de pools, pas des dosages individuels.
pred <- predict(ridge_cofactors, newdata = dat_cofactors$frequencies,
                metadata = metadata, clusters = clusters)
stopifnot(isTRUE(all.equal(pred, ridge_cofactors$fitted)))

# Un objet de données conserve également les metadata nécessaires :
predict(ridge_cofactors, newdata = dat_cofactors)
```

Le codage est réalisé **avant** l'agrégation individuelle. Avec les contrastes de traitement habituels de R, une colonne indicatrice devient la proportion d'individus de la catégorie correspondante dans chaque pool. Les références de l'exemple sont `F`, `L1`, `Nord` et `temoin` ; définir explicitement les niveaux avant l'ajustement permet de choisir ces références. Le schéma et les contrastes appris sont conservés pour la prédiction, même si une catégorie connue est absente des nouvelles observations. Un niveau inconnu est refusé.

Une formule telle que `~ age + I(age^2) + sexe:traitement` est également encodée au niveau individuel avant les moyennes. Ainsi, la colonne du carré contient la **moyenne des carrés**, et non le carré de l'âge moyen ; l'interaction suit la même règle. Pour des données déjà poolées, `as_pool_gwas_data(..., covariates = ~ ..., metadata = metadata_pools)` applique au contraire la formule directement aux lignes de pools. Des moyennes seules ne permettent pas de reconstruire les transformations individuelles.

Les termes directs `scale(age)` et `poly(age, 2)` conservent les paramètres appris pour la prédiction. Les transformations usuelles `I`, `log`, `log1p`, `sqrt`, `exp` et `abs` sont aussi prises en charge. Les fonctions personnalisées, les variables hors de `metadata` et les imbrications comme `I(scale(age))` ou `poly(scale(age), 2)` sont refusées ; détailler les variables plutôt qu'utiliser `~ .`, et conserver l'intercept de la formule. Les variables logiques gardent les niveaux `FALSE` et `TRUE` pour la prédiction. Un facteur déclaré avec un seul niveau est refusé explicitement.

Les matrices numériques existantes restent acceptées par `covariates` : elles sont considérées comme déjà encodées, sans colonne d'intercept. En prédiction, leurs noms de colonnes doivent correspondre à l'apprentissage. Pour une formule apprise sur des individus, fournir la metadata individuelle et `clusters` ; pour une formule apprise sur des pools, fournir la metadata de ces pools. Toutes les jointures reposent sur les identifiants. Les valeurs manquantes, les individus absents et les catégories inconnues sont des erreurs explicites ; aucun individu n'est écarté silencieusement.

Les objets construits avec une formule conservent la metadata et les appartenances aux pools nécessaires à la validation et à la prédiction, même avec `retain_individuals = FALSE`. Cette option contrôle la conservation des génotypes et phénotypes individuels, pas celle du schéma de cofacteurs.

Le rang doit être vérifié **après pooling** : un lot confondu avec le site, une catégorie absente ou des proportions identiques dans tous les pools peuvent rendre les colonnes redondantes. Le programme refuse alors l'ajustement et indique les colonnes concernées. `diagnostics$covariate_design$residual_df` décrit `K - rang(cbind(1, C))`, avec K pools et C les cofacteurs encodés ; ce n'est pas un nombre de degrés de liberté classique pour calculer des tests t sur la ridge. Le nombre d'individus ne compense pas un nombre insuffisant de pools.

**Cette première prise en charge des cofacteurs concerne les effets et les prédictions.** Les fonctions de calibration et de knockoffs n'ont pas acquis de modèle d'inférence ajusté sur ces cofacteurs. Les scores de `association_table(ridge_cofactors)` restent exploratoires, avec des p-valeurs manquantes sans calibration compatible ; les coefficients non nuls de sparse ne sont pas des découvertes significatives par eux-mêmes.

## Historique des modèles

| Date / version | Évolution du code |
|---|---|
| **24 novembre 2022** — [ba8e7db](https://github.com/madoui/PooGwaM/commit/ba8e7db) | L'ancien `R/selection.R` utilisait LASSO et la sélection de stabilité (`glmnet`, `stabs`), séparément pour deux traits. Cette procédure est distincte du sparse-group actuel. |
| **7 septembre 2023** — [06a68ae](https://github.com/madoui/PooGwaM/commit/06a68ae) | Introduction d'`adaptRidge()` avec mises à jour par gradient stochastique et pénalité dépendant des fréquences ; suppression de l'ancien module de sélection. |
| **14 septembre 2023** — [43f69a9](https://github.com/madoui/PooGwaM/commit/43f69a9) | Introduction du wrapper `poogwam()`, qui appelait `adaptRidge()` puis transformait les coefficients en p-valeurs. |
| **16 octobre 2023** — [a92744c](https://github.com/madoui/PooGwaM/commit/a92744c) | Ajout de `singValRidge()`, fondé sur une SVD et une pénalité commune `mean(alpha)`. Cette pénalité moyenne n'est pas l'actuelle option explicite `identity`. |
| **V2, livraison locale documentée le 25 septembre 2026 ; publication le 13 octobre 2026** — [89dbcef](https://github.com/madoui/PooGwaM/commit/89dbcef) | Introduction de la ridge exacte sur les pools, de l'option `identity` et du nouveau solveur sparse-group. Dosage `2F`, intercept/covariables non pénalisés et séparation entre ajustement et calibration des p-valeurs. |
| **V2.1, 14 octobre 2026** — [0c5c23a](https://github.com/madoui/PooGwaM/commit/0c5c23a) | Formulations ridge et sparse-group conservées. Validation des entrées renforcée, maxT studentisé, knockoffs expérimentaux et Manhattan plots ajoutés. |
| **V2.2.0, 15 octobre 2026** | Ajout de cofacteurs fixes par formule et metadata identifiées, avec encodage avant pooling et réutilisation du schéma appris pour les prédictions. Cette extension concerne l'ajustement et la prédiction ; les portées d'inférence précédentes sont conservées. |

Dans la V1 archivée (0.1.51), `adaptRidge()` utilisait `X=ZF` et des poids `2*colMeans(X)*(1-colMeans(X))`, calculés à partir des fréquences moyennes. La V2 utilise `2F` et la moyenne pondérée des hétérozygoties **par pool** pour `within_pool`. Les coefficients V1/V2 ne sont donc pas directement interchangeables. La V2 a également supprimé la transformation automatique des coefficients en p-valeurs : l'ajustement fournit des effets et des scores, puis une procédure d'inférence explicite est nécessaire.

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

**`maxT="raw"` reste le défaut de V2.2.** La promotion du mode studentisé attend l'étude confirmatoire de niveau et de puissance prévue au plan. Les deux modes gardent la portée **`global_null_only`** : la studentisation ne prouve ni contrôle fort de FWER sous alternatives partielles, ni FDR par variant. Les p-valeurs marginales et les maxima bruts restent disponibles. B permutations donnent une résolution minimale de `1/(B+1)`.

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
