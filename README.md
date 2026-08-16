# Diver (*Gavia* spp.) displacement around offshore wind farms in the German Bight

Code and input data for the joint-likelihood Bayesian LGCP analysis of diver
redistribution around offshore wind farms (OWFs) in the German North Sea,
comparing a **before** period (2001–2008) and an **after** period (2017–2021).

Manuscript: JEMA-D-25-08275, *Journal of Environmental Management*.

---

## Contents

```
R/joint_likelihood_analysis.R       Full analysis
R/impact_distance_and_habitat_loss.R   The 2022 distance calculation, for provenance
R/diagnose_*.R                      Verification scripts (see below)
data/                               Input data
outputs/                            Figures, tables and cached fits (created on run)
```

### Input data

`data/diver_owf_data.RData` holds every input the analysis reads:

| Object | Description |
|---|---|
| `counts_5km` | Detection-corrected diver counts aggregated to the nodes of a 5 km mesh, spring 2001–2021, by year (`phase`) and survey method |
| `owf_polygons` | OWF footprints with construction/operation dates and cluster labels, repeated per year |
| `prediction_pxl` | Prediction grid (150 × 150 at 2.26 × 2.08 km, 6089 cells) |
| `prediction_mask` | Outline of the prediction area |
| `hd_mask` | Diver main concentration area, BMU (2009) |
| `spa_mask` | SPA Eastern German Bight |
| `mesh_5km` | SPDE triangulation used for fitting |

`data/pred_change_975.tif` is the change field of the previous version of the
analysis, archived so that the comparison in the response to reviewers can be
reproduced.

All spatial data are in UTM zone 32N with **kilometre** units
(`+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs`). Coordinates are already
in that system, so the CRS is assigned rather than transformed.

### Verification scripts

| Script | Question it answers |
|---|---|
| `diagnose_vs_original_tif.R` | Does the current measurement code reproduce the earlier result when applied to the earlier change field? |
| `diagnose_distance_surfaces.R` | How much of the difference is `terra::distance` vs `distanceto::distance_raster`? |
| `diagnose_inla_mode.R` | How much is INLA's change of default computational mode from `classic` to `compact`? |

---

## Requirements

- R ≥ 4.3
- [INLA](https://www.r-inla.org/download-install) (not on CRAN)
- `inlabru` ≥ 2.12 (CRAN), `fmesher`, `sf`, `terra`, `ggplot2`, `dplyr`,
  `viridis`, `RColorBrewer`

```r
install.packages(c("inlabru", "fmesher", "sf", "terra", "ggplot2",
                   "dplyr", "viridis", "RColorBrewer"))
install.packages("INLA",
                 repos = c(getOption("repos"),
                           INLA = "https://inla.r-inla-download.org/R/stable"),
                 dep = TRUE)
```

## Running

Open `gavia-owf-displacement.Rproj` (or `setwd()` to the repository root) and
source `R/joint_likelihood_analysis.R`.

The model fits take hours and are cached under `outputs/`; delete
`outputs/fit_*.rds` and `outputs/change_fit.rds` to refit. Predictions are always
recomputed, so a re-run after changing the prediction geometry cannot return a
stale result.

Posterior sampling is seeded through INLA (`seed = SEED` on `predict()` and
`generate()`), not through `set.seed()`, which does not reach
`inla.posterior.sample`. With a non-zero seed inlabru also forces
single-threaded sampling, so results are deterministic.

---

## Model

Counts are modelled with a negative-binomial likelihood and a log link, with
effort (`area`) as exposure. The two periods enter a **single joint
likelihood** sharing one latent spatial field:

```
before:  log E[N] = spde_before(s) + Intercept              , E = area / c_det
after:   log E[N] = spde_before(s) + spde_change(s) + Intercept_after , E = area
```

`spde_change(s)` is therefore the additive log-scale contrast between periods,
estimated within the joint posterior rather than by differencing two separate
fits. Both fields use PC priors on the Matérn range and marginal standard
deviation. `c_det` is the detection correction applied to the conventional
visual surveys of the before period.

## Displacement distances

The change field is summarised through the posterior probability of a decrease
at each location, p(s) = P(δ(s) < 0), where δ(s) = log₁₀ exp(spde_change(s)).
For a threshold *t*, the affected zone is A(t) = { s : p(s) ≥ t }, so that

- *t* = 0.975 is the 95% significance contour, i.e. the upper limit of the 95%
  credible interval lies below zero — a conservative lower bound on the extent
  of the effect;
- *t* = 0.50 is the zero contour.

Measurement:

- OWFs assigned to the same cluster are **dissolved** into a single polygon, so
  overlapping influence zones are handled by construction.
- The affected zone is restricted to the **connected component adjacent to that
  cluster**, so distances are never measured against an unrelated area of change
  elsewhere in the domain.
- A distance surface is computed at 1 km resolution giving, for every cell, its
  shortest Euclidean distance to the nearest OWF of the cluster. The reported
  distances are the values of that surface **along the boundary of the affected
  component**: each measurement is one point of the contour, valued by how far it
  lies from the wind farm.
- Distances are summarised by their mean and standard deviation; the full
  distribution is the violin plot of Fig. 6.

The extent is reported across the full range of thresholds rather than at a
single cut, together with the proportion of each contour that runs along the
boundary of the prediction area — a zone that does not close inside the study
area is not a measurement of effect range.

A radial summary of the change field is also produced. It is **not** an effect
range: total abundance in the region was conserved, so the distance at which the
radial average returns to zero marks where local loss ceases to outweigh
redistribution gains, and depends on where suitable receiving habitat lies.

---

## Outputs

| File | Content |
|---|---|
| `fig3_density_2001_2008.png`, `fig4_density_2017_2021.png` | Posterior mean density, both periods, common scale |
| `fig5_change.png` | Log₁₀ change with evidence isolines for decrease and increase |
| `fig6_effect_distance.png` | Effect distance by cluster |
| `fig7_zero_contour_posterior.png` | Posterior of the zero-contour distance |
| `fig8_distance_by_threshold.png` | Extent against evidential requirement |
| `fig9_radial_effect_profile.png` | Radial summary of the change field |
| `effect_distances_*.csv`, `effect_distance_tests.csv` | Distances, summaries and the cluster comparison |
| `effect_distance_by_threshold.csv` | Extent and area at each threshold |
| `zero_contour_posterior_*.csv` | Posterior draws, summary and North–South contrast |
| `effect_range_summary.csv`, `radial_profile.csv` | Radial summary |
| `displaced_individuals.csv` | Abundance change within each affected zone |
| `comparison_vs_original_tif.csv`, `distance_surface_comparison.csv`, `inla_mode_comparison.csv` | Verification against the previous version |
| `sessionInfo.txt` | Package versions for the reproducibility statement |

---

## Licence

| | |
|---|---|
| Code (`R/`), documentation, project files | [MIT](LICENSE) |
| Data (`data/`) | [CC BY 4.0](data/LICENSE-DATA.md) |

## Data availability

The aerial survey counts are aggregated to the nodes of the spatial mesh and
were first released with:

> Vilela, R., Burger, C., Diederichs, A., Bachl, F.E., Szostek, L., Freund, A.,
> Braasch, A., Bellebaum, J., Beckers, B., Piper, W. and Nehls, G. (2021).
> Use of an INLA latent Gaussian modeling approach to assess bird population
> changes due to the development of offshore wind farms.
> *Frontiers in Marine Science*, 8, 701332.
> <https://doi.org/10.3389/fmars.2021.701332>

Please cite that paper alongside the present manuscript when using these data.
Full provenance is in [`data/LICENSE-DATA.md`](data/LICENSE-DATA.md).

<!-- TODO: add the DOI of the present paper once published. -->

## Contact

Raúl Vilela — <rvp@duck.com>
