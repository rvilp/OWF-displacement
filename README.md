# Diver (*Gavia* spp.) displacement around offshore wind farms in the German Bight

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22117556.svg)](https://doi.org/10.5281/zenodo.22117556)

Code and input data for the joint-likelihood Bayesian LGCP analysis of diver
redistribution around offshore wind farms (OWFs) in the German North Sea,
comparing a **before** period (2001–2008) and an **after** period (2017–2021).

Manuscript: JEMA-D-25-08275, *Journal of Environmental Management*.

---

## Contents

```
R/          Analysis scripts, listed by destination below
data/       Input data
outputs/    Figures, tables and cached fits (created on run)
```

Run `joint_likelihood_analysis.R` first: it writes the cached model fits that
every other script reads.

Note that the density figures are numbered from 3 in `outputs/` and from 4 in
the manuscript, because the schematic was added later as Fig. 3. The table
below gives the correspondence.

### Manuscript

| Script | Produces | Appears as |
|---|---|---|
| `R/fig_method_schematic.R` | `figS_method_schematic.png` | Fig. 3 |
| `R/joint_likelihood_analysis.R` | `fig3_density_2001_2008.png` | Fig. 4 |
| | `fig4_density_2017_2021.png` | Fig. 5 |
| | `fig5_change.png` | Fig. 6 |
| | `fig6_effect_distance.png` | Fig. 7 |
| | `effect_distance_by_threshold.csv` | Table 1 |
| | `effect_distances_summary.csv` | Results §3.4, displacement distances |
| | `displaced_individuals.csv` | Results §3.4, change in abundance |
| | `detection_offsets_applied.csv` | Methods §2.3, detection rates |
| `R/model_validation.R` | `model_validation.csv`, `figS_pit_histogram.png` | Results, model diagnostics |

### Appendix A

| Script | Produces | Appears as |
|---|---|---|
| `R/fig_survey_effort.R` | `figS_survey_effort.png` | Fig. A.1 |
| `R/sensitivity_baseline_period.R` | `sensitivity_baseline_summary.csv`, `figS_baseline_sensitivity.png` | A.2 (noise floor), A.3, Table A.1, Fig. A.2 |
| `R/sensitivity_baseline_effort_control.R` | `sensitivity_effort_vs_years.csv`, `figS_effort_control.png` | A.4, Table A.2, Fig. A.3 |
| `R/sensitivity_after_period.R` | `sensitivity_after_summary.csv` | A.2 (noise floor), A.5, Table A.3 |
| `R/sensitivity_both_periods.R` | `sensitivity_both_periods.csv`, `figS_both_periods_sensitivity.png` | A.6, A.7, Tables A.4–A.5, Fig. A.4 |

The Monte Carlo noise floor in A.2 has no script of its own: it comes from the
control rows, in which the full data are refitted under a different sampling
seed, inside the baseline and after-period series.

### Verification, cited in the response to reviewers

| Script | Question it answers |
|---|---|
| `R/diagnose_vs_original_tif.R` | Does the current measurement code reproduce the earlier result when applied to the earlier change field? |
| `R/diagnose_distance_surfaces.R` | How much of the difference is `terra::distance` versus `distanceto::distance_raster`? |
| `R/diagnose_inla_mode.R` | How much is INLA's change of default computational mode from `classic` to `compact`? |
| `R/diagnose_sampling_reproducibility.R` | Is the posterior sampling deterministic once the seed is fixed? |

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

---

## Requirements

- R ≥ 4.3
- [INLA](https://www.r-inla.org/download-install) (not on CRAN)
- `inlabru` ≥ 2.12 (CRAN), `fmesher`, `sf`, `terra`, `ggplot2`, `dplyr`,
  `viridis`, `RColorBrewer`

The results reported in the manuscript were produced with R 4.5.2, INLA 25.10.19,
inlabru 2.13.0, fmesher 0.7.0, sf 1.1-0, terra 1.8-93, ggplot2 4.0.2 and
dplyr 1.2.1. The full session is written to `outputs/sessionInfo.txt` on every
run.

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

The model fits are cached under `outputs/`, although they are not expensive:
the mesh has 1000 nodes, so the joint latent field is roughly 2000-dimensional
with a sparse precision over 2666 observations, and INLA fits it in seconds.
A full run, predictions and figures included, is a matter of minutes. The cache file names
carry `MODEL_TAG`, which is bumped whenever the model changes, so a formula
change cannot silently pick up an old fit; delete
`outputs/fit_*_<tag>.rds` and `outputs/change_fit_<tag>.rds` to force a refit
anyway. Predictions are always recomputed, so a re-run after changing the
prediction geometry cannot return a stale result.

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
before:  log E[N] = spde_before(s) + Intercept                       , E = area × d(m)
after:   log E[N] = spde_before(s) + spde_change(s) + Intercept_after , E = area × d(m)
```

`spde_change(s)` is therefore the additive log-scale contrast between periods,
estimated within the joint posterior rather than by differencing two separate
fits. Both fields use PC priors on the Matérn range and marginal standard
deviation.

### Survey technique

`d(m)` is the detection rate of the technique that produced each observation,
so that exposure is the area searched scaled by the fraction of birds present
that the technique would have recorded. HiDef is the reference at 1.

| Technique | log effect vs HiDef | `d(m)` |
|---|---|---|
| HiDef | 0 | 1.000 |
| DAISI | −0.10 | 0.905 |
| APEM | −0.40 | 0.670 |
| Conventional (visual) | −0.22 | 0.803 |

The log effects are the medians of Table A-1 of Vilela et al. (2021); the
visual figure is what remains after the distance-sampling correction already
applied in `NHAT`.

They enter as a **known offset, not as a fitted coefficient**, because a fitted
technique term is not identifiable in this design. The after period aggregates
five years without a temporal term, and the techniques are segregated in both
time and space: APEM and the visual surveys were flown only in 2017–2018; DAISI
covers the north and the visual surveys the south, so the two never sample the
same mesh node; only 14 of APEM's 130 observations fall on a node-year that
HiDef also covered. Fitting the term returns technique confounded with year and
region, and reverses the sign of the DAISI and APEM effects relative to the
published values. Vilela et al. (2021) estimated them in a year-by-year
spatiotemporal model, where technique is separated from when and where each
technique was flown, so the estimate is taken from there.

## Displacement distances

The change field is summarised through the posterior probability of a decrease
at each location, p(s) = P(δ(s) < 0), where δ(s) = log₁₀ exp(spde_change(s)).
For a threshold *t*, the affected zone is A(t) = { s : p(s) ≥ t }, so that

- *t* = 0.975 is the **high-evidence zone**: the whole 95% credible interval for
  the change lies below zero — a conservative lower bound on the extent of the
  effect;
- *t* = 0.50 is the **zero contour**.

Thresholds are written as probabilities (0.975, 0.95, 0.90) and never as
percentages, because "95%" would then refer both to a relaxed threshold and to
the width of the credible interval that defines the 0.975 zone. The credible
interval appears only once, in the definition above.

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
| `detection_offsets_applied.csv` | Detection rate applied to each technique, with effort by period |
| `comparison_vs_original_tif.csv`, `distance_surface_comparison.csv`, `inla_mode_comparison.csv` | Verification against the previous version |
| `sessionInfo.txt` | Package versions for the reproducibility statement |

---

## Licence

| | |
|---|---|
| Code (`R/`), documentation, project files | [MIT](LICENSE) |
| Data (`data/`) | [CC BY 4.0](data/LICENSE-DATA.md) |

## Data availability

The aerial survey counts are aggregated to the nodes of the spatial mesh. The
spring 2001–2018 records were first released with the paper below; the spring
2019–2021 records are released here for the first time, under the same licence.

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
