# Diver (*Gavia* spp.) displacement around offshore wind farms in the German Bight

Code and input data for the joint-likelihood Bayesian LGCP analysis of diver
redistribution around offshore wind farms (OWFs) in the German North Sea,
comparing a **before** period (2001–2008) and an **after** period (2017–2021).

Manuscript: JEMA-D-25-08275, *Journal of Environmental Management*.

---

## Contents

```
R/joint_likelihood_before_after.R   Full analysis, Sections 0–12
data/                               Input data (see below)
outputs/                            All figures, tables and cached fits (created on run)
```

### Input data

| File | Description |
|---|---|
| `ips_01_21_noFN10.rds` | Detection-corrected diver counts aggregated to a 5 km grid, spring 2001–2021, by year (`phase`) and survey method |
| `mesh_5km.rds` | SPDE triangulation of the study domain (~5 km resolution) |
| `farms18_spring_byyear_updt.rds` | OWF polygons with construction/operation dates and cluster labels |
| `mask_zee_full.*` | German EEZ mask used to clip the prediction grid |

All spatial data are in UTM zone 32N with **kilometre** units
(`+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs`).

`prediction_pixels.rds` (the original prediction grid) is not included; the
script rebuilds an equivalent grid from the mesh with `fmesher::fm_pixels()`,
clipped to the EEZ mask.

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
source `R/joint_likelihood_before_after.R`.

Sections 1–5 load and check the data and print diagnostics — run these first.
Section 6 fits the models and takes hours; fits are cached under `outputs/`, so
subsequent runs are fast. Set `REFRESH_FITS <- TRUE` to force a refit.

All posterior sampling is seeded (`SEED`), so results are reproducible.

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

Two distances are derived from the posterior of
δ(s) = log₁₀ exp(spde_change(s)):

1. **95% significance contour** — boundary of the region where
   `q0.975(δ(s)) < 0`, i.e. at least 97.5% posterior probability of a
   decrease. A conservative lower bound on the extent of the effect.
2. **Zero-effect contour** — boundary of the region where the change is
   negative, i.e. the transition between negative and neutral/positive change.
   A central estimate, comparable with displacement distances reported in
   earlier studies.

Measurement procedure (Section 7 of the script):

- OWFs assigned to the same cluster are **dissolved** into a single polygon, so
  overlapping influence zones are handled by construction.
- The affected region is restricted to the **connected component** adjacent to
  that cluster, so distances are never measured against an unrelated effect
  zone elsewhere in the domain.
- Measurement points are placed every 0.5 km along the dissolved cluster
  outline; each distance is the shortest Euclidean distance from that point to
  the contour.
- Points that do not lie inside the affected region are recorded as unaffected
  and excluded from the mean; their proportion is reported, which makes the
  asymmetry of the effect explicit.

**Uncertainty (Section 8).** For the zero-effect contour the entire measurement
is repeated inside each posterior draw, so the credible interval comes from the
posterior distribution of the distance rather than from a standard error over
measurement points along a single contour (those points are 0.5 km apart on a
smooth field and are not independent). The North-vs-South comparison is the
paired difference within each draw, reported as P(d_north > d_south).

The 95% significance contour is a property of the posterior, not of a single
draw, and so has no per-draw analogue; it is reported without a resampled
interval.

---

## Outputs

| File | Content |
|---|---|
| `pred_0108.png`, `pred_1721.png` | Posterior mean density, both periods, common scale |
| `int_change_0121.png` | Log₁₀ change with significance and zero contours |
| `fig6_displacement_distances_2026.png` | Distances to both contours, by cluster |
| `fig6b_contour_geometry_2026.png` | Measurement geometry |
| `fig7_zero_contour_posterior_2026.png` | Posterior of the zero-contour distance |
| `displacement_distances_*.csv` | Per-point distances and summaries |
| `zero_contour_posterior_*.csv` | Posterior draws, summary and North–South contrast |
| `displaced_individuals_2026.csv` | Abundance change within each affected zone |
| `grid_resolution_sensitivity_2026.csv` | Sensitivity to prediction grid resolution |
| `sessionInfo_2026.txt` | Package versions for the reproducibility statement |

---

## Data availability and licence

<!-- TODO before making the repository public:
     - confirm that the survey data in data/ may be redistributed
     - add a licence (e.g. MIT or Apache-2.0 for the code, CC-BY for the data)
     - add the DOI once the paper is published -->

## Contact

Raúl Vilela — <rvp@duck.com>
