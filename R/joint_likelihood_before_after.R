###############################################################################
#  Joint-likelihood before/after model - Gavia spp. displacement around OWFs
#  German Bight (JEMA-D-25-08275)
#
#  Original: Raul Vilela, 13.07.2021 / rev. Dec 2022
#            ("joint likelihood before after_dec22.R")
#  This version: 2026-08 working copy.
#
#  Changes relative to the Dec-2022 script:
#    (1) All paths made relative to the script folder. No setwd().
#    (2) Ported from the vendored inlabru 2.1.12.999 to current CRAN inlabru
#        (>= 2.12, tested against 2.14). Retired packages rgdal / rgeos / sp
#        replaced by sf / fmesher / terra.
#    (3) Dead code from the inlabru vignette removed (see NOTE-6 below).
#    (4) New Section 7: displacement distances measured to BOTH the 95%
#        significance contour (conservative lower bound, as in Fig. 6 of the
#        submitted MS) AND the zero-effect contour (central estimate),
#        as requested by Reviewer #3 / handling editor.
#    (5) New Section 8: the zero-effect distance is computed PER POSTERIOR
#        SAMPLE, so its credible interval and the North-vs-South comparison
#        come from the joint posterior instead of a t-test over hundreds of
#        autocorrelated points along a single contour. See REVIEW NOTE-13.
#    (6) Reproducibility: fixed seed, cached model fits, no hidden state.
#
#  HOW TO RUN: open gavia-owf-displacement.Rproj (or setwd() to the repo
#  root) and source this file. Sections 1-5 are cheap and print diagnostic
#  messages; check them before running Section 6, which fits the models and
#  takes hours. Fits are cached under outputs/, so a re-run is fast.
#
#  Tested with: R <version>, INLA <version>, inlabru <version>.
#  (fill in from outputs/sessionInfo_2026.txt after the first full run)
#
#  Assumptions (unchanged):
#    before period: 2001-2008   (phase <= 8)
#    after  period: 2017-2021   (phase > 16)
#    Merkur OWF (construction started summer 2017) sits inside an area already
#    impacted by Trianel, Alpha Ventus and Borkum Riffgrund, so its own
#    construction is not treated as an independent impact.
###############################################################################


# =============================================================================
# 0. Configuration
# =============================================================================

## --- Project layout ------------------------------------------------------
## The working directory must be the repository root. Opening
## gavia-owf-displacement.Rproj in RStudio does this automatically;
## otherwise setwd() to the folder that contains R/, data/ and outputs/.
DATA_DIR <- "data"
OUT_DIR  <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)

if (!dir.exists(DATA_DIR)) {
  stop("Working directory is not the repository root: '", DATA_DIR,
       "' not found. Open gavia-owf-displacement.Rproj or setwd() to the ",
       "repository root.", call. = FALSE)
}
dir.create(OUT_DIR, showWarnings = FALSE)

## --- Packages -----------------------------------------------------------
## inlabru is now loaded from CRAN. The vendored copy in inlabru.zip
## (v2.1.12.999, 2019) is NO LONGER USED - see NOTE-1.
library(INLA)          # >= 24.x
library(inlabru)       # >= 2.12  (CRAN)
library(fmesher)
library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)

stopifnot(utils::packageVersion("inlabru") >= "2.12.0")

sf_use_s2(FALSE)  # projected CRS in km; s2 must stay off

## --- CRS ----------------------------------------------------------------
## The data are stored in UTM32N with KILOMETRE units. The original script
## used "+init=epsg:32632 +proj=utm +zone=32 ... +units=km", which is
## self-contradictory under PROJ >= 6 (epsg:32632 implies metres, then
## overridden). We declare the km-unit CRS directly and ASSIGN it to the
## loaded objects - we do not transform them, because their coordinates are
## already in km-UTM32.
main_crs <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")

## --- Analysis settings ---------------------------------------------------
BEFORE_PHASE_MAX <- 8    # 2001-2008
AFTER_PHASE_MIN  <- 17   # 2017-2021

## Detection correction applied to the "before" (conventional visual) surveys.
## E is divided by this factor, i.e. effective effort is reduced by ~12%,
## equivalent to an average detection probability of 1/1.136364 = 0.88.
## NOTE-3: the Methods section states the mrds correction was staggered by
## group size and sea state; here a single scalar is applied. Either the text
## or the code needs to be reconciled before resubmission.
VISUAL_DETECTION_CORRECTION <- 1.136364

## SPDE PC priors.
## The baseline field and the change field now have SEPARATE prior objects.
## The default values are identical to the original script, so this changes
## nothing by itself - it just makes the prior sensitivity check a two-line
## edit instead of a rewrite. The change field plausibly operates at a
## shorter range than the background distribution; PRIOR_RANGE_CHANGE is the
## knob to test that.
PRIOR_SIGMA_BA     <- c(0.2, 0.01)   # P(sigma > 0.2) = 0.01
PRIOR_RANGE_BA     <- c(15, NA)      # range median 15 km
PRIOR_SIGMA_CHANGE <- PRIOR_SIGMA_BA
PRIOR_RANGE_CHANGE <- PRIOR_RANGE_BA

N_SAMPLES_CHANGE <- 500          # posterior samples for the summary maps
N_SAMPLES_DIST   <- 250          # posterior samples for the distance posterior
PRED_CELL_KM     <- 1.0          # prediction grid resolution (km)

## Reproducibility. Every posterior sampling call is seeded, otherwise the
## numbers in the paper cannot be reproduced even on this machine.
SEED <- 20260811

## Model fits are cached to disk. Set REFRESH_FITS <- TRUE to force a refit.
REFRESH_FITS <- FALSE

## The exploratory spatio-temporal model of Section 3 is not used for any
## result in the manuscript.
RUN_EXPLORATORY <- FALSE

## Search radius for the displacement measurement. Distances are measured
## inside a buffer of this size around each cluster; if the affected zone
## fills the whole buffer the measurement is flagged as right-censored
## instead of being silently reported as MAX_SEARCH_KM.
MAX_SEARCH_KM <- 50

## Cluster definitions used for the displacement distances.
## Overlapping / adjacent OWFs are DISSOLVED into a single cluster polygon
## before measuring distances (Reviewer #3, point ii).
CLUSTER_DEF <- list(
  north = c("West of Sylt", "Butendiek"),        # matched against farms$cluster
  south = c("North of Borkum")                   # <- CHECK against your data
)
## Fallback: if the cluster labels in the shapefile do not match the names used
## in the manuscript, define the south cluster by farm name instead:
SOUTH_FARM_NAMES <- c("Trianel Windpark Borkum", "Merkur Offshore",
                      "Alpha Ventus", "Borkum Riffgrund 1",
                      "Borkum Riffgrund 2")

OWF_BOUNDARY_SPACING_KM <- 0.5   # spacing of measurement points along OWF edge


# =============================================================================
# 1. Load inputs
# =============================================================================

## Plain ggplot2 on plain data frames: avoids any dependence on the inlabru
## gg() S3 methods, whose sf handling has changed across versions.
as_xyz <- function(x) {
  cbind(as.data.frame(sf::st_coordinates(x)), sf::st_drop_geometry(x)) %>%
    dplyr::rename(x = X, y = Y)
}

## Cache an expensive expression on disk. `expr` is lazy: it is only
## evaluated if the cache file is absent or REFRESH_FITS is TRUE.
cached <- function(file, expr, refresh = REFRESH_FITS) {
  if (!refresh && file.exists(file)) {
    message("Loading cached: ", file)
    return(readRDS(file))
  }
  value <- expr
  saveRDS(value, file)
  message("Cached: ", file)
  value
}

## Small helper: fail loudly and early on missing inputs.
read_rds_checked <- function(path) {
  if (!file.exists(path)) {
    stop("Missing input file: ", normalizePath(path, mustWork = FALSE),
         call. = FALSE)
  }
  readRDS(path)
}

## 5 km aggregated count data, 2001-2021, excluding FINO-10 (?) subset.
countdata <- read_rds_checked(p_data("ips_01_21_noFN10.rds"))

## Spatial mesh (built with the old inla.mesh code -> upgrade to fm_mesh_2d).
mesh_raw <- read_rds_checked(p_data("mesh_5km.rds"))
spatial_mesh <- fmesher::fm_as_mesh_2d(mesh_raw)

## OWF polygons, by year/phase.
farms_raw <- read_rds_checked(p_data("farms18_spring_byyear_updt.rds"))

## German EEZ mask (used to clip the prediction grid).
eez <- st_read(p_data("mask_zee_full.shp"), quiet = TRUE)

## --- Convert to sf and harmonise CRS ------------------------------------
## The .rds objects and mask_zee_full.shp are already in km-UTM32, so their
## CRS is ASSIGNED, not transformed. Anything stored in geographic
## coordinates (e.g. ger_EEZ.shp, which is plain lon/lat WGS84) is
## reprojected instead. Never blanket-assign - that was NOTE-9 in the
## original script.
to_sf_km <- function(x) {
  y <- if (inherits(x, "sf")) x else sf::st_as_sf(x)
  crs_in <- sf::st_crs(y)
  if (!is.na(crs_in) && isTRUE(sf::st_is_longlat(y))) {
    y <- sf::st_transform(y, main_crs)
  } else {
    sf::st_crs(y) <- main_crs
  }
  y
}

countdata <- to_sf_km(countdata)
farms     <- to_sf_km(farms_raw)
eez       <- to_sf_km(eez)

## Sanity check: the mesh and the data must live in the same coordinate space.
message("mesh bbox:  ", paste(round(apply(spatial_mesh$loc[, 1:2], 2, range), 1),
                              collapse = " "))
message("data bbox:  ", paste(round(st_bbox(countdata), 1), collapse = " "))


# =============================================================================
# 2. Prepare the modelling data
# =============================================================================

## Drop vertices with no survey coverage.
dat <- countdata %>%
  filter(!is.na(area), !is.na(NHAT), area > 0)

## Method indicators (kept for the exploratory model in Section 3).
dat <- dat %>%
  mutate(
    group_index          = match(phase, sort(unique(phase))),
    indicate_daisi       = as.integer(method == "DAISI"),
    indicate_apem        = as.integer(method == "APEM"),
    indicate_conventional = as.integer(method == "conventional"),
    ## NOTE-8: NHAT is a detection-corrected (non-integer) abundance estimate.
    ## Rounding is required by the nbinomial likelihood but discards part of
    ## the mrds correction, especially for small counts.
    NHAT = round(NHAT)
  )

message("phases present: ", paste(sort(unique(dat$phase)), collapse = ", "))
message("methods present: ", paste(unique(dat$method), collapse = ", "))

## Quick look
ggplot() +
  geom_sf(data = farms, fill = NA, colour = "grey40") +
  geom_sf(data = dat, aes(colour = NHAT / area), size = 0.4) +
  scale_colour_viridis_c() +
  facet_wrap(~ phase + method) +
  coord_sf() +
  theme_bw()


# =============================================================================
# 3. (Optional) Exploratory spatio-temporal model with method effects
#    Not used for the results in the manuscript. Kept for reference only.
# =============================================================================

if (RUN_EXPLORATORY) {

  matern_all <- INLA::inla.spde2.pcmatern(
    spatial_mesh,
    prior.sigma = c(0.1, 0.01),
    prior.range = c(45, NA)
  )

  ## NOTE-5: the original code used
  ##   is_daisi + is_apem + is_conventional + Intercept
  ## which is rank-deficient: the three indicators sum to 1 and are therefore
  ## perfectly collinear with the intercept. One level must be dropped
  ## (here: conventional = reference).
  cmp_all <- NHAT ~
    space_spde(geometry, model = matern_all) +
    is_daisi(indicate_daisi, model = "linear") +
    is_apem(indicate_apem,  model = "linear") +
    Intercept(1)

  fit_all <- bru(
    cmp_all,
    bru_obs("nbinomial", formula = NHAT ~ ., data = dat, E = dat$area),
    options = list(
      verbose = TRUE,
      control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
      control.inla    = list(int.strategy = "ccd")
    )
  )
  summary(fit_all)
}


# =============================================================================
# 4. Before / after subsets
# =============================================================================

before <- dat %>% filter(phase <= BEFORE_PHASE_MAX)
after  <- dat %>% filter(phase >= AFTER_PHASE_MIN)

stopifnot(nrow(before) > 0, nrow(after) > 0)
message("before: ", nrow(before), " rows, phases ",
        paste(range(before$phase), collapse = "-"),
        " | methods: ", paste(unique(before$method), collapse = ", "))
message("after:  ", nrow(after),  " rows, phases ",
        paste(range(after$phase), collapse = "-"),
        " | methods: ", paste(unique(after$method), collapse = ", "))

## OWF footprint for the AFTER period, dissolved by farm (the byyear file
## repeats each polygon once per year).
farms_after <- farms %>%
  filter(phase >= AFTER_PHASE_MIN) %>%
  group_by(Name, cluster) %>%
  summarise(.groups = "drop")

if (nrow(farms_after) == 0) {
  ## Fall back to the last available phase if the phase coding differs.
  farms_after <- farms %>%
    filter(phase == max(phase, na.rm = TRUE)) %>%
    group_by(Name, cluster) %>%
    summarise(.groups = "drop")
}
message("OWFs in the after period: ", nrow(farms_after))
print(sort(unique(as.character(farms_after$cluster))))


# =============================================================================
# 5. Prediction grid
# =============================================================================

## The original script read p_data("prediction_pixels.rds"), which is NOT present in
## this folder (NOTE-2). We rebuild an equivalent grid from the mesh, clipped
## to the EEZ mask. If you still have the original file, uncomment the first
## branch so the figures reproduce cell-for-cell.

if (file.exists(p_data("prediction_pixels.rds"))) {
  pxl_raw <- readRDS(p_data("prediction_pixels.rds"))
  pxl <- if (inherits(pxl_raw, "sf")) pxl_raw else sf::st_as_sf(pxl_raw)
  ## SpatialPixelsDataFrame may come back as polygons; we need cell centres.
  if (!all(as.character(st_geometry_type(pxl)) == "POINT")) {
    pxl <- suppressWarnings(st_centroid(pxl))
  }
  st_crs(pxl) <- main_crs
} else {
  message("prediction_pixels.rds not found - rebuilding grid from the mesh.")
  xr   <- range(spatial_mesh$loc[, 1])
  yr   <- range(spatial_mesh$loc[, 2])
  dims <- c(ceiling(diff(xr) / PRED_CELL_KM),
            ceiling(diff(yr) / PRED_CELL_KM))
  pxl  <- fmesher::fm_pixels(spatial_mesh, dims = dims,
                             mask = eez, format = "sf")
  st_crs(pxl) <- main_crs
}
message("prediction grid: ", nrow(pxl), " cells")

## Cell area (km^2) - needed as the exposure when predicting expected counts.
pxl_coords <- st_coordinates(pxl)
cell_dx <- min(diff(sort(unique(round(pxl_coords[, 1], 6)))))
cell_dy <- min(diff(sort(unique(round(pxl_coords[, 2], 6)))))
cell_area <- cell_dx * cell_dy
message("cell size: ", cell_dx, " x ", cell_dy, " km  (area ", cell_area, " km2)")


# =============================================================================
# 6. Models
# =============================================================================

matern_ba <- INLA::inla.spde2.pcmatern(
  spatial_mesh,
  prior.sigma = PRIOR_SIGMA_BA,
  prior.range = PRIOR_RANGE_BA
)

matern_change <- INLA::inla.spde2.pcmatern(
  spatial_mesh,
  prior.sigma = PRIOR_SIGMA_CHANGE,
  prior.range = PRIOR_RANGE_CHANGE
)

inla_opts <- list(
  verbose = TRUE,
  control.compute   = list(dic = TRUE, cpo = TRUE, waic = TRUE),
  control.inla      = list(int.strategy = "ccd"),
  control.predictor = list(compute = TRUE,
                           quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975))
)

# -----------------------------------------------------------------------------
# 6a. Separate before / after fits (Figs. 3 and 4)
# -----------------------------------------------------------------------------

cmp_b <- ~ space_spde(geometry, model = matern_ba) + Intercept(1)

fit_b <- cached(p_out("fit_before_2026.rds"), bru(
  cmp_b,
  bru_obs("nbinomial",
          formula = NHAT ~ space_spde + Intercept,
          data    = before,
          E       = before$area / VISUAL_DETECTION_CORRECTION),
  options = inla_opts
))

fit_a <- cached(p_out("fit_after_2026.rds"), bru(
  cmp_b,
  bru_obs("nbinomial",
          formula = NHAT ~ space_spde + Intercept,
          data    = after,
          E       = after$area),
  options = inla_opts
))

## --- Posterior mean density (individuals per km^2) ----------------------
## NOTE-4: the original script predicted a SINGLE negative-binomial draw per
## cell (vector_rnbinom) and then applied scale() separately to each period.
## That (a) adds observation noise to what is presented as a density surface,
## and (b) z-standardises the two periods independently, so Figs. 3 and 4 were
## not on a comparable scale despite a shared legend labelled [Ind. km^-2].
## Here we predict the latent intensity directly, on a common scale.

set.seed(SEED)
dens_before <- cached(p_out("dens_before_2026.rds"),
  predict(fit_b, pxl, ~ exp(space_spde + Intercept),
          n.samples = N_SAMPLES_CHANGE))
set.seed(SEED + 1)
dens_after  <- cached(p_out("dens_after_2026.rds"),
  predict(fit_a, pxl, ~ exp(space_spde + Intercept),
          n.samples = N_SAMPLES_CHANGE))

dens_limits <- range(c(dens_before$mean, dens_after$mean), na.rm = TRUE)

df_before <- as_xyz(dens_before)
df_after  <- as_xyz(dens_after)

density_map <- function(df, title, farm_lty) {
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = mean)) +
    scale_fill_viridis(option = "magma",
                       name = expression(paste("[Ind. km"^-2 * "]")),
                       limits = dens_limits) +
    geom_sf(data = farms_after, fill = NA, colour = "red",
            linetype = farm_lty) +
    coord_sf(expand = FALSE) +
    ggtitle(title) +
    theme_bw(base_size = 14)
}

p_before <- density_map(df_before, "2001-2008", 2)
p_after  <- density_map(df_after,  "2017-2021", 1)

ggsave(p_out("pred_0108.png"), p_before, width = 8.3, height = 5, dpi = 300)
ggsave(p_out("pred_1721.png"), p_after,  width = 8.3, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 6b. Joint-likelihood change model
# -----------------------------------------------------------------------------
## Both periods enter a single joint likelihood. spde_before is shared;
## spde_change is the additive spatial contrast (log-scale change).
## NOTE-6: the vignette boilerplate that used to sit here (df1/df2, lik1/lik2,
## jcmp with undefined x/y/z) has been deleted - it referenced variables that
## do not exist in this dataset and would have errored on execution.
## NOTE-7: the original after.likelihood passed `ips = pxl.all`. Integration
## points only apply to family = "cp" (Cox process); with family = "nbinomial"
## they are meaningless and are dropped here.

cmp_joint <- ~
  spde_before(geometry, model = matern_ba) +
  spde_change(geometry, model = matern_change) +
  Intercept(1) +
  Intercept_after(1)

obs_before <- bru_obs(
  "nbinomial",
  data    = before,
  formula = NHAT ~ spde_before + Intercept,
  E       = before$area / VISUAL_DETECTION_CORRECTION
)

obs_after <- bru_obs(
  "nbinomial",
  data    = after,
  formula = NHAT ~ spde_before + spde_change + Intercept_after,
  E       = after$area
)

## The predictor is linear in the latent field, so a single linearisation
## step is exact. Raise bru_max_iter if you add any nonlinear term.
change_fit <- cached(p_out("change_fit_2026.rds"), bru(
  cmp_joint, obs_before, obs_after,
  options = c(inla_opts, list(bru_max_iter = 1))
))

summary(change_fit)

## Posterior change field, log10 scale (log10(exp(spde_change)))
set.seed(SEED + 2)
change_intensity <- cached(p_out("change_intensity_2026.rds"),
  predict(change_fit, pxl,
          ~ log10(exp(spde_change)),
          n.samples = N_SAMPLES_CHANGE))

df_intens <- as_xyz(change_intensity)

## --- Figure 5 -----------------------------------------------------------
p_change <- ggplot() +
  geom_raster(data = df_intens, aes(x = x, y = y, fill = mean)) +
  geom_contour(aes(x = x, y = y, z = q0.975), data = df_intens,
               breaks = 0, colour = "red",  linetype = 2) +   # sig. decrease
  geom_contour(aes(x = x, y = y, z = q0.025), data = df_intens,
               breaks = 0, colour = "blue", linetype = 2) +   # sig. increase
  geom_contour(aes(x = x, y = y, z = mean),  data = df_intens,
               breaks = 0, colour = "grey30", linetype = 1) + # zero contour
  geom_sf(data = farms_after, fill = NA, colour = "black") +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(3, "RdBu"),
    limits  = c(-max(abs(df_intens$mean), na.rm = TRUE),
                 max(abs(df_intens$mean), na.rm = TRUE))
  ) +
  coord_sf() +
  ggtitle("Log10 change (2017-2021 vs 2001-2008)") +
  theme_minimal(base_size = 14)

ggsave(p_out("int_change_0121.png"), p_change, width = 8.3, height = 5, dpi = 300)


###############################################################################
# =============================================================================
# 7. DISPLACEMENT DISTANCES
#    (a) 95% significance contour  -> conservative lower bound (Fig. 6 as
#        submitted)
#    (b) ZERO-EFFECT contour       -> central estimate, comparable with
#        earlier displacement studies  [NEW, Reviewer #3]
# =============================================================================
###############################################################################

## Formal definition of the procedure (answers Reviewer #3 points i-iii):
##
##  (i)  THRESHOLD. The affected region A is the set of prediction cells
##       satisfying a criterion on the posterior of the change field
##       delta(s) = log10 exp(spde_change(s)):
##         - significance criterion : q0.975(delta(s)) < 0
##           i.e. at least 97.5% posterior probability of a decrease;
##         - zero-effect criterion  : E[delta(s)] < 0
##           i.e. the posterior mean change is negative.
##
##  (ii) MULTIPLE / OVERLAPPING OWFs. All OWFs assigned to a cluster are
##       dissolved (union) into a single polygon before measuring, so
##       overlapping influence zones are handled by construction. The
##       affected region is then restricted to the connected component of A
##       that intersects that cluster - distances are never measured to an
##       unrelated effect zone elsewhere in the domain.
##
## (iii) AGGREGATION. Measurement points are placed at regular
##       OWF_BOUNDARY_SPACING_KM intervals along the dissolved cluster
##       outline. For each point p the distance is the shortest Euclidean
##       distance from p to the boundary of the cluster's affected region.
##       Points that do not lie inside the affected region are recorded as
##       "not affected" and excluded from the mean, but their proportion is
##       reported, which makes the asymmetry of the effect explicit rather
##       than hiding it inside a single mean.

# -----------------------------------------------------------------------------
# 7.1 Helpers
# -----------------------------------------------------------------------------

## Rasterise an sf point grid column into a SpatRaster.
sf_points_to_rast <- function(x, column) {
  df <- data.frame(st_coordinates(x), z = st_drop_geometry(x)[[column]])
  names(df) <- c("x", "y", "z")
  terra::rast(df, type = "xyz", crs = st_crs(x)$wkt)
}

## Regular measurement points along a polygon outline.
boundary_points <- function(poly, spacing) {
  g <- st_geometry(poly)
  b <- st_boundary(g)
  b <- st_segmentize(b, dfMaxLength = spacing)
  b <- suppressWarnings(st_cast(b, "LINESTRING"))
  p <- suppressWarnings(st_cast(b, "POINT"))
  st_sf(geometry = p)
}

## Core routine. Returns one row per measurement point.
##
##   change_sf : sf point grid with the posterior summary columns
##   column    : "q0.975" (significance) or "mean" (zero effect)
##   owf       : sf polygon(s) of ONE cluster (will be dissolved)
##   label     : cluster label for the output
displacement_distances <- function(change_sf, column, owf, label,
                                   spacing = OWF_BOUNDARY_SPACING_KM) {

  owf_u <- st_sf(geometry = st_union(st_geometry(owf)))

  ## --- affected region A (criterion < 0) --------------------------------
  r <- sf_points_to_rast(change_sf, column)
  r_neg <- terra::classify(r, matrix(c(-Inf, 0, 1, 0, Inf, NA), ncol = 3,
                                     byrow = TRUE))

  A <- terra::as.polygons(r_neg) %>%
    st_as_sf() %>%
    st_set_crs(st_crs(change_sf)) %>%
    st_make_valid() %>%
    st_cast("POLYGON", warn = FALSE)

  ## keep the connected component(s) touching this cluster
  hit <- lengths(st_intersects(A, owf_u)) > 0
  if (!any(hit)) {
    warning("No affected region adjacent to cluster '", label,
            "' for criterion '", column, "'.")
    return(NULL)
  }
  A_c <- st_sf(geometry = st_union(st_geometry(A[hit, ])))

  ## --- measurement geometry ---------------------------------------------
  ## Use the smooth marching-squares isoline where available (avoids the
  ## staircase artefact of raster polygonisation); fall back to the
  ## polygonised boundary.
  edge <- tryCatch({
    ## maxcells must cover the whole grid, otherwise terra silently
    ## downsamples (default 1e5) and the isoline loses resolution.
    iso <- terra::as.contour(r, maxcells = terra::ncell(r), levels = 0) %>%
      st_as_sf() %>%
      st_set_crs(st_crs(change_sf)) %>%
      st_cast("LINESTRING", warn = FALSE)
    keep <- lengths(st_intersects(
      iso, st_buffer(st_boundary(A_c), 1.5 * max(terra::res(r))))) > 0
    if (any(keep)) st_union(st_geometry(iso[keep, ])) else st_boundary(A_c)
  }, error = function(e) st_boundary(A_c))

  ## --- measure -----------------------------------------------------------
  pts <- boundary_points(owf_u, spacing)
  inside <- lengths(st_intersects(pts, A_c)) > 0
  d <- as.numeric(st_distance(pts, edge))

  tibble::tibble(
    cluster   = label,
    criterion = column,
    point_id  = seq_len(nrow(pts)),
    affected  = inside,
    distance  = ifelse(inside, d, NA_real_),
    x = st_coordinates(pts)[, 1],
    y = st_coordinates(pts)[, 2]
  )
}

## Mean, SD, n and normal-approximation 95% CI on the mean.
##
## WARNING (REVIEW NOTE-13). The ci_lo/ci_hi columns treat the measurement
## points as independent replicates. They are not: consecutive points are
## OWF_BOUNDARY_SPACING_KM apart along a smooth contour, so the effective
## sample size is a small fraction of n_affected and this interval is far too
## narrow. It is reproduced here only because it is the interval reported in
## the submitted manuscript for the significance contour. For the zero-effect
## contour use the posterior interval from Section 8 instead.
summarise_distances <- function(d) {
  d %>%
    group_by(cluster, criterion) %>%
    summarise(
      n_points   = n(),
      n_affected = sum(affected),
      prop_affected = mean(affected),
      mean_km  = mean(distance, na.rm = TRUE),
      sd_km    = sd(distance,   na.rm = TRUE),
      median_km = median(distance, na.rm = TRUE),
      q25_km   = unname(quantile(distance, 0.25, na.rm = TRUE)),
      q75_km   = unname(quantile(distance, 0.75, na.rm = TRUE)),
      max_km   = max(distance, na.rm = TRUE),
      se       = sd_km / sqrt(n_affected),
      ci_lo    = mean_km - 1.96 * se,
      ci_hi    = mean_km + 1.96 * se,
      .groups  = "drop"
    )
}

## Hedges' g with the small-sample correction (no extra dependency).
hedges_g <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  d  <- (mean(x) - mean(y)) / sp
  J  <- 1 - 3 / (4 * (nx + ny) - 9)
  d * J
}

# -----------------------------------------------------------------------------
# 7.2 Build the cluster polygons
# -----------------------------------------------------------------------------

farms_after$cluster <- as.character(farms_after$cluster)

owf_north <- farms_after %>% filter(cluster %in% CLUSTER_DEF$north)

owf_south <- farms_after %>% filter(cluster %in% CLUSTER_DEF$south)
if (nrow(owf_south) == 0) {
  message("Cluster label '", paste(CLUSTER_DEF$south, collapse = "/"),
          "' not found - falling back to explicit farm names.")
  owf_south <- farms_after %>% filter(Name %in% SOUTH_FARM_NAMES)
}

stopifnot(nrow(owf_north) > 0, nrow(owf_south) > 0)
message("north cluster OWFs: ", paste(owf_north$Name, collapse = ", "))
message("south cluster OWFs: ", paste(owf_south$Name, collapse = ", "))

# -----------------------------------------------------------------------------
# 7.3 Compute both distance measures
# -----------------------------------------------------------------------------

dist_all <- bind_rows(
  ## (a) conservative lower bound - as in the submitted Fig. 6
  displacement_distances(change_intensity, "q0.975", owf_north, "North"),
  displacement_distances(change_intensity, "q0.975", owf_south, "South"),
  ## (b) NEW: central estimate at the zero-effect contour
  displacement_distances(change_intensity, "mean",   owf_north, "North"),
  displacement_distances(change_intensity, "mean",   owf_south, "South")
) %>%
  mutate(criterion = factor(
    dplyr::recode(criterion,
                  "q0.975" = "95% significance contour",
                  "mean"   = "Zero-effect contour"),
    levels = c("95% significance contour", "Zero-effect contour")
  ))

dist_summary <- summarise_distances(dist_all)
print(as.data.frame(dist_summary), digits = 3)

write.csv(dist_all,     p_out("displacement_distances_points_2026.csv"), row.names = FALSE)
write.csv(dist_summary, p_out("displacement_distances_summary_2026.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 7.4 North vs South, Welch t-test
# -----------------------------------------------------------------------------
## Restricted to the SIGNIFICANCE contour, purely to reproduce the number
## already in the manuscript (Welch t = 2.21, p = 0.027, Hedge's g = 0.22).
## It inherits the pseudo-replication problem described above and should not
## be extended to the zero-effect contour - Section 8 does that properly.

tests <- lapply("95% significance contour", function(cr) {
  n <- dist_all$distance[dist_all$criterion == cr & dist_all$cluster == "North"]
  s <- dist_all$distance[dist_all$criterion == cr & dist_all$cluster == "South"]
  tt <- t.test(n, s)                      # Welch
  data.frame(
    criterion = cr,
    mean_north = mean(n, na.rm = TRUE),
    mean_south = mean(s, na.rm = TRUE),
    welch_t = unname(tt$statistic),
    df      = unname(tt$parameter),
    p_value = tt$p.value,
    hedges_g = hedges_g(n, s)
  )
}) %>% bind_rows()

print(tests, digits = 3)
write.csv(tests, p_out("displacement_tests_2026.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 7.5 Figure 6 - violin plot, now with both contours side by side
# -----------------------------------------------------------------------------

p_violin <- ggplot(dist_all %>% filter(!is.na(distance)),
                   aes(x = cluster, y = distance, fill = criterion)) +
  geom_violin(position = position_dodge(0.8), alpha = 0.6,
              scale = "width", trim = TRUE) +
  geom_boxplot(position = position_dodge(0.8), width = 0.12,
               outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = c("95% significance contour" = "#B2182B",
                               "Zero-effect contour"      = "#2166AC"),
                    name = NULL) +
  labs(x = NULL, y = "Displacement distance (km)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(p_out("fig6_displacement_distances_2026.png"), p_violin,
       width = 7.5, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 7.6 Map of the measurement geometry (for the response to reviewers)
# -----------------------------------------------------------------------------

p_geom <- ggplot() +
  geom_raster(data = df_intens, aes(x = x, y = y, fill = mean)) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(3, "RdBu"),
                       name = "log10 change") +
  geom_contour(aes(x = x, y = y, z = mean), data = df_intens,
               breaks = 0, colour = "grey20", linewidth = 0.6) +
  geom_contour(aes(x = x, y = y, z = q0.975), data = df_intens,
               breaks = 0, colour = "red", linetype = 2, linewidth = 0.6) +
  geom_sf(data = owf_north, fill = NA, colour = "black") +
  geom_sf(data = owf_south, fill = NA, colour = "black") +
  geom_point(data = dist_all %>% filter(affected),
             aes(x = x, y = y), size = 0.3, colour = "black") +
  coord_sf() +
  ggtitle("Zero-effect (grey) vs 95% significance (red dashed) contours") +
  theme_minimal(base_size = 13)

ggsave(p_out("fig6b_contour_geometry_2026.png"), p_geom,
       width = 8.3, height = 5.5, dpi = 300)


###############################################################################
# =============================================================================
# 8. ZERO-EFFECT DISTANCE WITH FULL POSTERIOR UNCERTAINTY
# =============================================================================
###############################################################################
#
#  REVIEW NOTE-13. Section 7 derives the contour from a single posterior
#  summary field and then measures hundreds of points along it. Those points
#  are 0.5 km apart on a smooth field: they are not independent replicates,
#  so a mean +- 1.96 SE over them, or a Welch t-test between clusters,
#  describes the length of the contour rather than the uncertainty in the
#  displacement. That is what produces an implausibly tight interval such as
#  7.59 km [7.09, 8.09] and a marginal p with a small effect size.
#
#  Here the whole measurement is repeated INSIDE each posterior draw:
#
#     for j in 1..N_SAMPLES_DIST
#         delta_j(s)  <- one posterior sample of log10 exp(spde_change(s))
#         A_j         <- connected negative region of delta_j touching the cluster
#         d_j         <- mean distance from the OWF outline to the zero contour of A_j
#
#  The posterior of d is then summarised directly. The North-vs-South
#  comparison is the PAIRED difference d_north,j - d_south,j within each draw
#  (both clusters share the same latent field, so the draws must be paired),
#  reported as P(d_north > d_south) instead of a p-value.
#
#  This applies to the ZERO-EFFECT contour only. The 95% significance contour
#  is a property of the posterior, not of a single draw, so it has no
#  per-sample analogue; the Section 7 numbers for it stand as published.

# -----------------------------------------------------------------------------
# 8.1 Posterior samples of the change field
# -----------------------------------------------------------------------------

## Deliberately NOT cached: this matrix is n_cells x n_draws (easily >100 MB)
## and this folder is synced to Mega. It is seeded, so re-running reproduces
## it exactly; the expensive part (the model fit) is what gets cached.
set.seed(SEED + 3)
change_samples <- generate(change_fit, pxl,
                           ~ log10(exp(spde_change)),
                           n.samples = N_SAMPLES_DIST)

## generate() returns a matrix: rows = prediction cells, columns = samples.
change_samples <- as.matrix(change_samples)
stopifnot(nrow(change_samples) == nrow(pxl))
message("posterior sample matrix: ", nrow(change_samples), " cells x ",
        ncol(change_samples), " draws")

## Raster template covering the prediction grid, plus the mapping from pxl
## rows to raster cells. Built once; only the values change per draw.
grid_template <- sf_points_to_rast(change_intensity, "mean")
pxl_xy        <- st_coordinates(pxl)

# -----------------------------------------------------------------------------
# 8.2 Distance for one realisation of the field
# -----------------------------------------------------------------------------

## Returns the mean distance from the OWF outline to the zero contour of the
## connected negative region adjacent to that OWF, for a single field.
##   r     : SpatRaster of one draw, already cropped to the cluster window
##   owf_u : dissolved cluster polygon
##   pts    : pre-computed measurement points on the OWF outline
zero_distance_one <- function(r, owf_u, pts) {

  none <- list(mean_km = NA_real_, prop_affected = 0,
               prop_censored = NA_real_, has_zone = FALSE)

  vals <- terra::values(r)
  if (all(is.na(vals)) || !any(vals < 0, na.rm = TRUE)) return(none)

  r_neg <- terra::classify(r, matrix(c(-Inf, 0, 1, 0, Inf, NA),
                                     ncol = 3, byrow = TRUE))

  A <- try(suppressWarnings(
    st_cast(
      st_make_valid(
        st_set_crs(st_as_sf(terra::as.polygons(r_neg)), st_crs(pts))),
      "POLYGON")), silent = TRUE)
  if (inherits(A, "try-error") || nrow(A) == 0) return(none)

  hit <- lengths(st_intersects(A, owf_u)) > 0
  if (!any(hit)) return(none)
  A_c <- st_sf(geometry = st_union(st_geometry(A[hit, ])))

  ## Smooth isoline where possible, polygonised boundary as fallback.
  edge <- tryCatch({
    iso <- suppressWarnings(st_cast(
      st_set_crs(st_as_sf(terra::as.contour(r, maxcells = terra::ncell(r),
                                            levels = 0)), st_crs(pts)),
      "LINESTRING"))
    k <- lengths(st_intersects(
      iso, st_buffer(st_boundary(A_c), 1.5 * max(terra::res(r))))) > 0
    if (any(k)) st_union(st_geometry(iso[k, ])) else st_geometry(st_boundary(A_c))
  }, error = function(e) st_geometry(st_boundary(A_c)))

  inside <- lengths(st_intersects(pts, A_c)) > 0
  if (!any(inside)) return(none)

  d <- as.numeric(st_distance(pts, edge))
  d[!inside] <- NA_real_

  ## Right-censoring: the affected zone reaches the edge of the search window,
  ## so the true zero crossing is somewhere beyond MAX_SEARCH_KM.
  censored <- d >= 0.98 * MAX_SEARCH_KM

  list(mean_km       = mean(d, na.rm = TRUE),
       prop_affected = mean(inside),
       prop_censored = mean(censored, na.rm = TRUE),
       has_zone      = TRUE)
}

# -----------------------------------------------------------------------------
# 8.3 Loop over posterior draws, per cluster
# -----------------------------------------------------------------------------

posterior_displacement <- function(samples, cluster_owf, label) {

  owf_u <- st_sf(geometry = st_union(st_geometry(cluster_owf)))
  pts   <- boundary_points(owf_u, OWF_BOUNDARY_SPACING_KM)
  win   <- st_buffer(owf_u, MAX_SEARCH_KM)

  ## Crop once: the search window is much smaller than the full domain, which
  ## makes the per-draw polygonisation cheap AND removes any risk of measuring
  ## against an unrelated effect zone elsewhere in the German Bight.
  tmpl <- terra::crop(grid_template, terra::vect(st_geometry(win)))
  cid  <- terra::cellFromXY(tmpl, pxl_xy)   # NA for cells outside the window
  keep <- !is.na(cid)
  nc   <- terra::ncell(tmpl)

  message("Cluster ", label, ": ", nrow(pts), " measurement points, ",
          sum(keep), " cells in window, ", ncol(samples), " draws")

  res <- lapply(seq_len(ncol(samples)), function(j) {
    v <- rep(NA_real_, nc)
    v[cid[keep]] <- samples[keep, j]
    zero_distance_one(terra::setValues(tmpl, v), owf_u, pts)
  })

  tibble::tibble(
    cluster       = label,
    draw          = seq_along(res),
    mean_km       = vapply(res, function(z) z$mean_km,       numeric(1)),
    prop_affected = vapply(res, function(z) z$prop_affected, numeric(1)),
    prop_censored = vapply(res, function(z) z$prop_censored, numeric(1)),
    has_zone      = vapply(res, function(z) z$has_zone,      logical(1))
  )
}

post_north <- cached(p_out("post_dist_north_2026.rds"),
  posterior_displacement(change_samples, owf_north, "North"))
post_south <- cached(p_out("post_dist_south_2026.rds"),
  posterior_displacement(change_samples, owf_south, "South"))

post_dist <- bind_rows(post_north, post_south)
write.csv(post_dist, p_out("zero_contour_posterior_draws_2026.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 8.4 Posterior summary
# -----------------------------------------------------------------------------

post_summary <- post_dist %>%
  group_by(cluster) %>%
  summarise(
    n_draws         = n(),
    ## P(no negative zone adjacent to this cluster in a given draw). A large
    ## value means the displacement estimate is not well supported.
    p_no_zone       = mean(!has_zone),
    post_mean_km    = mean(mean_km,   na.rm = TRUE),
    post_median_km  = median(mean_km, na.rm = TRUE),
    ci_lo_km        = unname(quantile(mean_km, 0.025, na.rm = TRUE)),
    ci_hi_km        = unname(quantile(mean_km, 0.975, na.rm = TRUE)),
    ## Fraction of the OWF outline with a negative zone, and fraction of the
    ## measured points where the zero crossing lies beyond the search window.
    mean_prop_affected = mean(prop_affected, na.rm = TRUE),
    mean_prop_censored = mean(prop_censored, na.rm = TRUE),
    .groups = "drop"
  )

print(as.data.frame(post_summary), digits = 3)
write.csv(post_summary, p_out("zero_contour_posterior_summary_2026.csv"),
          row.names = FALSE)

if (any(post_summary$mean_prop_censored > 0.05, na.rm = TRUE)) {
  warning("Some draws are right-censored at MAX_SEARCH_KM = ", MAX_SEARCH_KM,
          " km. The zero contour extends beyond the search window; increase ",
          "MAX_SEARCH_KM or report the estimate as a lower bound.")
}

# -----------------------------------------------------------------------------
# 8.5 North vs South: paired posterior difference
# -----------------------------------------------------------------------------
## Both clusters are measured on the SAME draw, so the comparison is paired.

stopifnot(identical(post_north$draw, post_south$draw))
delta <- post_north$mean_km - post_south$mean_km

post_contrast <- data.frame(
  comparison        = "North - South (zero-effect contour)",
  n_draws_usable    = sum(!is.na(delta)),
  post_mean_diff_km = mean(delta,   na.rm = TRUE),
  post_median_diff  = median(delta, na.rm = TRUE),
  ci_lo_km          = unname(quantile(delta, 0.025, na.rm = TRUE)),
  ci_hi_km          = unname(quantile(delta, 0.975, na.rm = TRUE)),
  ## Posterior probability that displacement is larger in the north.
  ## This replaces the Welch p-value; there is no multiple-testing or
  ## pseudo-replication issue because each draw is one coherent realisation.
  p_north_gt_south  = mean(delta > 0, na.rm = TRUE)
)

print(post_contrast, digits = 3)
write.csv(post_contrast, p_out("zero_contour_posterior_contrast_2026.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 8.6 Figure: posterior of the zero-effect displacement
# -----------------------------------------------------------------------------

sig_means <- dist_summary %>%
  filter(criterion == "95% significance contour") %>%
  select(cluster, mean_km)

p_post <- ggplot(post_dist %>% filter(!is.na(mean_km)),
                 aes(x = mean_km, fill = cluster)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = sig_means, aes(xintercept = mean_km, colour = cluster),
             linetype = 2, linewidth = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c(North = "#2166AC", South = "#B2182B"),
                    name = NULL) +
  scale_colour_manual(values = c(North = "#2166AC", South = "#B2182B")) +
  labs(x = "Mean displacement distance to the zero-effect contour (km)",
       y = "Posterior density",
       caption = paste("Dashed lines: 95% significance-contour estimates",
                       "(conservative lower bound, Section 7)")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(p_out("fig7_zero_contour_posterior_2026.png"), p_post,
       width = 7.5, height = 5, dpi = 300)


# =============================================================================
# 9. (Optional) Number of individuals displaced within each affected zone
# =============================================================================
## Same machinery: integrate the before-period and after-period posterior mean
## density over the affected region of each cluster and take the difference.

if (TRUE) {

  affected_region <- function(change_sf, column, owf) {
    r <- sf_points_to_rast(change_sf, column)
    r_neg <- terra::classify(r, matrix(c(-Inf, 0, 1, 0, Inf, NA),
                                       ncol = 3, byrow = TRUE))
    A <- terra::as.polygons(r_neg) %>% st_as_sf() %>%
      st_set_crs(st_crs(change_sf)) %>% st_make_valid() %>%
      st_cast("POLYGON", warn = FALSE)
    owf_u <- st_union(st_geometry(owf))
    hit <- lengths(st_intersects(A, owf_u)) > 0
    if (!any(hit)) return(NULL)
    st_sf(geometry = st_union(st_geometry(A[hit, ])))
  }

  abundance_change <- function(region, label, criterion) {
    if (is.null(region)) return(NULL)
    idx <- lengths(st_intersects(pxl, region)) > 0
    data.frame(
      cluster   = label,
      criterion = criterion,
      area_km2  = sum(idx) * cell_area,
      n_before  = sum(dens_before$mean[idx]) * cell_area,
      n_after   = sum(dens_after$mean[idx])  * cell_area
    ) %>% mutate(n_displaced = n_before - n_after)
  }

  displaced <- bind_rows(
    abundance_change(affected_region(change_intensity, "q0.975", owf_north),
                     "North", "95% significance contour"),
    abundance_change(affected_region(change_intensity, "q0.975", owf_south),
                     "South", "95% significance contour"),
    abundance_change(affected_region(change_intensity, "mean", owf_north),
                     "North", "Zero-effect contour"),
    abundance_change(affected_region(change_intensity, "mean", owf_south),
                     "South", "Zero-effect contour")
  )

  print(displaced, digits = 4)
  write.csv(displaced, p_out("displaced_individuals_2026.csv"), row.names = FALSE)
}


# =============================================================================
# 10. Grid-resolution sensitivity
# =============================================================================
## The measured distance depends on the prediction grid: too coarse and the
## contour is quantised, too fine and nothing changes but the run time. This
## needs NO refit - only a re-prediction on a different grid - so it is cheap
## and worth reporting. Set to FALSE if you only want the main results.

RUN_GRID_SENSITIVITY <- TRUE

if (RUN_GRID_SENSITIVITY) {

  grid_check <- lapply(c(0.5, 2.0), function(res_km) {
    xr <- range(spatial_mesh$loc[, 1]); yr <- range(spatial_mesh$loc[, 2])
    g  <- fmesher::fm_pixels(
      spatial_mesh,
      dims = c(ceiling(diff(xr) / res_km), ceiling(diff(yr) / res_km)),
      mask = eez, format = "sf")
    st_crs(g) <- main_crs

    set.seed(SEED + 2)   # same seed as the main prediction
    ci <- predict(change_fit, g, ~ log10(exp(spde_change)),
                  n.samples = N_SAMPLES_CHANGE)

    bind_rows(
      displacement_distances(ci, "mean", owf_north, "North"),
      displacement_distances(ci, "mean", owf_south, "South")
    ) %>%
      group_by(cluster) %>%
      summarise(res_km = res_km,
                mean_km = mean(distance, na.rm = TRUE),
                prop_affected = mean(affected),
                .groups = "drop")
  }) %>% bind_rows()

  ## The main run, for comparison.
  grid_check <- bind_rows(
    grid_check,
    dist_summary %>%
      filter(criterion == "Zero-effect contour") %>%
      transmute(cluster, res_km = PRED_CELL_KM, mean_km, prop_affected)
  ) %>% arrange(cluster, res_km)

  print(as.data.frame(grid_check), digits = 3)
  write.csv(grid_check, p_out("grid_resolution_sensitivity_2026.csv"),
            row.names = FALSE)
}


# =============================================================================
# 11. Sensitivity checks that DO require a refit
# =============================================================================
## Not run by default - each is a full model fit (hours). Both are things
## Reviewer #3 asked for, so they are scaffolded here rather than described.
##
##  (a) Prior sensitivity. Set PRIOR_RANGE_CHANGE to e.g. c(5, NA) and
##      c(30, NA) at the top of the file, set REFRESH_FITS <- TRUE, and
##      re-run Sections 6-8 writing to a different cache prefix. The
##      question is whether the zero-contour distance tracks the prior
##      range; if it does, the range prior is doing the work and that has
##      to be said in the Discussion.
##
##  (b) Baseline robustness. Refit the "before" period dropping 2003 (the
##      strong aggregation reported in Results):
##        before_alt <- dat %>% filter(phase <= BEFORE_PHASE_MAX, phase != 3)
##      then re-run the joint model with obs_before built from before_alt.
##      Check that the displacement estimates are not driven by that one
##      season.


# =============================================================================
# 12. Session info (for the reproducibility statement)
# =============================================================================
writeLines(capture.output(sessionInfo()), p_out("sessionInfo_2026.txt"))


###############################################################################
#  REVIEW NOTES - issues found in the Dec-2022 script
#  (see the accompanying code_review.md for the full discussion)
#
#  NOTE-1  Three consecutive setwd() calls to three different machines, and
#          two devtools::load_all() of a vendored inlabru. Removed.
#  NOTE-2  p_data("prediction_pixels.rds") and "../modelos diverpop/count_model_3.R"
#          are referenced but are NOT in this folder. Grid is now rebuilt.
#  NOTE-3  Magic number 1.136364 applied as an exposure correction, with no
#          comment; conflicts with the Methods text (mrds staggered by group
#          size and sea state).
#  NOTE-4  Figs. 3/4 were built from a single rnbinom draw per cell and then
#          scale()d SEPARATELY per period, while sharing a legend labelled
#          [Ind. km^-2]. The two panels were therefore not comparable and the
#          units were wrong. Fixed.
#  NOTE-5  is_daisi + is_apem + is_conventional + Intercept is rank-deficient
#          (dummy-variable trap). Reference level now dropped.
#  NOTE-6  ~40 lines of inlabru-vignette boilerplate (df1/df2, lik1/lik2 with
#          undefined x, y, z) sat between the two models and would have
#          errored. Deleted.
#  NOTE-7  ips = pxl.all passed to a family = "nbinomial" observation model.
#          Only meaningful for family = "cp". Dropped.
#  NOTE-8  data$NHAT <- round(NHAT) discards part of the mrds detection
#          correction. Kept (required by nbinomial) but flagged.
#  NOTE-9  proj4string(farms) <- proj4string(before) ASSIGNED a CRS instead of
#          reprojecting. Now handled explicitly in to_sf_km().
#  NOTE-10 wf_distance()/gDistance computed distances to ALL farms including
#          those not yet built in the before period, and the result was never
#          used. Removed.
#  NOTE-11 n.samples = 100 for the density maps is too few for stable 2.5/97.5
#          quantiles. Raised to N_SAMPLES_CHANGE (500).
#  NOTE-12 The joint change model contains NO method covariate, although the
#          exploratory model does. This is the confounding the reviewers
#          raised; it is a modelling choice, not a bug, but code and Methods
#          text must agree.
#  NOTE-13 PSEUDO-REPLICATION. Measurement points 0.5 km apart along a smooth
#          contour were treated as independent replicates, both for the CI on
#          the mean displacement and for the Welch t-test between clusters.
#          The effective sample size is a small fraction of the nominal n, so
#          the published interval (+-0.5 km on a 7.6 km estimate) is far too
#          narrow and the cluster comparison (p = 0.027, g = 0.22) is the
#          signature of an inflated n. Section 8 replaces this with a
#          per-draw posterior for the zero-effect contour. The significance
#          contour has no per-draw analogue and is left as published.
#  NOTE-14 No seed was set before any posterior sampling call, so none of the
#          published numbers were exactly reproducible. Fixed (SEED).
###############################################################################
