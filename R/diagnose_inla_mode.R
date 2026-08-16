###############################################################################
#  Does inla.mode explain the shift from the 2022 fit?
#
#  INLA introduced "compact" as an internal reparameterisation and later made
#  it the default, replacing "classic". The 2022 analysis will have run under
#  "classic"; the current one runs under "compact" (inla.getOption confirms
#  it). The two are not "more" or "less" approximate - they are different
#  formulations of the same problem - but they give slightly different
#  numbers.
#
#  Reference points established so far, all with terra::distance on the
#  identical measurement code:
#
#      2022 field                 north 7.465   south 6.388
#      current field ("compact")  north 7.15    south 6.20
#
#  If refitting under "classic" lands near 7.465 / 6.388, inla.mode accounts
#  for the whole refit difference.
#
#  Costs one full refit of the joint model. Nothing already computed is
#  overwritten: the classic fit is cached separately.
#
#      source("R/diagnose_inla_mode.R")
###############################################################################

library(INLA); library(inlabru); library(fmesher)
library(sf); library(terra); library(dplyr)

DATA_DIR <- "data"; OUT_DIR <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)

sf_use_s2(FALSE)
main_crs <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")

BEFORE_PHASE_MAX <- 8
AFTER_PHASE_MIN  <- 17
VISUAL_DETECTION_CORRECTION <- 1.136364
PRIOR_SIGMA <- c(0.2, 0.01)
PRIOR_RANGE <- c(15, NA)
N_SAMPLES_DIST <- 1000
SEED <- 20260811L
DIST_CELL_KM <- 1
THR_SIGNIF <- 0.975
CLUSTER_NORTH <- c("West of Sylt", "Butendiek")
CLUSTER_SOUTH <- c("UMBO")

as_km <- function(x) {
  y <- if (inherits(x, "sf")) x else st_as_sf(x); st_crs(y) <- main_crs; y
}
D <- new.env(parent = emptyenv())
load(p_data("diver_owf_data.RData"), envir = D)

countdata    <- as_km(D$counts_5km)
farms        <- as_km(D$owf_polygons)
pxl          <- as_km(D$prediction_pxl)
pred_outline <- as_km(D$prediction_mask)
spatial_mesh <- fm_as_mesh_2d(D$mesh_5km)
fm_crs(spatial_mesh) <- main_crs

if (!all(as.character(st_geometry_type(pxl)) == "POINT")) {
  pxl <- suppressWarnings(st_centroid(pxl))
}
pxl <- suppressWarnings(st_filter(pxl, st_union(st_geometry(pred_outline))))
pxl_xy <- st_coordinates(pxl)

dat <- countdata %>%
  filter(!is.na(area), !is.na(NHAT), area > 0) %>%
  mutate(NHAT = round(NHAT))
before <- dat %>% filter(phase <= BEFORE_PHASE_MAX)
after  <- dat %>% filter(phase >= AFTER_PHASE_MIN)

owf <- farms %>% group_by(Name, cluster) %>% summarise(.groups = "drop")
owf_u <- list(
  north = st_sf(geometry = st_union(st_geometry(
    owf %>% filter(cluster %in% CLUSTER_NORTH)))),
  south = st_sf(geometry = st_union(st_geometry(
    owf %>% filter(cluster %in% CLUSTER_SOUTH))))
)

# ---------------------------------------------------------------- the refit
matern <- inla.spde2.pcmatern(spatial_mesh,
                              prior.sigma = PRIOR_SIGMA,
                              prior.range = PRIOR_RANGE)

cmp_joint <- ~
  spde_before(geometry, model = matern) +
  spde_change(geometry, model = matern) +
  Intercept(1) +
  Intercept_after(1)

inla_opts <- list(
  control.compute   = list(dic = TRUE, cpo = TRUE, waic = TRUE),
  control.inla      = list(int.strategy = "ccd"),
  control.predictor = list(compute = TRUE,
                           quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975))
)

CLASSIC_FIT <- p_out("change_fit_classic.rds")

if (file.exists(CLASSIC_FIT)) {
  cat("loading cached classic fit\n")
  change_fit_classic <- readRDS(CLASSIC_FIT)
} else {
  ## Set globally rather than through options(), so it applies regardless of
  ## whether inlabru plumbs the argument through in this version.
  old_mode <- INLA::inla.getOption("inla.mode")
  cat("inla.mode was:", old_mode, "-> switching to classic\n")
  INLA::inla.setOption(inla.mode = "classic")
  on.exit(INLA::inla.setOption(inla.mode = old_mode), add = TRUE)
  cat("inla.mode now:", INLA::inla.getOption("inla.mode"), "\n")

  change_fit_classic <- bru(
    cmp_joint,
    bru_obs("nbinomial", data = before,
            formula = NHAT ~ spde_before + Intercept,
            E = before$area / VISUAL_DETECTION_CORRECTION),
    bru_obs("nbinomial", data = after,
            formula = NHAT ~ spde_before + spde_change + Intercept_after,
            E = after$area),
    options = c(inla_opts, list(bru_max_iter = 1)))

  saveRDS(change_fit_classic, CLASSIC_FIT)
  INLA::inla.setOption(inla.mode = old_mode)
  cat("inla.mode restored to:", INLA::inla.getOption("inla.mode"), "\n")
}

# ---------------------------------------------------------------- measure
samples <- as.matrix(generate(change_fit_classic, pxl,
                              ~ log10(exp(spde_change)),
                              n.samples = N_SAMPLES_DIST, seed = SEED + 3L))
p_decrease <- rowMeans(samples < 0)

grid_template <- terra::rast(
  data.frame(x = pxl_xy[, 1], y = pxl_xy[, 2], z = p_decrease),
  type = "xyz", crs = main_crs$wkt)
cell_idx <- terra::cellFromXY(grid_template, pxl_xy)

zone_raster <- function(thr) {
  v <- rep(NA_real_, terra::ncell(grid_template))
  v[cell_idx] <- thr - p_decrease
  terra::setValues(grid_template, v)
}

affected_region <- function(r, owf_one) {
  r_neg <- terra::classify(r, matrix(c(-Inf, 0, 1, 0, Inf, NA),
                                     ncol = 3, byrow = TRUE))
  if (all(is.na(terra::values(r_neg)))) return(NULL)
  A <- suppressWarnings(st_cast(st_make_valid(
    st_set_crs(st_as_sf(terra::as.polygons(r_neg)), main_crs)), "POLYGON"))
  hit <- lengths(st_intersects(A, owf_one)) > 0
  if (!any(hit)) return(NULL)
  st_sf(geometry = st_union(st_geometry(A[hit, ])))
}

make_dist_raster <- function(owf_one) {
  d <- terra::rast(terra::ext(grid_template), resolution = DIST_CELL_KM,
                   crs = main_crs$wkt, vals = 0)
  terra::distance(d, terra::vect(st_geometry(owf_one))) / 1000
}

res <- bind_rows(lapply(names(owf_u), function(cl) {
  A_c <- affected_region(zone_raster(THR_SIGNIF), owf_u[[cl]])
  if (is.null(A_c)) return(NULL)
  d <- as.numeric(stats::na.omit(terra::extract(
    make_dist_raster(owf_u[[cl]]),
    terra::vect(st_geometry(st_boundary(A_c))))[[2]]))
  data.frame(cluster = cl, n = length(d),
             area_km2 = as.numeric(sum(st_area(A_c))),
             mean_km = mean(d), median_km = median(d))
}))

cat("\n=== inla.mode = \"classic\" ===\n")
print(as.data.frame(res), digits = 5, row.names = FALSE)

cat("\n=== reference points (all terra::distance, same measurement code) ===\n")
cat("2022 field                  north  N=256  area=1170.6  mean=7.465\n")
cat("current fit, compact        north  N=238  area=1138    mean=7.146\n")
cat("2022 field                  south  N=176  area= 879.1  mean=6.388\n")
cat("current fit, compact        south  N=164  area= 860    mean=6.201\n")
cat("\npublished (distanceto): north 7.59, south 6.62\n")
cat("\nIf the classic refit sits near 7.465 / 6.388, inla.mode accounts for\n")
cat("the refit difference. If it sits near the compact values, it does not.\n\n")

write.csv(res, p_out("inla_mode_comparison.csv"), row.names = FALSE)
