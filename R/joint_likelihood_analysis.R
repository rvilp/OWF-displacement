###############################################################################
#  Diver (Gavia spp.) displacement around offshore wind farms, German Bight
#  Joint-likelihood before/after LGCP model
#
#  JEMA-D-25-08275
#
#  Before period 2001-2008 (phase <= 8), after period 2017-2021 (phase >= 17).
#
#  This is the clean working version. The heavily annotated
#  joint_likelihood_before_after.R documents the port from the 2022 code and
#  the reasoning behind each change; keep it as the reference.
#
#  Reported results were produced with R 4.5.2, INLA 25.10.19, inlabru 2.13.0,
#  fmesher 0.7.0, sf 1.1-0, terra 1.8-93. The full session is written to
#  outputs/sessionInfo.txt at the end of the run.
#
#  Model fits are cached in outputs/ because they take hours. Predictions are
#  always recomputed - they take minutes and depend on the prediction grid.
#  The cache file names carry MODEL_TAG, so a change to the model formula
#  cannot pick up an old fit; to force a refit anyway, delete
#  outputs/fit_*_<tag>.rds and outputs/change_fit_<tag>.rds.
###############################################################################

# =============================================================================
# 0. Configuration
# =============================================================================

DATA_DIR <- "data"
OUT_DIR  <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)
dir.create(OUT_DIR, showWarnings = FALSE)

library(INLA)
library(inlabru)
library(fmesher)
library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)

sf_use_s2(FALSE)

## All inputs are stored in UTM zone 32N with KILOMETRE units. Coordinates are
## already in that system, so the CRS is assigned, not transformed.
main_crs <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")

BEFORE_PHASE_MAX <- 8
AFTER_PHASE_MIN  <- 17

## ---------------------------------------------------------------------------
## Detection correction by survey technique
##
## Every observation's exposure is effort x the detection rate of the
## technique that produced it, so that all techniques are placed on a common
## scale before the periods are compared. HiDef is the reference: it is the
## technique for which a 100% detection rate has been established.
##
## The values are the medians of Table A-1 of Vilela et al. (2021), which
## reports the technique effects on the LOG scale relative to HiDef:
##
##       technique   log     SD     Q0.025   Q0.975     exp(log)
##       DAISI      -0.10   0.08    -0.26     0.06        0.905
##       APEM       -0.40   0.21    -0.82     0.01        0.670
##       VISUAL     -0.22   0.05    -0.31    -0.13        0.803
##
## The visual figure is the residual technique effect that remains AFTER the
## distance-sampling correction, which is already applied in NHAT.
##
## These enter as a KNOWN OFFSET rather than as a fitted coefficient. A
## fitted technique term is not identifiable in this design: the after period
## aggregates five years with no temporal term, and the techniques are
## segregated in both time and space. APEM and the visual surveys were flown
## only in 2017-2018; DAISI covers the north (y >= 5996) and the visual
## surveys the south (y <= 5983), so the two never sample the same mesh node;
## only 14 of APEM's 130 observations fall on a node-year that HiDef also
## covered. Fitting the term returns technique confounded with year and
## region - it reverses the sign of DAISI and APEM relative to the values
## above. Vilela et al. (2021) estimated these effects in a year-by-year
## spatiotemporal model, where technique is separated from when and where it
## was flown, so the estimate is taken from there.
DETECTION <- c(HiDef        = 1,
               DAISI        = exp(-0.10),
               APEM         = exp(-0.40),
               conventional = exp(-0.22))
## ---------------------------------------------------------------------------

## SPDE PC priors
PRIOR_SIGMA <- c(0.2, 0.01)   # P(sigma > 0.2) = 0.01
PRIOR_RANGE <- c(15, NA)      # median range 15 km

## Tag written into the cached fit names. Bump it whenever the model changes,
## so that an old cache can never be picked up silently.
##   m2  fitted technique indicators in the after period
##   m3  technique as a known detection offset in both periods
MODEL_TAG <- "m3"

N_SAMPLES_MAP  <- 500   # posterior samples for the summary maps
N_SAMPLES_DIST <- 10000  # posterior samples for the distance posterior
SEED <- 20260811L  # integer: passed to INLA's sampler, not to set.seed()

## OWF clusters, matched against farms$cluster. The southern cluster is
## labelled "UMBO" in the shapefile; the manuscript calls it "North of
## Borkum".
CLUSTER_NORTH <- c("West of Sylt", "Butendiek")
CLUSTER_SOUTH <- c("UMBO")

## Resolution of the distance-to-OWF raster used to measure the effect
## distance (impact_distance_and_habitat_loss.R used distance_raster with
## cellsize = 1).
DIST_CELL_KM <- 1

MAP_MARGIN_KM <- 12


# =============================================================================
# 1. Load data
# =============================================================================

## Convert to sf and put everything on the same CRS.
as_km <- function(x) {
  y <- if (inherits(x, "sf")) x else st_as_sf(x)
  st_crs(y) <- main_crs
  y
}

## All inputs live in one bundle. Objects, in load order:
##   counts_5km      detection-corrected counts on the 5 km mesh nodes
##   owf_polygons    OWF footprints, repeated per year
##   prediction_pxl  prediction grid (the one used for the published figures)
##   prediction_mask outline of the prediction area
##   hd_mask         BMU (2009) diver main concentration area
##   spa_mask        SPA Eastern German Bight
##   mesh_5km        SPDE mesh used for fitting
D <- new.env(parent = emptyenv())
load(p_data("diver_owf_data.RData"), envir = D)

countdata    <- as_km(D$counts_5km)
farms        <- as_km(D$owf_polygons)
pxl          <- as_km(D$prediction_pxl)
pred_outline <- as_km(D$prediction_mask)
conc_area    <- as_km(D$hd_mask)
spa_area     <- as_km(D$spa_mask)

spatial_mesh <- fm_as_mesh_2d(D$mesh_5km)
fm_crs(spatial_mesh) <- main_crs

## Prediction grid: cell centres, clipped to the prediction mask.
if (!all(as.character(st_geometry_type(pxl)) == "POINT")) {
  pxl <- suppressWarnings(st_centroid(pxl))
}
pxl <- suppressWarnings(st_filter(pxl, st_union(st_geometry(pred_outline))))

pxl_xy    <- st_coordinates(pxl)
cell_dx   <- min(diff(sort(unique(round(pxl_xy[, 1], 6)))))
cell_dy   <- min(diff(sort(unique(round(pxl_xy[, 2], 6)))))
cell_area <- cell_dx * cell_dy

message("prediction grid: ", nrow(pxl), " cells of ",
        round(cell_dx, 3), " x ", round(cell_dy, 3), " km")


# =============================================================================
# 2. Prepare the modelling data
# =============================================================================

dat <- countdata %>%
  filter(!is.na(area), !is.na(NHAT), area > 0) %>%
  mutate(NHAT = round(NHAT))

## Detection rate of the technique that produced each observation. Anything
## not listed in DETECTION is an error rather than a silent NA in E, which
## would propagate into the fit unnoticed.
dat$method_chr <- as.character(dat$method)
unknown <- setdiff(unique(dat$method_chr), names(DETECTION))
if (length(unknown)) {
  stop("no detection rate for technique(s): ", paste(unknown, collapse = ", "),
       "\nAdd them to DETECTION in the configuration block.", call. = FALSE)
}
dat$det <- unname(DETECTION[dat$method_chr])

## Effective effort: area actually searched, scaled by how much of what was
## there the technique would have seen.
dat$E_eff <- dat$area * dat$det

before <- dat %>% filter(phase <= BEFORE_PHASE_MAX)
after  <- dat %>% filter(phase >= AFTER_PHASE_MIN)

message("before: ", nrow(before), " rows | after: ", nrow(after), " rows")
message("methods before: ", paste(sort(unique(before$method_chr)), collapse = ", "))
message("methods after : ", paste(sort(unique(after$method_chr)),  collapse = ", "))

## Effort by technique and period, at the detection rate applied to each.
tech_applied <- dat %>%
  st_drop_geometry() %>%
  mutate(period = ifelse(phase <= BEFORE_PHASE_MAX, "before",
                  ifelse(phase >= AFTER_PHASE_MIN,  "after", "excluded"))) %>%
  filter(period != "excluded") %>%
  group_by(period, technique = method_chr, detection = det) %>%
  summarise(n = n(), effort_km2 = sum(area), .groups = "drop") %>%
  arrange(desc(period), technique)
print(as.data.frame(tech_applied), digits = 4, row.names = FALSE)
write.csv(tech_applied, p_out("detection_offsets_applied.csv"), row.names = FALSE)

## OWF footprints: dissolve the by-year repetitions into one polygon per farm.
owf <- farms %>%
  group_by(Name, cluster) %>%
  summarise(.groups = "drop")

message("OWF clusters available: ",
        paste(sort(unique(as.character(owf$cluster))), collapse = ", "))

owf_north <- owf %>% filter(cluster %in% CLUSTER_NORTH)
owf_south <- owf %>% filter(cluster %in% CLUSTER_SOUTH)


# =============================================================================
# 3. Models
# =============================================================================

matern <- inla.spde2.pcmatern(spatial_mesh,
                              prior.sigma = PRIOR_SIGMA,
                              prior.range = PRIOR_RANGE)

inla_opts <- list(
  control.compute   = list(dic = TRUE, cpo = TRUE, waic = TRUE),
  control.inla      = list(int.strategy = "ccd"),
  control.predictor = list(compute = TRUE,
                           quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975))
)

## Survey technique does not appear in any linear predictor: it is carried by
## the exposure E, as a known detection offset (see DETECTION above).

## Fits are cached; delete the file to refit.
fit_cached <- function(file, expr) {
  if (file.exists(file)) {
    message("loading fit: ", file)
    return(readRDS(file))
  }
  value <- expr
  saveRDS(value, file)
  value
}

# --- 3a. Separate period fits, for the density maps --------------------------

cmp_before <- ~ space_spde(geometry, model = matern) + Intercept(1)
cmp_after  <- ~ space_spde(geometry, model = matern) + Intercept(1)

fit_b <- fit_cached(p_out(paste0("fit_before_", MODEL_TAG, ".rds")), bru(
  cmp_before,
  bru_obs("nbinomial", formula = NHAT ~ space_spde + Intercept,
          data = before, E = before$E_eff),
  options = inla_opts))

fit_a <- fit_cached(p_out(paste0("fit_after_", MODEL_TAG, ".rds")), bru(
  cmp_after,
  bru_obs("nbinomial", formula = NHAT ~ space_spde + Intercept,
          data = after, E = after$E_eff),
  options = inla_opts))

# --- 3b. Joint-likelihood change model ---------------------------------------
## Both periods share one latent field; spde_change is the log-scale contrast.

cmp_joint <- ~
  spde_before(geometry, model = matern) +
  spde_change(geometry, model = matern) +
  Intercept(1) +
  Intercept_after(1)

change_fit <- fit_cached(p_out(paste0("change_fit_", MODEL_TAG, ".rds")), bru(
  cmp_joint,
  bru_obs("nbinomial",
          data = before,
          formula = NHAT ~ spde_before + Intercept,
          E = before$E_eff),
  bru_obs("nbinomial",
          data = after,
          formula = NHAT ~ spde_before + spde_change + Intercept_after,
          E = after$E_eff),
  options = c(inla_opts, list(bru_max_iter = 1))))

summary(change_fit)


# =============================================================================
# 4. Predictions
# =============================================================================

## REPRODUCIBILITY. set.seed() does NOT control this: the sampling happens
## inside INLA::inla.posterior.sample, which has its own generator. Both
## predict() and generate() take a `seed` argument that is passed through to
## it. With seed != 0 inlabru also forces num.threads = "1:1:1", which makes
## the draws deterministic at the cost of running single-threaded.
##
## Without this, repeated runs of this script differ: two runs gave a mean
## effect distance of 7.136 and 7.090 km for the northern cluster.

## Technique is carried entirely by the exposure, so nothing has to be held
## fixed here: both periods are already on the HiDef scale, i.e. these are
## densities as a technique with 100% detection would have recorded them.
dens_before <- predict(fit_b, pxl, ~ exp(space_spde + Intercept),
                       n.samples = N_SAMPLES_MAP, seed = SEED)
dens_after  <- predict(fit_a, pxl, ~ exp(space_spde + Intercept),
                       n.samples = N_SAMPLES_MAP, seed = SEED + 1L)
change_intensity <- predict(change_fit, pxl, ~ log10(exp(spde_change)),
                            n.samples = N_SAMPLES_MAP, seed = SEED + 2L)

## Posterior samples of the change field, used for the evidence-threshold
## contours in Section 5 and for the distance posteriors in Sections 7-9.
change_samples <- as.matrix(
  generate(change_fit, pxl, ~ log10(exp(spde_change)),
           n.samples = N_SAMPLES_DIST, seed = SEED + 3L))

## Posterior probability of a decrease, per prediction cell.
p_decrease <- rowMeans(change_samples < 0)

message("P(decrease): min ", round(min(p_decrease), 3),
        " | median ", round(median(p_decrease), 3),
        " | max ", round(max(p_decrease), 3))

## Raster template over the prediction grid, and a helper to drop a vector of
## per-cell values onto it.
grid_template <- terra::rast(
  data.frame(x = pxl_xy[, 1], y = pxl_xy[, 2], z = change_intensity$mean),
  type = "xyz", crs = st_crs(pxl)$wkt)
cell_idx <- terra::cellFromXY(grid_template, pxl_xy)

draw_raster <- function(v) {
  vals <- rep(NA_real_, terra::ncell(grid_template))
  vals[cell_idx] <- v
  terra::setValues(grid_template, vals)
}

## THE zone definition, used everywhere: cells whose posterior probability of
## a decrease reaches the threshold. affected_region() keeps negative cells,
## so the field handed to it is thr - p(s).
zone_raster <- function(thr) draw_raster(thr - p_decrease)

## The two thresholds that are reported individually.
CRIT_HIGH <- "High-evidence zone (P >= 0.975)"
CRIT_ZERO <- "Zero contour (P = 0.50)"

THR_SIGNIF <- 0.975   # high-evidence zone; equals a 95% CI entirely below 0
THR_ZERO   <- 0.50    # zero-effect contour


# =============================================================================
# 5. Maps
# =============================================================================

as_xyz <- function(x) {
  cbind(as.data.frame(st_coordinates(x)), st_drop_geometry(x)) %>%
    rename(x = X, y = Y)
}
df_before <- as_xyz(dens_before)
df_after  <- as_xyz(dens_after)
df_change <- as_xyz(change_intensity)

bb <- st_bbox(pxl)
map_xlim <- c(bb[["xmin"]] - MAP_MARGIN_KM, bb[["xmax"]] + MAP_MARGIN_KM)
map_ylim <- c(bb[["ymin"]] - MAP_MARGIN_KM, bb[["ymax"]] + MAP_MARGIN_KM)

## Geometry stays projected; only the graticule is drawn in degrees.
coord_map <- function() {
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE,
           datum = st_crs(4326))
}

reference_layers <- function() {
  list(
    geom_sf(data = spa_area,     fill = NA, colour = "green4", linewidth = 0.6),
    geom_sf(data = conc_area,    fill = NA, colour = "grey60",  linewidth = 0.6),
    geom_sf(data = pred_outline, fill = NA, colour = "grey40", linewidth = 0.25)
  )
}

dens_limits <- range(c(dens_before$mean, dens_after$mean), na.rm = TRUE)

## PALETTE  "ylorrd" | "rocket" | "ylorbr" | "oranges" | "ylgnbu"
PALETTE     <- "ylorbr"
DENS_TRANS  <- "sqrt"
OWF_COL     <- "black"          # con rampa cálida, el rojo deja de distinguirse
DENS_BREAKS <- c(0, 0.5, 1, 2, 4, 8, 12)

dens_fill <- function() {
  brks <- DENS_BREAKS[DENS_BREAKS >= dens_limits[1] & DENS_BREAKS <= dens_limits[2]]
  args <- list(name   = expression(paste("[Ind. km"^-2 * "]")),
               limits = dens_limits, breaks = brks,
               trans  = DENS_TRANS, na.value = "grey96")
  if (PALETTE == "rocket")
    do.call(scale_fill_viridis_c, c(args, list(option = "rocket", direction = -1)))
  else
    do.call(scale_fill_distiller,
            c(args, list(direction = 1,
                         palette = switch(PALETTE,
                                          ylorrd = "YlOrRd", ylorbr = "YlOrBr",
                                          oranges = "Oranges", ylgnbu = "YlGnBu"))))
}
ggsave(p_out("fig3_density_2001_2008.png"),
       density_map(df_before, "2001-2008", 2),
       width = 8.3, height = 5.5, dpi = 300)

ggsave(p_out("fig4_density_2017_2021.png"),
       density_map(df_after, "2017-2021", 1),
       width = 8.3, height = 5.5, dpi = 300)

## Change map: solid = decrease, dashed = increase, one isoline per
## increase, grey solid = zero-effect contour.
change_max <- max(abs(df_change$mean), na.rm = TRUE)

## Evidence isolines. Both directions are drawn at the same set of
## thresholds: solid for a decrease, dashed for an increase (the areas the
## birds redistributed into). Colour is a DISCRETE scale over the thresholds
## rather than a continuous ramp, so each line reads as one evidence level;
## line type carries the direction, which keeps it to a single colour scale.
MAP_THRESHOLDS <- c(0.975, 0.95, 0.90, 0.85)

df_change$p_decrease <- p_decrease
df_change$p_increase <- 1 - p_decrease

## Labelled as probabilities, NOT as percentages. "95%" would collide with
## the 95% credible interval used to define the 0.975 zone, which is the
## single most confusing thing a reader can meet in this analysis.
## Written out explicitly rather than with format(): format() applies a common
## width to the whole vector and would render 0.95 as "0.950", and
## as.character() would render 0.90 as "0.9". Neither matches Table 2.
thr_labels <- c("0.975", "0.95", "0.90", "0.85")
stopifnot(length(thr_labels) == length(MAP_THRESHOLDS))

## Isoline colours, dark to light with the threshold.
##
## The fill is a red-blue diverging scale, the OWF outlines are black, the
## concentration area grey and the SPA green, which leaves purple as the only
## hue that is not already carrying meaning. Dark to light also puts the
## strongest line where the fill is most saturated: 0.975 is the innermost
## contour, sitting on deep red, while 0.85 is the outermost and sits on
## near-white. A reversed viridis ramp was used before and gave 0.975 a yellow
## line on a red background - the least legible line was the one the paper
## actually reports.
THRESHOLD_COLOURS <- c("#3F007D", "#6A51A3", "#9E7FD0", "#C6A9E0")

p_change <- ggplot() +
  geom_raster(data = df_change, aes(x = x, y = y, fill = mean)) +
  geom_contour(data = df_change,
               aes(x = x, y = y, z = p_decrease,
                   colour = factor(after_stat(level), levels = MAP_THRESHOLDS),
                   linetype = "Decrease"),
               breaks = MAP_THRESHOLDS, linewidth = 0.55) +
  geom_contour(data = df_change,
               aes(x = x, y = y, z = p_increase,
                   colour = factor(after_stat(level), levels = MAP_THRESHOLDS),
                   linetype = "Increase"),
               breaks = MAP_THRESHOLDS, linewidth = 0.55) +
  scale_colour_manual(values = setNames(THRESHOLD_COLOURS,
                                        as.character(MAP_THRESHOLDS)),
                      drop = FALSE, labels = thr_labels,
                      name = "P(change)",
                      guide = guide_legend(order = 2, override.aes =
                                             list(linetype = "solid",
                                                  linewidth = 0.9))) +
  scale_linetype_manual(values = c(Decrease = "solid", Increase = "22"),
                        name = NULL,
                        guide = guide_legend(order = 3)) +
  reference_layers() +
  geom_sf(data = owf, fill = NA, colour = "black") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(3, "RdBu"),
                       limits = c(-change_max, change_max),
                       name = expression(paste("log"[10], " change")),
                       guide = guide_colourbar(order = 1)) +
  coord_map() +
  labs(x = "Longitude", y = "Latitude") +
  ggtitle("Change 2017-2021 vs 2001-2008") +
  theme_bw(base_size = 14)

ggsave(p_out("fig5_change.png"), p_change,
       width = 8.3, height = 5.5, dpi = 300)


# =============================================================================
# 6. Effect distance
# =============================================================================
#
#  Follows impact_distance_and_habitat_loss.R:
#
#    dist_farms_north <- distance_raster(farms_north_sf, cellsize = 1, ...)
#    cells_north      <- raster::extract(dist_farms_north, bound_north)
#
#  i.e. build a raster of distance to the nearest OWF in the cluster, then
#  read its values along the boundary of the affected zone. Each measurement
#  is therefore one cell of the contour, valued by how far that point of the
#  contour lies from the wind farm - the outward extent of the effect.
#
#  Two changes from the original. The affected component is selected
#  automatically as the one adjacent to the cluster, instead of the
#  hand-picked p[6] / p[3]; and the same routine is applied to the
#  zero contour as well as the high-evidence zone.

## Rasterise one column of an sf point grid.
sf_to_rast <- function(x, column) {
  df <- data.frame(st_coordinates(x), z = st_drop_geometry(x)[[column]])
  names(df) <- c("x", "y", "z")
  terra::rast(df, type = "xyz", crs = st_crs(x)$wkt)
}

## The affected zone adjacent to a cluster: cells where the criterion is
## negative, dissolved, split into connected components, keeping the
## component(s) that touch the cluster. This replaces the hand-picked
## p[6] / p[3] of the original script.
affected_region <- function(r, owf_u) {
  r_neg <- terra::classify(r, matrix(c(-Inf, 0, 1, 0, Inf, NA),
                                     ncol = 3, byrow = TRUE))
  if (all(is.na(terra::values(r_neg)))) return(NULL)
  A <- suppressWarnings(
    st_cast(st_make_valid(st_set_crs(st_as_sf(terra::as.polygons(r_neg)),
                                     st_crs(owf_u))), "POLYGON"))
  hit <- lengths(st_intersects(A, owf_u)) > 0
  if (!any(hit)) return(NULL)
  st_sf(geometry = st_union(st_geometry(A[hit, ])))
}

## Raster of distance to the nearest OWF of a cluster.
##
## UNITS. terra::distance() returns metres even when the CRS is in
## kilometres, so the result is divided by 1000 to keep every distance in
## this script in km. Sanity check: the northern cluster should come out at
## roughly 7.5 km for the high-evidence zone.
make_dist_raster <- function(cluster_owf, ext_sf, cellsize = DIST_CELL_KM) {
  e <- terra::ext(terra::vect(st_geometry(ext_sf)))
  r <- terra::rast(e, resolution = cellsize, crs = st_crs(ext_sf)$wkt, vals = 0)
  terra::distance(r, terra::vect(st_geometry(cluster_owf))) / 1000
}

## Distance from every cell to the edge of the prediction area. Used to flag
## stretches of contour that merely follow the boundary of the domain: where
## the affected zone runs off the edge of the study area, the "distance to
## the contour" is really the distance to the edge of the map, not an
## ecological effect range.
EDGE_TOL_KM <- 2 * DIST_CELL_KM

## Distances along the boundary of the affected zone adjacent to a cluster:
## one value per cell of the distance surface that the boundary crosses.
## Returns NULL if there is no affected zone next to that cluster.
contour_distances <- function(r_field, owf_u, dist_rast) {
  A_c <- affected_region(r_field, owf_u)
  if (is.null(A_c)) return(NULL)
  b <- terra::vect(st_geometry(st_boundary(A_c)))
  d <- terra::extract(dist_rast,  b)[[2]]
  e <- terra::extract(dist_to_edge, b)[[2]]
  ok <- !is.na(d)
  if (!any(ok)) return(NULL)
  data.frame(distance = as.numeric(d[ok]),
             on_edge  = as.numeric(e[ok]) <= EDGE_TOL_KM)
}

owf_north_u <- st_sf(geometry = st_union(st_geometry(owf_north)))
owf_south_u <- st_sf(geometry = st_union(st_geometry(owf_south)))

dist_rast_north <- make_dist_raster(owf_north, pxl)
dist_rast_south <- make_dist_raster(owf_south, pxl)
dist_to_edge    <- make_dist_raster(st_boundary(pred_outline), pxl)

effect_distances <- function(thr) {
  r <- zone_raster(thr)
  bind_rows(
    cbind(cluster = "north", threshold = thr,
          contour_distances(r, owf_north_u, dist_rast_north)),
    cbind(cluster = "south", threshold = thr,
          contour_distances(r, owf_south_u, dist_rast_south))
  )
}

dist_all <- bind_rows(effect_distances(THR_SIGNIF), effect_distances(THR_ZERO)) %>%
  mutate(criterion = factor(
    ifelse(threshold == THR_SIGNIF, CRIT_HIGH,
                                    CRIT_ZERO),
    levels = c(CRIT_HIGH, CRIT_ZERO)))

## ci_lo / ci_hi reproduce t.test(distance)$conf.int as used in the original
## script. Note that the contour cells are contiguous and strongly
## autocorrelated, so this interval is narrower than the real uncertainty;
## Section 7 gives a posterior interval for the zero-effect contour.
## prop_on_edge is the fraction of the contour that runs along the boundary
## of the prediction area. If it is not close to zero, the affected zone does
## not close inside the study area and the mean is a lower bound on the true
## extent rather than a measurement of it. mean_km_interior excludes those
## cells.
dist_summary <- dist_all %>%
  group_by(cluster, criterion) %>%
  summarise(n            = n(),
            prop_on_edge = mean(on_edge),
            mean_km      = mean(distance),
            sd_km        = sd(distance),
            median_km    = median(distance),
            max_km       = max(distance),
            ci_lo        = t.test(distance)$conf.int[1],
            ci_hi        = t.test(distance)$conf.int[2],
            mean_km_interior = mean(distance[!on_edge]),
            n_interior       = sum(!on_edge),
            .groups   = "drop")

print(as.data.frame(dist_summary), digits = 3)
write.csv(dist_all,     p_out("effect_distances_points.csv"), row.names = FALSE)
write.csv(dist_summary, p_out("effect_distances_summary.csv"), row.names = FALSE)

## Welch t-test and Hedges' g, north vs south, per criterion.
hedges_g <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / sp * (1 - 3 / (4 * (nx + ny) - 9))
}

dist_tests <- lapply(levels(dist_all$criterion), function(cr) {
  n <- dist_all$distance[dist_all$criterion == cr & dist_all$cluster == "north"]
  s <- dist_all$distance[dist_all$criterion == cr & dist_all$cluster == "south"]
  tt <- t.test(n, s)
  data.frame(criterion = cr, mean_north = mean(n), mean_south = mean(s),
             welch_t = unname(tt$statistic), df = unname(tt$parameter),
             p_value = tt$p.value, hedges_g = hedges_g(n, s))
}) %>% bind_rows()

print(dist_tests, digits = 3)
write.csv(dist_tests, p_out("effect_distance_tests.csv"), row.names = FALSE)

# --- Figure 6 ----------------------------------------------------------------
## Reproduces the published figure, which was produced by
## impact_distance_and_habitat_loss.R as
##
##     ggbetweenstats(data = combined, x = area, y = distance) +
##       labs(x = "Area", y = "Distance") + theme_minimal(base_size = 16)
##
## ggbetweenstats defaults give the violin outline, the inner boxplot, the
## jittered points in the Dark2 palette, the mean marker with its labelled
## leader line, and the "N = " counts under each category. Its default test
## is Welch's t with Hedge's g, which is where the values quoted in the
## manuscript come from. Set results.subtitle = FALSE to drop the statistics
## from the subtitle if they are reported in the text instead.

fig6_data <- dist_all %>%
  filter(criterion == CRIT_HIGH) %>%
  mutate(area = factor(ifelse(cluster == "north", "North", "South"),
                       levels = c("North", "South")))

if (requireNamespace("ggstatsplot", quietly = TRUE)) {

  p_violin <- ggstatsplot::ggbetweenstats(
    data = fig6_data,
    x    = area,
    y    = distance
  ) +
    labs(x = "Area", y = "Distance [km]") +
    theme_minimal(base_size = 16)

} else {

  ## Same layout without the package: violin outline, inner boxplot, jittered
  ## points in the Dark2 colours, mean marker with a labelled leader line and
  ## the sample size under each category.
  message("ggstatsplot not installed - drawing the equivalent layout by hand.")

  dark2 <- c(North = "#1B9E77", South = "#D95F02")

  fig6_stats <- fig6_data %>%
    group_by(area) %>%
    summarise(n = n(), mean = mean(distance), .groups = "drop") %>%
    mutate(label = paste0("Mean = ", format(round(mean, 2), nsmall = 2)))

  x_labels <- setNames(paste0(fig6_stats$area, "\nN = ", fig6_stats$n),
                       fig6_stats$area)

  p_violin <- ggplot(fig6_data, aes(x = area, y = distance)) +
    geom_violin(fill = NA, colour = "grey20", linewidth = 0.4,
                scale = "width", trim = TRUE) +
    geom_boxplot(fill = NA, colour = "grey20", width = 0.45,
                 linewidth = 0.3, outlier.shape = NA) +
    geom_jitter(aes(colour = area), width = 0.11, height = 0,
                size = 1.5, alpha = 0.75, show.legend = FALSE) +
    geom_point(data = fig6_stats, aes(y = mean), size = 3.2,
               colour = "grey15") +
    geom_segment(data = fig6_stats,
                 aes(x = as.numeric(area), xend = as.numeric(area) + 0.34,
                     y = mean, yend = mean),
                 linetype = 3, linewidth = 0.35, colour = "grey25") +
    geom_label(data = fig6_stats,
               aes(x = as.numeric(area) + 0.35, y = mean, label = label),
               hjust = 0, size = 3.4, label.size = 0.25,
               label.padding = unit(0.15, "lines")) +
    scale_colour_manual(values = dark2) +
    scale_x_discrete(labels = x_labels, expand = expansion(add = 0.65)) +
    labs(x = "Area", y = "Distance [km]") +
    theme_minimal(base_size = 16)
}

ggsave(p_out("fig6_effect_distance.png"), p_violin,
       width = 9, height = 5.5, dpi = 300)


# =============================================================================
# 7. Zero-effect distance with posterior uncertainty
# =============================================================================
#
#  The contour cells are contiguous and not independent, so a t-test interval
#  over them understates the uncertainty. Here the whole measurement is
#  repeated inside each posterior draw and the interval comes from the
#  posterior of the mean distance. North vs South is the paired difference
#  within each draw.
#
#  Applies to the zero contour only: the high-evidence zone is a
#  property of the posterior, not of a single draw.

posterior_displacement <- function(samples, owf_u, dist_rast, label) {
  message("cluster ", label, ": ", ncol(samples), " draws")
  res <- vapply(seq_len(ncol(samples)), function(j) {
    d <- contour_distances(draw_raster(samples[, j]), owf_u, dist_rast)
    if (is.null(d)) c(NA_real_, 0, NA_real_)
    else c(mean(d$distance), 1, mean(d$on_edge))
  }, numeric(3))

  data.frame(cluster      = label,
             draw         = seq_len(ncol(samples)),
             mean_km      = res[1, ],
             has_zone     = res[2, ] == 1,
             prop_on_edge = res[3, ])
}

post_north <- posterior_displacement(change_samples, owf_north_u,
                                     dist_rast_north, "north")
post_south <- posterior_displacement(change_samples, owf_south_u,
                                     dist_rast_south, "south")
post_dist  <- bind_rows(post_north, post_south)

post_summary <- post_dist %>%
  group_by(cluster) %>%
  summarise(n_draws        = n(),
            ## P(no negative zone adjacent to this cluster in a draw)
            p_no_zone      = mean(!has_zone),
            post_mean_km   = mean(mean_km, na.rm = TRUE),
            post_median_km = median(mean_km, na.rm = TRUE),
            ci_lo_km       = unname(quantile(mean_km, 0.025, na.rm = TRUE)),
            ci_hi_km       = unname(quantile(mean_km, 0.975, na.rm = TRUE)),
            mean_prop_on_edge = mean(prop_on_edge, na.rm = TRUE),
            .groups = "drop")

print(as.data.frame(post_summary), digits = 3)
write.csv(post_dist,    p_out("zero_contour_posterior_draws.csv"), row.names = FALSE)
write.csv(post_summary, p_out("zero_contour_posterior_summary.csv"), row.names = FALSE)

delta <- post_north$mean_km - post_south$mean_km

post_contrast <- data.frame(
  comparison        = "north - south (zero-effect contour)",
  n_draws_usable    = sum(!is.na(delta)),
  post_mean_diff_km = mean(delta, na.rm = TRUE),
  ci_lo_km          = unname(quantile(delta, 0.025, na.rm = TRUE)),
  ci_hi_km          = unname(quantile(delta, 0.975, na.rm = TRUE)),
  p_north_gt_south  = mean(delta > 0, na.rm = TRUE)
)

print(post_contrast, digits = 3)
write.csv(post_contrast, p_out("zero_contour_posterior_contrast.csv"),
          row.names = FALSE)

p_post <- ggplot(post_dist %>% filter(!is.na(mean_km)),
                 aes(x = mean_km, fill = cluster)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = dist_summary %>%
               filter(criterion == CRIT_HIGH),
             aes(xintercept = mean_km, colour = cluster),
             linetype = 2, linewidth = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c(north = "#2166AC", south = "#B2182B"),
                    name = NULL) +
  scale_colour_manual(values = c(north = "#2166AC", south = "#B2182B")) +
  labs(x = "Mean distance to the zero-effect contour (km)",
       y = "Posterior density",
       caption = "Dashed lines: 95% significance-contour estimates") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(p_out("fig7_zero_contour_posterior.png"), p_post,
       width = 7.5, height = 5, dpi = 300)


# =============================================================================
# 8. Effect distance as a function of the evidence threshold
# =============================================================================
#
#  Reviewer #3 asked for a central estimate to complement the conservative
#  95% significance bound, and proposed the zero contour. With posterior
#  samples the two are the endpoints of one continuum: for each cell,
#
#      p(s) = P( delta(s) < 0 )
#
#  and the affected zone at threshold t is { s : p(s) >= t }. Then
#
#      t = 0.975  is the high-evidence zone used in Fig. 7
#      t = 0.50   is the zero contour
#
#  Reporting the series makes the trade-off explicit - how far the effect
#  reaches as the evidential requirement is relaxed - without having to
#  interpret the t = 0.50 endpoint, where by construction the evidence for a
#  decrease is nil.
#
#  prop_on_edge says how much of each contour merely follows the boundary of
#  the prediction area. A zone that does not close inside the study area is
#  not a measurement of effect range.

EVIDENCE_THRESHOLDS <- c(0.975, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50)

## Cells with p_decrease >= t are the affected zone. affected_region() keeps
## cells whose value is negative, so the field passed to it is t - p(s).
threshold_distances <- function(thr) {
  r <- zone_raster(thr)
  north <- contour_distances(r, owf_north_u, dist_rast_north)
  south <- contour_distances(r, owf_south_u, dist_rast_south)
  bind_rows(
    if (!is.null(north)) cbind(cluster = "north", threshold = thr, north),
    if (!is.null(south)) cbind(cluster = "south", threshold = thr, south)
  )
}

thr_all <- bind_rows(lapply(EVIDENCE_THRESHOLDS, threshold_distances))

## Area of the affected zone at each threshold, for context.
threshold_area <- function(thr) {
  r <- zone_raster(thr)
  bind_rows(
    data.frame(cluster = "north", threshold = thr,
               area_km2 = {
                 a <- affected_region(r, owf_north_u)
                 if (is.null(a)) NA_real_ else
                   sum(lengths(st_intersects(pxl, a)) > 0) * cell_area
               }),
    data.frame(cluster = "south", threshold = thr,
               area_km2 = {
                 a <- affected_region(r, owf_south_u)
                 if (is.null(a)) NA_real_ else
                   sum(lengths(st_intersects(pxl, a)) > 0) * cell_area
               })
  )
}
thr_area <- bind_rows(lapply(EVIDENCE_THRESHOLDS, threshold_area))

thr_summary <- thr_all %>%
  group_by(cluster, threshold) %>%
  summarise(n            = n(),
            prop_on_edge = mean(on_edge),
            mean_km      = mean(distance),
            sd_km        = sd(distance),
            median_km    = median(distance),
            max_km       = max(distance),
            mean_km_interior = mean(distance[!on_edge]),
            .groups = "drop") %>%
  left_join(thr_area, by = c("cluster", "threshold")) %>%
  arrange(cluster, desc(threshold))

print(as.data.frame(thr_summary), digits = 3)
write.csv(thr_summary, p_out("effect_distance_by_threshold.csv"),
          row.names = FALSE)

## Figure: effect distance against evidential requirement. Points where a
## large share of the contour follows the domain edge are hollow, to show
## where the estimate stops being a measurement.
p_thr <- ggplot(thr_summary,
                aes(x = threshold, y = mean_km, colour = cluster)) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = prop_on_edge > 0.1), size = 2.6, fill = "white") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21),
                     labels = c("closes inside study area",
                                ">10% of contour on domain edge"),
                     name = NULL) +
  scale_x_reverse(breaks = EVIDENCE_THRESHOLDS) +
  scale_colour_manual(values = c(north = "#2166AC", south = "#B2182B"),
                      name = NULL) +
  labs(x = "Posterior probability of a decrease required",
       y = "Mean effect distance (km)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(p_out("fig8_distance_by_threshold.png"), p_thr,
       width = 7.5, height = 5.5, dpi = 300)


# =============================================================================
# 9. Radial effect profile and effect range
# =============================================================================
#
#  The contour thresholds of Sections 6 and 8 answer a per-location question:
#  "was there a decrease here?". They do not produce an interval around the
#  reported distance, because the distance is a summary over the set of cells
#  that meet the criterion, not an estimated parameter.
#
#  This section estimates the quantity the paper is actually after: the range
#  over which the cluster depresses density. For every posterior draw and
#  every distance band from the cluster, the mean change within the band is
#  computed. The effect range of that draw is the distance at which the
#  profile returns to its reference level. Across draws this yields a
#  posterior distribution of the RANGE ITSELF, so that
#
#      posterior median  -> central estimate, comparable with the distance
#                           curves of Mendel et al., Heinanen et al.,
#                           Garthe et al. and Vilela et al. (2020)
#      5th percentile    -> "the effect extends at least X km, with 95%
#                           posterior probability" (one-sided, exact)
#
#  Two reference levels are reported. Zero is the natural null. The far-field
#  level is the mean change well away from the cluster, which is the analogue
#  of intersecting the distance curve with mean density as in Vilela et al.
#  (2020) and is robust to any overall shift in the change field.

RADIAL_BAND_KM   <- 2    # width of the distance bands
RADIAL_MAX_KM    <- 40   # how far out to profile
FARFIELD_FROM_KM <- 30   # distances beyond this define the far-field level

band_centres <- seq(RADIAL_BAND_KM / 2, RADIAL_MAX_KM - RADIAL_BAND_KM / 2,
                    by = RADIAL_BAND_KM)

## Distance from every prediction cell to each cluster.
d_north_cell <- terra::extract(dist_rast_north, pxl_xy)[, 1]
d_south_cell <- terra::extract(dist_rast_south, pxl_xy)[, 1]

## First outward crossing of a reference level, linearly interpolated between
## band centres. NA if the profile starts above the reference (no depression
## at the cluster) or never returns to it within RADIAL_MAX_KM (censored).
first_crossing <- function(y, ref) {
  ok <- which(!is.na(y))
  if (!length(ok) || y[ok[1]] >= ref) return(NA_real_)
  for (i in seq_along(y)[-1]) {
    if (is.na(y[i]) || is.na(y[i - 1])) next
    if (y[i - 1] < ref && y[i] >= ref) {
      f <- (ref - y[i - 1]) / (y[i] - y[i - 1])
      return(band_centres[i - 1] + f * (band_centres[i] - band_centres[i - 1]))
    }
  }
  NA_real_
}

radial_analysis <- function(d_cell, d_other, label) {
  ## Cells are assigned to the nearer cluster, so the two profiles do not
  ## share cells and neither cluster is credited with the other's effect.
  own  <- d_cell < d_other
  band <- cut(d_cell, breaks = seq(0, RADIAL_MAX_KM, by = RADIAL_BAND_KM),
              labels = FALSE, include.lowest = TRUE)

  idx_band <- lapply(seq_along(band_centres),
                     function(b) which(own & !is.na(band) & band == b))
  idx_far  <- which(own & d_cell > FARFIELD_FROM_KM)

  message("cluster ", label, ": ",
          paste(vapply(idx_band, length, integer(1)), collapse = "/"),
          " cells per band, ", length(idx_far), " far-field cells")

  prof <- vapply(seq_len(ncol(change_samples)), function(j) {
    v <- change_samples[, j]
    c(vapply(idx_band, function(ix) if (length(ix)) mean(v[ix]) else NA_real_,
             numeric(1)),
      mean(v[idx_far]))
  }, numeric(length(band_centres) + 1))

  nb  <- length(band_centres)
  far <- prof[nb + 1, ]

  range_zero <- vapply(seq_len(ncol(prof)),
                       function(j) first_crossing(prof[1:nb, j], 0), numeric(1))
  range_far  <- vapply(seq_len(ncol(prof)),
                       function(j) first_crossing(prof[1:nb, j], far[j]), numeric(1))

  ## Posterior summary of the profile itself, for the figure.
  prof_summary <- data.frame(
    cluster  = label,
    distance = band_centres,
    mean     = apply(prof[1:nb, , drop = FALSE], 1, mean, na.rm = TRUE),
    lo       = apply(prof[1:nb, , drop = FALSE], 1, quantile, 0.025, na.rm = TRUE),
    hi       = apply(prof[1:nb, , drop = FALSE], 1, quantile, 0.975, na.rm = TRUE)
  )

  list(
    draws = data.frame(cluster = label, draw = seq_along(range_zero),
                       range_zero = range_zero, range_far = range_far,
                       farfield = far),
    profile = prof_summary
  )
}

rad_north <- radial_analysis(d_north_cell, d_south_cell, "north")
rad_south <- radial_analysis(d_south_cell, d_north_cell, "south")

rad_draws   <- bind_rows(rad_north$draws,   rad_south$draws)
rad_profile <- bind_rows(rad_north$profile, rad_south$profile)

## p_censored: draws in which the profile never returns to the reference
## within RADIAL_MAX_KM. p_no_effect: draws with no depression at the cluster
## at all. Both must be small for the range to be a meaningful estimate.
summarise_range <- function(x, ref_label) {
  rad_draws %>%
    group_by(cluster) %>%
    summarise(
      reference    = ref_label,
      n_draws      = n(),
      p_undefined  = mean(is.na(.data[[x]])),
      median_km    = median(.data[[x]], na.rm = TRUE),
      mean_km      = mean(.data[[x]], na.rm = TRUE),
      ## one-sided 95% lower bound: the effect extends at least this far
      q05_km       = unname(quantile(.data[[x]], 0.05, na.rm = TRUE)),
      ci_lo_km     = unname(quantile(.data[[x]], 0.025, na.rm = TRUE)),
      ci_hi_km     = unname(quantile(.data[[x]], 0.975, na.rm = TRUE)),
      .groups = "drop")
}

rad_summary <- bind_rows(
  summarise_range("range_zero", "zero"),
  summarise_range("range_far",  "far-field level")
)

print(as.data.frame(rad_summary), digits = 3)
write.csv(rad_draws,   p_out("effect_range_posterior_draws.csv"), row.names = FALSE)
write.csv(rad_summary, p_out("effect_range_summary.csv"), row.names = FALSE)
write.csv(rad_profile, p_out("radial_profile.csv"), row.names = FALSE)

## Figure: radial profile with 95% credible band, and the estimated range.
rad_lines <- rad_summary %>% filter(reference == "far-field level")

p_radial <- ggplot(rad_profile, aes(x = distance)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = cluster), alpha = 0.25) +
  geom_line(aes(y = mean, colour = cluster), linewidth = 0.8) +
  geom_vline(data = rad_lines, aes(xintercept = median_km, colour = cluster),
             linetype = 2, linewidth = 0.6, show.legend = FALSE) +
  geom_vline(data = rad_lines, aes(xintercept = q05_km, colour = cluster),
             linetype = 3, linewidth = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c(north = "#2166AC", south = "#B2182B"), name = NULL) +
  scale_colour_manual(values = c(north = "#2166AC", south = "#B2182B"), name = NULL) +
  labs(x = "Distance from OWF cluster (km)",
       y = expression(paste("Mean ", log[10], " change")),
       caption = paste("Dashed: posterior median effect range.",
                       "Dotted: 5th percentile (effect extends at least this far,",
                       "95% posterior probability).")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(p_out("fig9_radial_effect_profile.png"), p_radial,
       width = 7.5, height = 5.5, dpi = 300)


# =============================================================================
# 10. Individuals displaced within each affected zone
# =============================================================================
#
#  Abundance in the after period minus abundance in the before period, summed
#  over the affected zone. This is what the Methods section describes and
#  corresponds to the `diff <- aft - bef` / cellStats(sum) route in
#  impact_distance_and_habitat_loss.R.
#
#  Note that the other route in that script, sum(exp(pred_change_975)), does
#  not give individuals: the raster is on a log10 scale, so the multiplicative
#  change would be 10^x, and it is a ratio in any case, not a count.

displaced_individuals <- function(thr, owf_u, label, criterion) {
  A_c <- affected_region(zone_raster(thr), owf_u)
  if (is.null(A_c)) return(NULL)

  idx      <- lengths(st_intersects(pxl, A_c)) > 0
  area_km2 <- sum(idx) * cell_area
  n_before <- sum(dens_before$mean[idx]) * cell_area
  n_after  <- sum(dens_after$mean[idx])  * cell_area

  data.frame(cluster     = label,
             criterion   = criterion,
             area_km2    = area_km2,
             n_before    = n_before,
             n_after     = n_after,
             n_displaced = n_before - n_after,
             density_change_per_km2 = (n_after - n_before) / area_km2)
}

displaced <- bind_rows(
  displaced_individuals(THR_SIGNIF, owf_north_u, "north", CRIT_HIGH),
  displaced_individuals(THR_SIGNIF, owf_south_u, "south", CRIT_HIGH),
  displaced_individuals(THR_ZERO,   owf_north_u, "north", CRIT_ZERO),
  displaced_individuals(THR_ZERO,   owf_south_u, "south", CRIT_ZERO)
)

print(displaced, digits = 4)
write.csv(displaced, p_out("displaced_individuals.csv"), row.names = FALSE)



#------------
#Density population stats before-after in the SPA area
spa_stats <- function(d) {
  pts <- st_as_sf(d, coords = c("x", "y"), crs = main_crs)
  v <- d$mean[lengths(st_intersects(pts, spa_area)) > 0]
  c(mean = mean(v, na.rm = TRUE), median = median(v, na.rm = TRUE),
    max = max(v, na.rm = TRUE), cells = sum(!is.na(v)))
}
round(rbind(before = spa_stats(df_before), after = spa_stats(df_after)), 2)


#---------
zstat <- function(d, z) {
  pts <- st_as_sf(d, coords = c("x","y"), crs = main_crs)
  v <- d$mean[lengths(st_intersects(pts, z)) > 0]
  c(mean = mean(v, na.rm=TRUE), median = median(v, na.rm=TRUE), max = max(v, na.rm=TRUE))
}
round(rbind(
  spa_before   = zstat(df_before, spa_area),   spa_after   = zstat(df_after, spa_area),
  mca_before   = zstat(df_before, conc_area),  mca_after   = zstat(df_after, conc_area),
  study_before = c(mean(df_before$mean, na.rm=TRUE), median(df_before$mean, na.rm=TRUE), max(df_before$mean, na.rm=TRUE)),
  study_after  = c(mean(df_after$mean,  na.rm=TRUE), median(df_after$mean,  na.rm=TRUE), max(df_after$mean,  na.rm=TRUE))
), 2)


# -------------

#Inside North of helgoland
unique(owf$cluster)

hel <- owf %>% filter(cluster == "Helgoland") %>% st_geometry() %>% st_union() %>% st_sf()
rbind(before = zstat(df_before, st_buffer(hel, 5)),
      after  = zstat(df_after,  st_buffer(hel, 5)))   # 5 km alrededor, el CRS ya está en km

# =============================================================================
# 11. Session info
# =============================================================================

writeLines(capture.output(sessionInfo()), p_out("sessionInfo.txt"))
