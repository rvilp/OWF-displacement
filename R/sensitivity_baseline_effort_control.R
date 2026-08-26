###############################################################################
#  Control for the baseline sensitivity: is it the years, or the effort?
#
#  Dropping baseline years moves the estimate downwards, and in the northern
#  cluster the size of the shift tracks the amount of survey effort removed
#  almost exactly (r = 0.975):
#
#      excluding 2003     -8.6% effort    -0.28 km   (within sampling noise)
#      excluding 2004    -26.8% effort    -1.30 km
#      high-effort only  -33.0% effort    -2.42 km
#      excluding 2008    -40.2% effort    -2.11 km
#      early 2001-2004   -46.5% effort    -2.55 km
#      late 2005-2008    -53.5% effort    -3.13 km
#
#  Two mechanisms could produce this and they have opposite implications.
#
#    (a) PRECISION. The affected zone is defined by an evidential threshold,
#        P(decrease) >= 0.975. Less baseline data widens the posterior, fewer
#        cells clear the threshold, the zone shrinks and the mean distance
#        falls. On this reading the reported distance is a lower bound that
#        behaves exactly as a lower bound should, and the sensitivity is a
#        property of the criterion rather than instability in the signal.
#
#    (b) COMPOSITION. Individual baseline years genuinely describe different
#        diver distributions, and the result depends on which ones are in.
#        On this reading the reported distance is contingent on the baseline
#        chosen, which is the Reviewer's concern.
#
#  This script separates them. It keeps ALL SIX baseline years and randomly
#  removes observations to reach the same effort levels, with several
#  replicates per level so the result does not depend on which rows are
#  dropped. If the year-dropping variants fall on the curve traced by random
#  subsampling, the mechanism is (a). If they fall below it, part of the
#  effect is (b) and must be reported as such.
#
#  Run sensitivity_baseline_period.R first.
#
#      source("R/sensitivity_baseline_effort_control.R")
###############################################################################

library(INLA); library(inlabru); library(fmesher)
library(sf); library(terra); library(dplyr); library(ggplot2)

DATA_DIR <- "data"; OUT_DIR <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)

sf_use_s2(FALSE)
main_crs <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")

BEFORE_PHASE_MAX <- 8
AFTER_PHASE_MIN  <- 17
DETECTION <- c(HiDef = 1, DAISI = exp(-0.10),
               APEM = exp(-0.40), conventional = exp(-0.22))
PRIOR_SIGMA <- c(0.2, 0.01); PRIOR_RANGE <- c(15, NA)
N_SAMPLES_DIST <- 10000        # raise together with the main script
SEED <- 20260811L
DIST_CELL_KM <- 1; THR_SIGNIF <- 0.975; EDGE_TOL_KM <- 2 * DIST_CELL_KM
CLUSTER_NORTH <- c("West of Sylt", "Butendiek"); CLUSTER_SOUTH <- c("UMBO")

## Effort fractions matching the year-dropping variants, and replicates per
## level. Row selection is seeded, so the whole thing is reproducible.
EFFORT_LEVELS <- c(0.914, 0.732, 0.670, 0.598, 0.535, 0.465)
N_REPLICATES  <- 5

as_km <- function(x) { y <- if (inherits(x,"sf")) x else st_as_sf(x); st_crs(y) <- main_crs; y }
D <- new.env(parent = emptyenv()); load(p_data("diver_owf_data.RData"), envir = D)

countdata    <- as_km(D$counts_5km); farms <- as_km(D$owf_polygons)
pxl          <- as_km(D$prediction_pxl); pred_outline <- as_km(D$prediction_mask)
spatial_mesh <- fm_as_mesh_2d(D$mesh_5km); fm_crs(spatial_mesh) <- main_crs

if (!all(as.character(st_geometry_type(pxl)) == "POINT"))
  pxl <- suppressWarnings(st_centroid(pxl))
pxl <- suppressWarnings(st_filter(pxl, st_union(st_geometry(pred_outline))))
pxl_xy <- st_coordinates(pxl)

dat <- countdata %>% filter(!is.na(area), !is.na(NHAT), area > 0) %>%
  mutate(NHAT = round(NHAT))
dat$E_eff <- dat$area * unname(DETECTION[as.character(dat$method)])
before_full <- dat %>% filter(phase <= BEFORE_PHASE_MAX)
after       <- dat %>% filter(phase >= AFTER_PHASE_MIN)
total_effort <- sum(before_full$area)

owf <- farms %>% group_by(Name, cluster) %>% summarise(.groups = "drop")
owf_u <- list(north = st_sf(geometry = st_union(st_geometry(
                 owf %>% filter(cluster %in% CLUSTER_NORTH)))),
              south = st_sf(geometry = st_union(st_geometry(
                 owf %>% filter(cluster %in% CLUSTER_SOUTH)))))
owf_raw <- list(north = owf %>% filter(cluster %in% CLUSTER_NORTH),
                south = owf %>% filter(cluster %in% CLUSTER_SOUTH))

grid_template <- terra::rast(data.frame(x = pxl_xy[,1], y = pxl_xy[,2], z = 0),
                             type = "xyz", crs = main_crs$wkt)
cell_idx  <- terra::cellFromXY(grid_template, pxl_xy)
cell_area <- prod(terra::res(grid_template))
draw_raster <- function(v) { vv <- rep(NA_real_, terra::ncell(grid_template))
  vv[cell_idx] <- v; terra::setValues(grid_template, vv) }

affected_region <- function(r, o) {
  rn <- terra::classify(r, matrix(c(-Inf,0,1, 0,Inf,NA), ncol = 3, byrow = TRUE))
  if (all(is.na(terra::values(rn)))) return(NULL)
  A <- suppressWarnings(st_cast(st_make_valid(
    st_set_crs(st_as_sf(terra::as.polygons(rn)), main_crs)), "POLYGON"))
  hit <- lengths(st_intersects(A, o)) > 0
  if (!any(hit)) return(NULL)
  st_sf(geometry = st_union(st_geometry(A[hit, ])))
}
make_dist_raster <- function(g) {
  r <- terra::rast(terra::ext(terra::vect(st_geometry(pxl))),
                   resolution = DIST_CELL_KM, crs = main_crs$wkt, vals = 0)
  terra::distance(r, terra::vect(st_geometry(g))) / 1000
}
dist_rast <- lapply(owf_raw, make_dist_raster)
dist_to_edge <- make_dist_raster(st_boundary(pred_outline))

measure <- function(p_dec, cl) {
  A <- affected_region(draw_raster(THR_SIGNIF - p_dec), owf_u[[cl]])
  if (is.null(A)) return(NULL)
  b <- terra::vect(st_geometry(st_boundary(A)))
  d <- terra::extract(dist_rast[[cl]], b)[[2]]
  e <- terra::extract(dist_to_edge,    b)[[2]]
  ok <- !is.na(d); if (!any(ok)) return(NULL)
  data.frame(cluster = cl, n = sum(ok), mean_km = mean(d[ok]),
             prop_on_edge = mean(as.numeric(e[ok]) <= EDGE_TOL_KM),
             area_km2 = sum(lengths(st_intersects(pxl, A)) > 0) * cell_area)
}

matern <- inla.spde2.pcmatern(spatial_mesh, prior.sigma = PRIOR_SIGMA,
                              prior.range = PRIOR_RANGE)
cmp <- ~ spde_before(geometry, model = matern) +
         spde_change(geometry, model = matern) + Intercept(1) + Intercept_after(1)
opts <- list(control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
             control.inla = list(int.strategy = "ccd"),
             control.predictor = list(compute = TRUE,
               quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)))

## Drop whole observations at random until the retained effort is at or below
## the target fraction. Rows are the natural unit: they are what a shorter
## survey season would remove.
subsample <- function(frac, rep) {
  set.seed(1000L * rep + round(frac * 1000))
  ord <- sample(nrow(before_full))
  keep <- ord[cumsum(before_full$area[ord]) <= frac * total_effort]
  before_full[sort(keep), ]
}

run_one <- function(frac, rep) {
  bef <- subsample(frac, rep)
  f <- p_out(sprintf("sens_eff_%03d_r%d.rds", round(frac * 100), rep))
  if (file.exists(f)) fit <- readRDS(f) else {
    fit <- bru(cmp,
      bru_obs("nbinomial", data = bef, formula = NHAT ~ spde_before + Intercept,
              E = bef$E_eff),
      bru_obs("nbinomial", data = after,
              formula = NHAT ~ spde_before + spde_change + Intercept_after,
              E = after$E_eff),
      options = c(opts, list(bru_max_iter = 1)))
    saveRDS(fit, f)
  }
  s <- as.matrix(generate(fit, pxl, ~ log10(exp(spde_change)),
                          n.samples = N_SAMPLES_DIST, seed = SEED + 3L))
  p_dec <- rowMeans(s < 0)
  out <- bind_rows(lapply(names(owf_u), function(cl) measure(p_dec, cl)))
  if (is.null(out) || !nrow(out)) return(NULL)
  out$effort_frac <- frac; out$replicate <- rep
  out$effort_km2 <- sum(bef$area); out$n_rows <- nrow(bef)
  cat(sprintf("  frac %.3f rep %d : %s\n", frac, rep,
              paste(sprintf("%s %.2f km", out$cluster, out$mean_km), collapse = " | ")))
  out
}

cat("Random effort subsampling, all six baseline years retained\n")
res <- bind_rows(lapply(EFFORT_LEVELS, function(fr) {
  cat(sprintf("\neffort fraction %.3f\n", fr))
  bind_rows(lapply(seq_len(N_REPLICATES), function(r) run_one(fr, r)))
}))

write.csv(res, p_out("sensitivity_effort_control.csv"), row.names = FALSE)

summ <- res %>% group_by(cluster, effort_frac) %>%
  summarise(n_rep = n(), mean_km = mean(mean_km), sd_km = sd(mean_km),
            lo = min(mean_km), hi = max(mean_km), .groups = "drop")
cat("\n=========== random subsampling, by effort level ===========\n")
print(as.data.frame(summ), digits = 4, row.names = FALSE)

## Overlay the year-dropping variants, if they are available.
yr_file <- p_out("sensitivity_baseline_summary.csv")
if (file.exists(yr_file)) {
  yr <- read.csv(yr_file) %>%
    filter(!tag %in% c("ref_seed2")) %>%
    mutate(effort_frac = baseline_effort_km2 / total_effort)
  p <- ggplot() +
    geom_ribbon(data = summ, aes(x = effort_frac, ymin = lo, ymax = hi,
                                 fill = cluster), alpha = 0.2) +
    geom_line(data = summ, aes(x = effort_frac, y = mean_km, colour = cluster),
              linewidth = 0.8) +
    geom_point(data = yr, aes(x = effort_frac, y = mean_km, colour = cluster),
               shape = 17, size = 3) +
    scale_colour_manual(values = c(north = "#2166AC", south = "#B2182B"), name = NULL) +
    scale_fill_manual(values = c(north = "#2166AC", south = "#B2182B"), guide = "none") +
    labs(x = "Fraction of baseline survey effort retained",
         y = "Mean effect distance (km)",
         caption = paste("Lines and bands: all six baseline years retained,",
                         "effort removed at random (range over replicates).",
                         "\nTriangles: whole baseline years dropped.")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  ggsave(p_out("figS_effort_control.png"), p, width = 8, height = 5.5, dpi = 300)

  cat("\n=========== year-dropping vs random subsampling ===========\n")
  ## rule = 1, not 2. The reference sits at effort_frac = 1, outside the
  ## subsampling range, and rule = 2 would silently clamp it to the value at
  ## the highest level sampled (0.914) and report a spurious excess for a
  ## variant that removed nothing at all.
  cmpn <- yr %>% rowwise() %>%
    mutate(random_at_same_effort = {
      s <- summ[summ$cluster == cluster, ]
      approx(s$effort_frac, s$mean_km, xout = effort_frac, rule = 1)$y }) %>%
    ungroup() %>%
    mutate(excess_km = mean_km - random_at_same_effort) %>%
    filter(!is.na(random_at_same_effort)) %>%
    select(variant, cluster, effort_frac, mean_km, random_at_same_effort, excess_km)
  print(as.data.frame(cmpn), digits = 4, row.names = FALSE)
  cat("\nexcess_km near zero => the shift is explained by how much evidence was\n")
  cat("removed, not by which years. Systematically negative => the identity of\n")
  cat("the excluded years matters beyond the loss of effort, and that has to be\n")
  cat("reported.\n")
  write.csv(cmpn, p_out("sensitivity_effort_vs_years.csv"), row.names = FALSE)
}
