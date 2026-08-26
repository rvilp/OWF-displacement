###############################################################################
#  Sensitivity of the result to the composition of the baseline period
#
#  Varying the baseline changes the data realisation rather than the model, 
#  so what this measures is interannualvariability in diver distribution, 
#  not sensitivity of the joint-likelihood formulation. 
#  It is nonetheless cheap to run and informative either way, so it is reported.
#
#  The baseline holds six spring seasons: 2001-2005 and 2008 (no data for
#  2006 or 2007). Effort is very unevenly spread across them:
#
#      2001  3 616 km2      2004  15 406 km2
#      2002  6 759 km2      2005   3 625 km2
#      2003  4 952 km2      2008  23 053 km2
#
#  2008 alone is a third of the baseline effort and 2004 another fifth, so
#  dropping either is the strongest available perturbation. 2003 is the year
#  carrying the strong aggregation referred to in the Results.
#
#  CONTROL ROW. The last variant repeats the reference configuration with a
#  different sampling seed. It changes no data at all, so the spread between
#  it and the reference is the Monte Carlo noise of the measurement. A
#  variant that moves the estimate by less than that has not moved it.
#
#      source("R/sensitivity_baseline_period.R")
###############################################################################

library(INLA); library(inlabru); library(fmesher)
library(sf); library(terra); library(dplyr); library(ggplot2)

DATA_DIR <- "data"; OUT_DIR <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)

sf_use_s2(FALSE)
main_crs <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")

## ---- must match joint_likelihood_analysis.R exactly ------------------------
BEFORE_PHASE_MAX <- 8
AFTER_PHASE_MIN  <- 17
DETECTION <- c(HiDef = 1, DAISI = exp(-0.10),
               APEM = exp(-0.40), conventional = exp(-0.22))
PRIOR_SIGMA <- c(0.2, 0.01)
PRIOR_RANGE <- c(15, NA)
N_SAMPLES_DIST <- 10000
SEED <- 20260811L
DIST_CELL_KM <- 1
THR_SIGNIF <- 0.975
EDGE_TOL_KM <- 2 * DIST_CELL_KM
CLUSTER_NORTH <- c("West of Sylt", "Butendiek")
CLUSTER_SOUTH <- c("UMBO")
## ---------------------------------------------------------------------------

## phase = year - 2000. Baseline years present: 1,2,3,4,5,8.
VARIANTS <- list(
  list(tag = "ref",        label = "Reference (2001-2008, all)",   keep = c(1,2,3,4,5,8), seed = SEED),
  list(tag = "no2003",     label = "Excluding 2003",               keep = c(1,2,4,5,8),   seed = SEED),
  list(tag = "no2008",     label = "Excluding 2008",               keep = c(1,2,3,4,5),   seed = SEED),
  list(tag = "no2004",     label = "Excluding 2004",               keep = c(1,2,3,5,8),   seed = SEED),
  list(tag = "early",      label = "Early baseline (2001-2004)",   keep = c(1,2,3,4),     seed = SEED),
  list(tag = "late",       label = "Late baseline (2005-2008)",    keep = c(5,8),         seed = SEED),
  list(tag = "higheffort", label = "High-effort years only (2004, 2008)", keep = c(4,8),  seed = SEED),
  list(tag = "ref_seed2",  label = "Reference, different seed (control)", keep = c(1,2,3,4,5,8),
       seed = SEED + 977L)
)

# ---------------------------------------------------------------- data
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
dat$method_chr <- as.character(dat$method)
stopifnot(all(dat$method_chr %in% names(DETECTION)))
dat$E_eff <- dat$area * unname(DETECTION[dat$method_chr])

after <- dat %>% filter(phase >= AFTER_PHASE_MIN)

owf <- farms %>% group_by(Name, cluster) %>% summarise(.groups = "drop")
owf_u <- list(
  north = st_sf(geometry = st_union(st_geometry(
    owf %>% filter(cluster %in% CLUSTER_NORTH)))),
  south = st_sf(geometry = st_union(st_geometry(
    owf %>% filter(cluster %in% CLUSTER_SOUTH)))))
owf_raw <- list(north = owf %>% filter(cluster %in% CLUSTER_NORTH),
                south = owf %>% filter(cluster %in% CLUSTER_SOUTH))

# ------------------------------------------------- measurement, as in Sec. 6
grid_template <- terra::rast(
  data.frame(x = pxl_xy[, 1], y = pxl_xy[, 2], z = 0),
  type = "xyz", crs = main_crs$wkt)
cell_idx  <- terra::cellFromXY(grid_template, pxl_xy)
cell_area <- prod(terra::res(grid_template))

draw_raster <- function(v) {
  vals <- rep(NA_real_, terra::ncell(grid_template))
  vals[cell_idx] <- v
  terra::setValues(grid_template, vals)
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

## NOTE the extent: the bounding box of the prediction POINTS, exactly as
## make_dist_raster() in joint_likelihood_analysis.R. Using the raster extent
## instead shifts the 1 km grid by half a prediction cell and moves the
## reported mean by about 0.1 km.
make_dist_raster <- function(cluster_owf) {
  e <- terra::ext(terra::vect(st_geometry(pxl)))
  r <- terra::rast(e, resolution = DIST_CELL_KM, crs = main_crs$wkt, vals = 0)
  terra::distance(r, terra::vect(st_geometry(cluster_owf))) / 1000
}
dist_rast <- lapply(owf_raw, make_dist_raster)
dist_to_edge <- make_dist_raster(st_boundary(pred_outline))

measure <- function(p_decrease, cl) {
  A_c <- affected_region(draw_raster(THR_SIGNIF - p_decrease), owf_u[[cl]])
  if (is.null(A_c)) return(NULL)
  b <- terra::vect(st_geometry(st_boundary(A_c)))
  d <- terra::extract(dist_rast[[cl]], b)[[2]]
  e <- terra::extract(dist_to_edge,    b)[[2]]
  ok <- !is.na(d)
  if (!any(ok)) return(NULL)
  data.frame(cluster = cl, n = sum(ok),
             prop_on_edge = mean(as.numeric(e[ok]) <= EDGE_TOL_KM),
             mean_km = mean(d[ok]), sd_km = sd(d[ok]),
             area_km2 = sum(lengths(st_intersects(pxl, A_c)) > 0) * cell_area)
}

# ---------------------------------------------------------------- the loop
matern <- inla.spde2.pcmatern(spatial_mesh, prior.sigma = PRIOR_SIGMA,
                              prior.range = PRIOR_RANGE)
cmp_joint <- ~ spde_before(geometry, model = matern) +
               spde_change(geometry, model = matern) +
               Intercept(1) + Intercept_after(1)
inla_opts <- list(
  control.compute   = list(dic = TRUE, cpo = TRUE, waic = TRUE),
  control.inla      = list(int.strategy = "ccd"),
  control.predictor = list(compute = TRUE,
                           quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)))

run_variant <- function(v) {
  bef <- dat %>% filter(phase %in% v$keep)
  cat("\n=====================================================\n")
  cat(v$label, "\n  years: ", paste(2000 + sort(v$keep), collapse = ", "),
      "\n  rows: ", nrow(bef), "  effort: ", round(sum(bef$area)), " km2\n", sep = "")

  f <- p_out(paste0("sens_", v$tag, ".rds"))
  if (file.exists(f)) {
    cat("  loading cached fit\n"); fit <- readRDS(f)
  } else {
    t0 <- Sys.time()
    fit <- bru(cmp_joint,
      bru_obs("nbinomial", data = bef,   formula = NHAT ~ spde_before + Intercept,
              E = bef$E_eff),
      bru_obs("nbinomial", data = after,
              formula = NHAT ~ spde_before + spde_change + Intercept_after,
              E = after$E_eff),
      options = c(inla_opts, list(bru_max_iter = 1)))
    saveRDS(fit, f)
    cat("  fitted in ", round(difftime(Sys.time(), t0, units = "mins"), 1),
        " min\n", sep = "")
  }

  s <- as.matrix(generate(fit, pxl, ~ log10(exp(spde_change)),
                          n.samples = N_SAMPLES_DIST, seed = v$seed + 3L))
  p_dec <- rowMeans(s < 0)

  out <- bind_rows(lapply(names(owf_u), function(cl) measure(p_dec, cl)))
  if (!is.null(out) && nrow(out)) {
    out$variant <- v$label; out$tag <- v$tag
    out$n_years <- length(v$keep)
    out$baseline_effort_km2 <- sum(bef$area)
    print(as.data.frame(out[, c("cluster", "n", "mean_km", "sd_km",
                                "prop_on_edge", "area_km2")]),
          digits = 4, row.names = FALSE)
  }
  out
}

res <- bind_rows(lapply(VARIANTS, run_variant)) %>%
  select(variant, tag, cluster, n_years, baseline_effort_km2,
         n, prop_on_edge, mean_km, sd_km, area_km2)

write.csv(res, p_out("sensitivity_baseline.csv"), row.names = FALSE)

# ---------------------------------------------------------------- report
ref   <- res %>% filter(tag == "ref")   %>% select(cluster, ref_km = mean_km)
ctrl  <- res %>% filter(tag == "ref_seed2") %>% select(cluster, ctrl_km = mean_km)
noise <- left_join(ref, ctrl, by = "cluster") %>%
  mutate(mc_noise_km = abs(ctrl_km - ref_km))

tab <- res %>% left_join(ref, by = "cluster") %>%
  left_join(noise %>% select(cluster, mc_noise_km), by = "cluster") %>%
  mutate(diff_km = mean_km - ref_km,
         beyond_noise = abs(diff_km) > mc_noise_km)

cat("\n\n================= SUMMARY =================\n")
print(as.data.frame(tab %>%
        select(variant, cluster, n, mean_km, diff_km, mc_noise_km,
               beyond_noise, area_km2)), digits = 4, row.names = FALSE)
cat("\nMonte Carlo noise floor (reference vs reference, different seed):\n")
print(as.data.frame(noise), digits = 4, row.names = FALSE)
cat("\nA variant whose diff_km is smaller than mc_noise_km has not moved the\n")
cat("estimate by more than the measurement's own sampling noise.\n")
write.csv(tab, p_out("sensitivity_baseline_summary.csv"), row.names = FALSE)

p <- ggplot(tab, aes(x = reorder(variant, mean_km), y = mean_km, colour = cluster)) +
  geom_hline(data = ref, aes(yintercept = ref_km, colour = cluster),
             linetype = 2, linewidth = 0.5, show.legend = FALSE) +
  geom_point(size = 2.6) +
  coord_flip() +
  scale_colour_manual(values = c(north = "#2166AC", south = "#B2182B"), name = NULL) +
  labs(x = NULL, y = "Mean effect distance (km)",
       caption = "Dashed: reference baseline. Bottom row repeats the reference with a different seed.") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")

ggsave(p_out("figS_baseline_sensitivity.png"), p, width = 8, height = 5, dpi = 300)
cat("\nwritten: outputs/sensitivity_baseline.csv, _summary.csv,",
    "figS_baseline_sensitivity.png\n")
