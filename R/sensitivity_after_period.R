###############################################################################
#  The same test applied to the AFTER period
#
#  Dropping baseline years moves the estimate downwards, and the size of the
#  shift tracks the survey effort removed NEAR THE CLUSTER CONCERNED much
#  better than the total effort removed (south r = 0.85 against r = 0.45).
#  That points at spatial coverage rather than interannual variability: the
#  baseline is unevenly covered, with 2001 and 2002 entirely northern,
#  2005 entirely non-northern, and only 2004 and 2008 covering the whole area.
#
#  After-period effort by latitude band (km2):
#
#      year   north  middle  south   total
#      2017    4592    2574   2320    9486
#      2018    4701    3808   3237   11745
#      2019    3056    2914    685    6655
#      2020     664    1573    470    2708
#      2021     427     749      0    1177
#
#  N_SAMPLES_DIST must match whatever the main script uses. 
#
#      source("R/sensitivity_after_period.R")
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
N_SAMPLES_DIST <- 10000          # raise together with the main script
SEED <- 20260811L
DIST_CELL_KM <- 1; THR_SIGNIF <- 0.975; EDGE_TOL_KM <- 2 * DIST_CELL_KM
CLUSTER_NORTH <- c("West of Sylt", "Butendiek"); CLUSTER_SOUTH <- c("UMBO")

## phase = year - 2000. After-period years present: 17..21.
VARIANTS <- list(
  list(tag = "a_ref",     label = "Reference (2017-2021, all)",  keep = 17:21, seed = SEED),
  list(tag = "a_no2017",  label = "Excluding 2017",              keep = c(18,19,20,21), seed = SEED),
  list(tag = "a_no2018",  label = "Excluding 2018",              keep = c(17,19,20,21), seed = SEED),
  list(tag = "a_no2019",  label = "Excluding 2019",              keep = c(17,18,20,21), seed = SEED),
  list(tag = "a_notail",  label = "Excluding 2020-2021 (small years)", keep = c(17,18,19), seed = SEED),
  list(tag = "a_first2",  label = "First two years only (2017-2018)",  keep = c(17,18),   seed = SEED),
  list(tag = "a_last3",   label = "Last three years only (2019-2021)", keep = c(19,20,21), seed = SEED),
  list(tag = "a_ref_s2",  label = "Reference, different seed (control)", keep = 17:21, seed = SEED + 977L)
)

as_km <- function(x) { y <- if (inherits(x,"sf")) x else st_as_sf(x); st_crs(y) <- main_crs; y }
D <- new.env(parent = emptyenv()); load(p_data("diver_owf_data.RData"), envir = D)

countdata <- as_km(D$counts_5km); farms <- as_km(D$owf_polygons)
pxl <- as_km(D$prediction_pxl); pred_outline <- as_km(D$prediction_mask)
spatial_mesh <- fm_as_mesh_2d(D$mesh_5km); fm_crs(spatial_mesh) <- main_crs
if (!all(as.character(st_geometry_type(pxl)) == "POINT"))
  pxl <- suppressWarnings(st_centroid(pxl))
pxl <- suppressWarnings(st_filter(pxl, st_union(st_geometry(pred_outline))))
pxl_xy <- st_coordinates(pxl)

dat <- countdata %>% filter(!is.na(area), !is.na(NHAT), area > 0) %>%
  mutate(NHAT = round(NHAT))
dat$E_eff <- dat$area * unname(DETECTION[as.character(dat$method)])
before <- dat %>% filter(phase <= BEFORE_PHASE_MAX)

owf <- farms %>% group_by(Name, cluster) %>% summarise(.groups = "drop")
owf_u <- list(north = st_sf(geometry = st_union(st_geometry(
                 owf %>% filter(cluster %in% CLUSTER_NORTH)))),
              south = st_sf(geometry = st_union(st_geometry(
                 owf %>% filter(cluster %in% CLUSTER_SOUTH)))))
owf_raw <- list(north = owf %>% filter(cluster %in% CLUSTER_NORTH),
                south = owf %>% filter(cluster %in% CLUSTER_SOUTH))

## Effort within LOCAL_KM of each cluster, so the summary can report how much
## of each cluster's own after-period coverage a variant removes. This is the
## quantity that predicted the baseline result, so it is measured directly
## here rather than approximated by latitude bands.
LOCAL_KM <- 30
local_effort <- function(rows, cl) {
  if (!nrow(rows)) return(0)
  d <- as.numeric(st_distance(rows, owf_u[[cl]]))
  sum(rows$area[d <= LOCAL_KM])
}

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
  data.frame(cluster = cl, n = sum(ok), mean_km = mean(d[ok]), sd_km = sd(d[ok]),
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

after_full <- dat %>% filter(phase >= AFTER_PHASE_MIN)
loc_ref <- sapply(names(owf_u), function(cl) local_effort(after_full, cl))

run_variant <- function(v) {
  aft <- dat %>% filter(phase %in% v$keep)
  cat("\n=====================================================\n")
  cat(v$label, "\n  years: ", paste(2000 + sort(v$keep), collapse = ", "),
      "\n  rows: ", nrow(aft), "  effort: ", round(sum(aft$area)), " km2\n", sep = "")

  f <- p_out(paste0("sens_", v$tag, ".rds"))
  if (file.exists(f)) { cat("  loading cached fit\n"); fit <- readRDS(f) } else {
    t0 <- Sys.time()
    fit <- bru(cmp,
      bru_obs("nbinomial", data = before, formula = NHAT ~ spde_before + Intercept,
              E = before$E_eff),
      bru_obs("nbinomial", data = aft,
              formula = NHAT ~ spde_before + spde_change + Intercept_after,
              E = aft$E_eff),
      options = c(opts, list(bru_max_iter = 1)))
    saveRDS(fit, f)
    cat("  fitted in ", round(difftime(Sys.time(), t0, units = "secs")), " s\n", sep = "")
  }

  s <- as.matrix(generate(fit, pxl, ~ log10(exp(spde_change)),
                          n.samples = N_SAMPLES_DIST, seed = v$seed + 3L))
  out <- bind_rows(lapply(names(owf_u), function(cl) {
    m <- measure(rowMeans(s < 0), cl)
    if (is.null(m)) return(NULL)
    m$local_effort_kept <- local_effort(aft, cl) / loc_ref[[cl]]
    m
  }))
  if (!is.null(out) && nrow(out)) {
    out$variant <- v$label; out$tag <- v$tag
    out$after_effort_km2 <- sum(aft$area)
    print(as.data.frame(out[, c("cluster","n","mean_km","local_effort_kept",
                                "prop_on_edge","area_km2")]),
          digits = 4, row.names = FALSE)
  }
  out
}

res <- bind_rows(lapply(VARIANTS, run_variant))
write.csv(res, p_out("sensitivity_after.csv"), row.names = FALSE)

ref  <- res %>% filter(tag == "a_ref")    %>% select(cluster, ref_km = mean_km)
ctrl <- res %>% filter(tag == "a_ref_s2") %>% select(cluster, ctrl_km = mean_km)
noise <- left_join(ref, ctrl, by = "cluster") %>%
  mutate(mc_noise_km = abs(ctrl_km - ref_km))

tab <- res %>% left_join(ref, by = "cluster") %>%
  left_join(noise %>% select(cluster, mc_noise_km), by = "cluster") %>%
  mutate(diff_km = mean_km - ref_km, beyond_noise = abs(diff_km) > mc_noise_km)

cat("\n\n================= SUMMARY (after period) =================\n")
print(as.data.frame(tab %>% select(variant, cluster, n, mean_km, diff_km,
                                   local_effort_kept, mc_noise_km, beyond_noise,
                                   area_km2)), digits = 4, row.names = FALSE)
write.csv(tab, p_out("sensitivity_after_summary.csv"), row.names = FALSE)

for (cl in unique(tab$cluster)) {
  s <- tab %>% filter(cluster == cl, tag != "a_ref_s2")
  if (nrow(s) > 2) cat(sprintf(
    "\n%s: correlation between local effort retained and the shift: r = %.3f\n",
    cl, cor(s$local_effort_kept, s$diff_km)))
}
cat("\nIf the 'Excluding 2020-2021' row sits within the noise while the rows\n")
cat("that remove a large share of a cluster's own coverage do not, the after\n")
cat("period behaves like the baseline and the explanation is coverage, not\n")
cat("the particular years chosen.\n")
