###############################################################################
#  Varying both periods at once, and measuring WHERE the zone is rather than
#  only how far it reaches
#
#  The one-period-at-a-time tests established two things.
#
#    Baseline: all variants but one lower the estimate, and the size of the
#      shift tracks the survey effort removed near the cluster concerned. The
#      exception is the southern cluster with 2003 excluded, which rises by
#      0.50 km. Removing whole years is systematically GENTLER than removing
#      the same effort at random, so what drives it is how much evidence
#      remains, not which years are in.
#
#    After period: not monotone. Removing 2020-2021, which is 0.5% of the
#      northern and 2.4% of the southern local effort, moves the estimate by
#      0.2 km. But removing 2018 raises the southern estimate by 2.73 km
#      since without 2018 only a peripheral patch still clears the threshold. 
#      Reducing the data can therefore move the estimate in either direction.
#
#  So the mean distance responds to how much data remain in each period, and
#  in the after period it does so non-monotonically. A zone of decline adjacent 
#  to both clusters was detected in every one of the sixteen variants. 
#  This script measures that directly.
#  Two additions over the previous scripts:
#
#    1. Both periods vary together, in a 3 x 3 factorial.
#
#    2. For each cell of the factorial the affected zone is compared with the
#       reference zone by JACCARD OVERLAP, |A n B| / |A u B|. That answers
#       "is the effect in the same place?" rather than "is the number the
#       same?". Also reported: what fraction of the reference zone is
#       recovered, and the distance between the two zone centroids.
#
#  Nine fits of a few seconds each, cached under outputs/sens_both_*.rds.
#
#      source("R/sensitivity_both_periods.R")
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
N_SAMPLES_DIST <- 10000
SEED <- 20260811L
DIST_CELL_KM <- 1; THR_SIGNIF <- 0.975; EDGE_TOL_KM <- 2 * DIST_CELL_KM
CLUSTER_NORTH <- c("West of Sylt", "Butendiek"); CLUSTER_SOUTH <- c("UMBO")

## Three levels per period, chosen to span what the single-period runs showed:
## the full data, the single most damaging exclusion, and the harshest split.
BEFORE_LEVELS <- list(
  full    = list(lab = "2001-2008 (all)", keep = c(1,2,3,4,5,8)),
  no2008  = list(lab = "no 2008",         keep = c(1,2,3,4,5)),
  late    = list(lab = "2005-2008 only",  keep = c(5,8))
)
AFTER_LEVELS <- list(
  full    = list(lab = "2017-2021 (all)", keep = 17:21),
  no2018  = list(lab = "no 2018",         keep = c(17,19,20,21)),
  last3   = list(lab = "2019-2021 only",  keep = c(19,20,21))
)

# ---------------------------------------------------------------- data
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

## Zone comparison. Areas are computed by counting prediction cells rather
## than with st_area, so that they are on the same footing as everything
## else reported and are not affected by polygon simplification.
cells_in <- function(A) if (is.null(A)) integer(0) else
  which(lengths(st_intersects(pxl, A)) > 0)

zone_stats <- function(A, cl) {
  if (is.null(A)) return(data.frame(
    cluster = cl, detected = FALSE, n = NA_integer_, mean_km = NA_real_,
    prop_on_edge = NA_real_, area_km2 = 0))
  b <- terra::vect(st_geometry(st_boundary(A)))
  d <- terra::extract(dist_rast[[cl]], b)[[2]]
  e <- terra::extract(dist_to_edge,    b)[[2]]
  ok <- !is.na(d)
  data.frame(cluster = cl, detected = TRUE, n = sum(ok),
             mean_km = if (any(ok)) mean(d[ok]) else NA_real_,
             prop_on_edge = if (any(ok)) mean(as.numeric(e[ok]) <= EDGE_TOL_KM) else NA_real_,
             area_km2 = length(cells_in(A)) * cell_area)
}

matern <- inla.spde2.pcmatern(spatial_mesh, prior.sigma = PRIOR_SIGMA,
                              prior.range = PRIOR_RANGE)
cmp <- ~ spde_before(geometry, model = matern) +
         spde_change(geometry, model = matern) + Intercept(1) + Intercept_after(1)
opts <- list(control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
             control.inla = list(int.strategy = "ccd"),
             control.predictor = list(compute = TRUE,
               quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)))

fit_cell <- function(bk, ak, tag) {
  bef <- dat %>% filter(phase %in% bk)
  aft <- dat %>% filter(phase %in% ak)
  f <- p_out(paste0("sens_both_", tag, ".rds"))
  if (file.exists(f)) { fit <- readRDS(f) } else {
    t0 <- Sys.time()
    fit <- bru(cmp,
      bru_obs("nbinomial", data = bef, formula = NHAT ~ spde_before + Intercept,
              E = bef$E_eff),
      bru_obs("nbinomial", data = aft,
              formula = NHAT ~ spde_before + spde_change + Intercept_after,
              E = aft$E_eff),
      options = c(opts, list(bru_max_iter = 1)))
    saveRDS(fit, f)
    cat(sprintf("    fitted in %.0f s\n", difftime(Sys.time(), t0, units = "secs")))
  }
  s <- as.matrix(generate(fit, pxl, ~ log10(exp(spde_change)),
                          n.samples = N_SAMPLES_DIST, seed = SEED + 3L))
  list(p_dec = rowMeans(s < 0),
       before_effort = sum(bef$area), after_effort = sum(aft$area))
}

# ---------------------------------------------------------------- reference
cat("reference cell (all data)\n")
r0 <- fit_cell(BEFORE_LEVELS$full$keep, AFTER_LEVELS$full$keep, "full_full")
ref_zone  <- lapply(names(owf_u),
                    function(cl) affected_region(draw_raster(THR_SIGNIF - r0$p_dec), owf_u[[cl]]))
names(ref_zone) <- names(owf_u)
ref_cells <- lapply(ref_zone, cells_in)
ref_stats <- bind_rows(lapply(names(owf_u), function(cl) zone_stats(ref_zone[[cl]], cl)))
print(as.data.frame(ref_stats), digits = 4, row.names = FALSE)

# ---------------------------------------------------------------- factorial
res <- list()
for (bn in names(BEFORE_LEVELS)) for (an in names(AFTER_LEVELS)) {
  tag <- paste0(bn, "_", an)
  cat(sprintf("\nbefore = %-16s after = %-16s\n",
              BEFORE_LEVELS[[bn]]$lab, AFTER_LEVELS[[an]]$lab))
  fc <- fit_cell(BEFORE_LEVELS[[bn]]$keep, AFTER_LEVELS[[an]]$keep, tag)
  rr <- draw_raster(THR_SIGNIF - fc$p_dec)

  for (cl in names(owf_u)) {
    A <- affected_region(rr, owf_u[[cl]])
    st <- zone_stats(A, cl)
    cc <- cells_in(A); rc <- ref_cells[[cl]]
    inter <- length(intersect(cc, rc)); uni <- length(union(cc, rc))

    st$before <- BEFORE_LEVELS[[bn]]$lab
    st$after  <- AFTER_LEVELS[[an]]$lab
    st$tag <- tag
    st$before_effort_km2 <- fc$before_effort
    st$after_effort_km2  <- fc$after_effort
    st$jaccard        <- if (uni) inter / uni else NA_real_
    st$prop_of_ref    <- if (length(rc)) inter / length(rc) else NA_real_
    st$centroid_shift_km <- if (is.null(A) || is.null(ref_zone[[cl]])) NA_real_ else
      as.numeric(st_distance(st_centroid(st_union(A)),
                             st_centroid(st_union(ref_zone[[cl]]))))
    res[[length(res) + 1]] <- st
    cat(sprintf("    %-6s detected=%s  mean %5.2f km  area %6.0f  jaccard %.2f  of_ref %.2f  shift %.1f km\n",
                cl, st$detected, st$mean_km, st$area_km2,
                st$jaccard, st$prop_of_ref, st$centroid_shift_km))
  }
}

out <- bind_rows(res) %>%
  select(before, after, cluster, detected, n, mean_km, area_km2,
         jaccard, prop_of_ref, centroid_shift_km, prop_on_edge,
         before_effort_km2, after_effort_km2, tag)
write.csv(out, p_out("sensitivity_both_periods.csv"), row.names = FALSE)

cat("\n\n================= FACTORIAL SUMMARY =================\n")
print(as.data.frame(out %>% select(before, after, cluster, detected, mean_km,
                                   area_km2, jaccard, prop_of_ref,
                                   centroid_shift_km)),
      digits = 3, row.names = FALSE)

cat("\n--- the claim being tested ---\n")
cat("detected      : is a zone of decline found next to this cluster at all?\n")
cat("jaccard       : overlap with the reference zone, |A n B| / |A u B|\n")
cat("prop_of_ref   : share of the reference zone recovered\n")
cat("centroid_shift: how far the zone moved\n\n")
for (cl in unique(out$cluster)) {
  s <- out %>% filter(cluster == cl)
  cat(sprintf("%s: detected in %d of %d cells | mean distance %.1f-%.1f km | jaccard %.2f-%.2f | centroid shift %.1f-%.1f km\n",
              cl, sum(s$detected), nrow(s),
              min(s$mean_km, na.rm = TRUE), max(s$mean_km, na.rm = TRUE),
              min(s$jaccard, na.rm = TRUE), max(s$jaccard, na.rm = TRUE),
              min(s$centroid_shift_km, na.rm = TRUE),
              max(s$centroid_shift_km, na.rm = TRUE)))
}
cat("\nIf the zone is always detected and the centroid barely moves while the\n")
cat("mean distance ranges widely, then the location of the effect is robust\n")
cat("and its measured extent is not - which is what the manuscript should say.\n")

# ---------------------------------------------------------------- figure
p <- ggplot(out %>% filter(detected),
            aes(x = after, y = before, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f\n%.1f km", jaccard, mean_km)),
            size = 3.1, colour = "grey15") +
  facet_wrap(~ cluster) +
  scale_fill_gradient(low = "#F7F7F7", high = "#3F007D", limits = c(0, 1),
                      name = "Jaccard overlap\nwith reference zone") +
  labs(x = "After period", y = "Baseline period",
       caption = "Cell label: Jaccard overlap above, mean effect distance below.") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "right")

ggsave(p_out("figS_both_periods_sensitivity.png"), p,
       width = 8.5, height = 4.5, dpi = 300)
cat("\nwritten: outputs/sensitivity_both_periods.csv,",
    "figS_both_periods_sensitivity.png\n")
