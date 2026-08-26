###############################################################################
#  Methods figure: the affected-zone criterion and the effect distance
#
#  Panel A  two real prediction cells, showing when a location is counted as
#           affected: the posterior of the change at that cell, the 95%
#           credible interval, and the probability of a decrease.
#  Panel B  the northern cluster, showing where the distances are measured
#           from and to.
#
#  Uses the cached joint fit; no refitting.
#
#      source("R/fig_method_schematic.R")
###############################################################################

library(INLA); library(inlabru); library(fmesher)
library(sf); library(terra); library(dplyr); library(ggplot2)

DATA_DIR <- "data"; OUT_DIR <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)

sf_use_s2(FALSE)
main_crs <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")

N_SAMPLES_DIST <- 1000
SEED           <- 20260811L
DIST_CELL_KM   <- 1
THR_SIGNIF     <- 0.975
CLUSTER_NORTH  <- c("West of Sylt", "Butendiek")

## Colours, deliberately few
COL_ZONE  <- "#B2182B"   # area of significant decline
COL_FARM  <- "grey25"    # wind farms
COL_ARROW <- "#2166AC"   # measured distances

as_km <- function(x) {
  y <- if (inherits(x, "sf")) x else st_as_sf(x); st_crs(y) <- main_crs; y
}

D <- new.env(parent = emptyenv())
load(p_data("diver_owf_data.RData"), envir = D)

farms        <- as_km(D$owf_polygons)
pxl          <- as_km(D$prediction_pxl)
pred_outline <- as_km(D$prediction_mask)

if (!all(as.character(st_geometry_type(pxl)) == "POINT")) {
  pxl <- suppressWarnings(st_centroid(pxl))
}
pxl <- suppressWarnings(st_filter(pxl, st_union(st_geometry(pred_outline))))
pxl_xy <- st_coordinates(pxl)

owf <- farms %>% group_by(Name, cluster) %>% summarise(.groups = "drop")
owf_u <- st_sf(geometry = st_union(st_geometry(
  owf %>% filter(cluster %in% CLUSTER_NORTH))))

# ---------------------------------------------------------------- posterior
## Must match MODEL_TAG in joint_likelihood_analysis.R.
MODEL_TAG  <- "m3"
fit_file   <- p_out(paste0("change_fit_", MODEL_TAG, ".rds"))
if (!file.exists(fit_file)) {
  stop("cached fit not found: ", fit_file,
       "\nRun R/joint_likelihood_analysis.R first, or set MODEL_TAG to match.",
       call. = FALSE)
}
change_fit <- readRDS(fit_file)

samples <- as.matrix(generate(change_fit, pxl, ~ log10(exp(spde_change)),
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
  A <- suppressWarnings(st_cast(st_make_valid(
    st_set_crs(st_as_sf(terra::as.polygons(r_neg)), main_crs)), "POLYGON"))
  hit <- lengths(st_intersects(A, owf_one)) > 0
  st_sf(geometry = st_union(st_geometry(A[hit, ])))
}

A_c <- affected_region(zone_raster(THR_SIGNIF), owf_u)

## The extent must match make_dist_raster() in joint_likelihood_analysis.R
## exactly, i.e. the bounding box of the prediction POINTS. Using
## ext(grid_template) instead is wrong: terra::rast(type = "xyz") centres
## cells on the points, so that extent is half a cell larger on every side.
## Half a cell here is 1.130 km in x and 1.040 km in y, which shifts a 1 km
## grid by 0.13 and 0.04 km, changes which cells the boundary crosses, and
## moved this figure to 7.41 km over 257 cells against the reported 7.31 over
## 252.
dist_rast <- terra::distance(
  terra::rast(terra::ext(terra::vect(st_geometry(pxl))),
              resolution = DIST_CELL_KM, crs = main_crs$wkt, vals = 0),
  terra::vect(st_geometry(owf_u))) / 1000

# =============================================================================
# Panel A: two real cells
# =============================================================================
## Two cells that straddle the threshold closely, both near the northern
## cluster. Picking one at 99% and one at 85% would show two arbitrary
## points rather than the boundary, and would leave the reader guessing
## where the cut actually falls. Cells just either side of 97.5% make the
## criterion visible: the two posteriors look almost the same, and the only
## difference is whether the interval clears zero.
near <- which(as.numeric(st_distance(pxl, owf_u)) < 25)
pick <- function(target) near[which.min(abs(p_decrease[near] - target))]
i_yes <- pick(0.985)
i_no  <- pick(0.960)

panel_a_data <- function(i, label) {
  d  <- density(samples[i, ], n = 512)
  ci <- quantile(samples[i, ], c(0.025, 0.975))
  data.frame(x = d$x, y = d$y, panel = label,
             lo = ci[1], hi = ci[2],
             p = mean(samples[i, ] < 0))
}

lab_yes <- sprintf("Counted as affected\nP(decrease) = %.1f%%, interval clears zero",
                   100 * mean(samples[i_yes, ] < 0))
lab_no  <- sprintf("Not counted\nP(decrease) = %.1f%%, interval crosses zero",
                   100 * mean(samples[i_no, ] < 0))

dA <- bind_rows(panel_a_data(i_yes, lab_yes), panel_a_data(i_no, lab_no)) %>%
  mutate(panel = factor(panel, levels = c(lab_yes, lab_no)))
ciA <- dA %>% group_by(panel) %>%
  summarise(lo = first(lo), hi = first(hi), top = max(y), .groups = "drop")

pA <- ggplot(dA, aes(x = x, y = y)) +
  geom_area(data = filter(dA, x < 0), fill = COL_ZONE, alpha = 0.18) +
  geom_line(linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
  geom_segment(data = ciA, inherit.aes = FALSE,
               aes(x = lo, xend = hi, y = -0.08 * top, yend = -0.08 * top),
               linewidth = 0.9) +
  geom_segment(data = ciA, inherit.aes = FALSE,
               aes(x = lo, xend = lo, y = -0.14 * top, yend = -0.02 * top),
               linewidth = 0.9) +
  geom_segment(data = ciA, inherit.aes = FALSE,
               aes(x = hi, xend = hi, y = -0.14 * top, yend = -0.02 * top),
               linewidth = 0.9) +
  geom_text(data = ciA, inherit.aes = FALSE,
            aes(x = (lo + hi) / 2, y = -0.26 * top, label = "95% CI"),
            size = 3.1, vjust = 1) +
  facet_wrap(~ panel, scales = "free") +
  labs(x = expression(paste("change in density, ", log[10], " scale")),
       y = "posterior density", tag = "A",
       subtitle = "a cell is counted as affected when P(decrease) reaches 97.5%") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 10))

# =============================================================================
# Panel B: where the distances are measured
# =============================================================================
win <- st_bbox(st_buffer(A_c, 8))

## The cells of the distance surface crossed by the boundary of the zone
bnd <- terra::vect(st_geometry(st_boundary(A_c)))
ex  <- terra::extract(dist_rast, bnd, cells = TRUE)
ex  <- ex[!is.na(ex[[2]]), ]
pts <- as.data.frame(terra::xyFromCell(dist_rast, ex$cell))
pts$d <- ex[[2]]
pts_sf <- st_as_sf(pts, coords = c("x", "y"), crs = main_crs)

## Three example measurements, spread around the zone
ang <- atan2(pts$y - mean(pts$y), pts$x - mean(pts$x))
sel <- sapply(c(-pi/2, pi/6, 5*pi/6), function(a) which.min(abs(ang - a)))
arrows_sf <- do.call(rbind, lapply(sel, function(k)
  st_sf(d = pts$d[k], geometry = st_nearest_points(pts_sf[k, ], owf_u))))
arrow_lab <- st_as_sf(data.frame(
  lab = sprintf("%.1f km", arrows_sf$d)),
  geometry = st_centroid(st_geometry(arrows_sf)))

## Faint distance contours, for context
dc <- terra::crop(dist_rast, terra::ext(win[c(1, 3, 2, 4)]))
dcl <- st_as_sf(terra::as.contour(dc, levels = c(5, 10, 15, 20)))
st_crs(dcl) <- main_crs

mean_d <- mean(pts$d)

## The figure recomputes the measurement in order to draw it, so it can drift
## away from the reported table if it is run against a different cached fit.
## Refuse to write a figure that would contradict the manuscript.
ref_file <- p_out("effect_distances_summary.csv")
if (file.exists(ref_file)) {
  ref <- read.csv(ref_file)
  ref <- ref[ref$cluster == "north" & grepl("95%", ref$criterion), ]
  if (nrow(ref) == 1) {
    if (abs(mean_d - ref$mean_km) > 0.01 || nrow(pts) != ref$n) {
      stop(sprintf(paste0(
        "figure disagrees with outputs/effect_distances_summary.csv\n",
        "  figure : mean %.3f km, n = %d\n",
        "  table  : mean %.3f km, n = %d\n",
        "Both must come from the same fit. Re-run ",
        "R/joint_likelihood_analysis.R, or check MODEL_TAG above."),
        mean_d, nrow(pts), ref$mean_km, ref$n), call. = FALSE)
    }
    message(sprintf("panel B agrees with the reported table: %.2f km, n = %d",
                    mean_d, nrow(pts)))
  }
}

pB <- ggplot() +
  geom_sf(data = dcl, colour = "grey75", linewidth = 0.3, linetype = 3) +
  geom_sf(data = A_c, fill = COL_ZONE, alpha = 0.12,
          colour = COL_ZONE, linewidth = 0.7) +
  geom_sf(data = owf_u, fill = COL_FARM, colour = COL_FARM) +
  geom_sf(data = pts_sf, size = 0.7, colour = "black") +
  geom_sf(data = arrows_sf, colour = COL_ARROW, linewidth = 0.7) +
  geom_sf_text(data = arrow_lab, aes(label = lab), colour = COL_ARROW,
               size = 3.1, nudge_y = 1.6) +
  coord_sf(xlim = win[c("xmin", "xmax")], ylim = win[c("ymin", "ymax")],
           expand = FALSE, datum = st_crs(4326)) +
  labs(x = NULL, y = NULL, tag = "B",
       subtitle = sprintf(
         "each point on the edge is valued by its distance to the cluster; mean = %.1f km (n = %d)",
         mean_d, nrow(pts))) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_line(colour = "grey92"),
        plot.subtitle = element_text(size = 9.5))

# =============================================================================
ggsave(p_out("figS_method_panelA.png"), pA, width = 7.2, height = 3.0, dpi = 400)
ggsave(p_out("figS_method_panelB.png"), pB, width = 7.2, height = 5.6, dpi = 400)

if (requireNamespace("patchwork", quietly = TRUE)) {
  comb <- patchwork::wrap_plots(pA, pB, ncol = 1, heights = c(1, 1.85))
  ggsave(p_out("figS_method_schematic.png"), comb,
         width = 7.2, height = 8.6, dpi = 400)
  cat("written: outputs/figS_method_schematic.png\n")
} else {
  cat("patchwork not installed; the two panels were saved separately.\n")
}

cat(sprintf("\npanel A cells: p(decrease) = %.3f and %.3f\n",
            p_decrease[i_yes], p_decrease[i_no]))
cat(sprintf("panel B: %d contour cells, mean %.2f km\n", nrow(pts), mean_d))
