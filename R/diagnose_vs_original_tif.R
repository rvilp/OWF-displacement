###############################################################################
#  Definitive comparison against the 2022 result
#
#  pred_change_975.tif holds q0.975 of the change field from the 2022 fit, on
#  a 150 x 150 grid at 2.2595 x 2.0806 km with 6089 cells carrying data - the
#  same grid as prediction_pxl.rds, cell for cell. There is therefore no
#  rasterisation step between the two pipelines.
#
#  This runs the CURRENT measurement code on the ORIGINAL field. If it
#  returns the published 7.59 km / N = 256 / 1170.6 km2, then the measurement
#  code is equivalent and the whole difference comes from the field itself,
#  i.e. from refitting. If it returns something else, part of the difference
#  is in the measurement.
#
#  Put pred_change_975.tif in data/ and run:
#      source("R/diagnose_vs_original_tif.R")
###############################################################################

library(sf); library(terra); library(dplyr)

DATA_DIR <- "data"; OUT_DIR <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)

sf_use_s2(FALSE)
main_crs <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")

TIF          <- "pred_change_975.tif"
DIST_CELL_KM <- 1
CLUSTER_NORTH <- c("West of Sylt", "Butendiek")
CLUSTER_SOUTH <- c("UMBO")

if (!file.exists(p_data(TIF))) {
  stop("Put ", TIF, " in data/ first.", call. = FALSE)
}

as_km <- function(x) {
  y <- if (inherits(x, "sf")) x else st_as_sf(x); st_crs(y) <- main_crs; y
}
D <- new.env(parent = emptyenv())
load(p_data("diver_owf_data.RData"), envir = D)

farms <- as_km(D$owf_polygons)
owf   <- farms %>% group_by(Name, cluster) %>% summarise(.groups = "drop")
owf_north_u <- st_sf(geometry = st_union(st_geometry(
  owf %>% filter(cluster %in% CLUSTER_NORTH))))
owf_south_u <- st_sf(geometry = st_union(st_geometry(
  owf %>% filter(cluster %in% CLUSTER_SOUTH))))

## The original field. Values are q0.975 of the log10 change; the affected
## zone is where that is below zero, exactly as in
## impact_distance_and_habitat_loss.R (subset(decr, pred_change_975 <= 0)).
r <- terra::rast(p_data(TIF))
terra::crs(r) <- main_crs$wkt
cell_area <- prod(terra::res(r))

cat("\ngrid            :", terra::ncol(r), "x", terra::nrow(r), "\n")
cat("resolution      :", paste(round(terra::res(r), 4), collapse = " x "), "km\n")
cat("cells with data :", sum(!is.na(terra::values(r))), "\n")
cat("cell area       :", round(cell_area, 4), "km2\n")

# ---------------------------------------------------------------- same code
affected_region <- function(r, owf_u) {
  r_neg <- terra::classify(r, matrix(c(-Inf, 0, 1, 0, Inf, NA),
                                     ncol = 3, byrow = TRUE))
  if (all(is.na(terra::values(r_neg)))) return(NULL)
  A <- suppressWarnings(st_cast(st_make_valid(
    st_set_crs(st_as_sf(terra::as.polygons(r_neg)), st_crs(owf_u))), "POLYGON"))
  hit <- lengths(st_intersects(A, owf_u)) > 0
  if (!any(hit)) return(NULL)
  st_sf(geometry = st_union(st_geometry(A[hit, ])))
}

make_dist_raster <- function(cluster_owf, template) {
  d <- terra::rast(terra::ext(template), resolution = DIST_CELL_KM,
                   crs = main_crs$wkt, vals = 0)
  terra::distance(d, terra::vect(st_geometry(cluster_owf))) / 1000
}

measure <- function(owf_u, label) {
  A_c <- affected_region(r, owf_u)
  if (is.null(A_c)) return(NULL)
  dr <- make_dist_raster(owf_u, r)
  d  <- as.numeric(stats::na.omit(
    terra::extract(dr, terra::vect(st_geometry(st_boundary(A_c))))[[2]]))
  data.frame(cluster = label,
             n_contour_cells = length(d),
             area_km2 = as.numeric(sum(st_area(A_c))),
             mean_km  = mean(d),
             median_km = median(d),
             max_km   = max(d))
}

res <- bind_rows(measure(owf_north_u, "north"), measure(owf_south_u, "south"))

cat("\n=== CURRENT CODE ON THE 2022 FIELD ===\n")
print(as.data.frame(res), digits = 4, row.names = FALSE)

cat("\n=== PUBLISHED (impact_distance_and_habitat_loss.R) ===\n")
cat("north  N = 256  area = 1170.6 km2  mean = 7.59 km\n")
cat("south  N = 176  area =  879.1 km2  mean = 6.62 km\n")

cat("\n=== CURRENT CODE ON THE CURRENT FIT ===\n")
cat("north  N = 238  area = 1138   km2  mean = 7.15 km\n")
cat("south  N = 164  area =  860   km2  mean = 6.20 km\n")

cat("\nIf the first block matches the published block, the measurement code is\n")
cat("equivalent and the difference is entirely in the fitted field.\n\n")

write.csv(res, p_out("comparison_vs_original_tif.csv"), row.names = FALSE)
