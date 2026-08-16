###############################################################################
#  terra::distance() vs distanceto::distance_raster()
#
#  Running the current measurement code on the 2022 change field reproduced
#  the affected zone exactly (N = 256 / 176, area 1170.577 / 879.108 km2) but
#  gave a mean effect distance of 7.47 / 6.39 km against the published
#  7.59 / 6.62. With the zone and the sampled cells identical, the only
#  remaining difference is the distance surface itself.
#
#  This compares the two implementations directly:
#    - grid geometry (extent, resolution, origin)
#    - values at the same coordinates
#    - value inside a wind farm polygon (a definitional check: distance to
#      the polygon as a filled area is 0 inside, distance to its boundary is
#      not)
#    - the resulting mean effect distance over the identical 2022 zone
#
#      source("R/diagnose_distance_surfaces.R")
###############################################################################

library(sf); library(terra); library(dplyr)

DATA_DIR <- "data"; OUT_DIR <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)

sf_use_s2(FALSE)
main_crs <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")

CLUSTER_NORTH <- c("West of Sylt", "Butendiek")
CLUSTER_SOUTH <- c("UMBO")
DIST_CELL_KM  <- 1

has_dt <- requireNamespace("distanceto", quietly = TRUE)
if (!has_dt) {
  cat("\ndistanceto is not installed. Install it with\n",
      "  install.packages('distanceto')\n",
      "or, if it has been archived on CRAN,\n",
      "  remotes::install_github('mdsumner/distanceto')\n\n",
      "Without it this script can only report the terra side.\n\n")
}

as_km <- function(x) {
  y <- if (inherits(x, "sf")) x else st_as_sf(x); st_crs(y) <- main_crs; y
}
D <- new.env(parent = emptyenv())
load(p_data("diver_owf_data.RData"), envir = D)

farms <- as_km(D$owf_polygons)
owf   <- farms %>% group_by(Name, cluster) %>% summarise(.groups = "drop")
owf_u <- list(
  north = st_sf(geometry = st_union(st_geometry(
    owf %>% filter(cluster %in% CLUSTER_NORTH)))),
  south = st_sf(geometry = st_union(st_geometry(
    owf %>% filter(cluster %in% CLUSTER_SOUTH))))
)

r <- terra::rast(p_data("pred_change_975.tif"))
terra::crs(r) <- main_crs$wkt
bb <- st_bbox(c(xmin = terra::xmin(r), xmax = terra::xmax(r),
                ymin = terra::ymin(r), ymax = terra::ymax(r)), crs = main_crs)

## The 2022 affected zone, reproduced exactly by the current code.
affected_region <- function(r, owf_one) {
  r_neg <- terra::classify(r, matrix(c(-Inf, 0, 1, 0, Inf, NA),
                                     ncol = 3, byrow = TRUE))
  A <- suppressWarnings(st_cast(st_make_valid(
    st_set_crs(st_as_sf(terra::as.polygons(r_neg)), main_crs)), "POLYGON"))
  hit <- lengths(st_intersects(A, owf_one)) > 0
  if (!any(hit)) return(NULL)
  st_sf(geometry = st_union(st_geometry(A[hit, ])))
}

## --- the two surfaces -------------------------------------------------------
surface_terra <- function(owf_one) {
  d <- terra::rast(terra::ext(r), resolution = DIST_CELL_KM,
                   crs = main_crs$wkt, vals = 0)
  terra::distance(d, terra::vect(st_geometry(owf_one))) / 1000
}

surface_dt <- function(owf_one) {
  distanceto::distance_raster(owf_one, cellsize = DIST_CELL_KM, extent = bb)
}

out <- list()

for (cl in names(owf_u)) {
  cat("\n=====================================================\n")
  cat("CLUSTER:", cl, "\n")
  cat("=====================================================\n")

  A_c <- affected_region(r, owf_u[[cl]])
  bnd <- terra::vect(st_geometry(st_boundary(A_c)))

  st <- surface_terra(owf_u[[cl]])
  d_terra <- as.numeric(stats::na.omit(terra::extract(st, bnd)[[2]]))

  cat("\nterra::distance\n")
  cat("  grid       :", terra::ncol(st), "x", terra::nrow(st),
      " res", paste(round(terra::res(st), 4), collapse = " x "), "\n")
  cat("  contour n  :", length(d_terra), "\n")
  cat("  mean       : ", sprintf("%.4f km\n", mean(d_terra)))
  cat("  inside OWF : ", sprintf("%.4f\n", mean(terra::extract(
        st, terra::vect(st_geometry(owf_u[[cl]])))[[2]], na.rm = TRUE)),
      "  (0 means the polygon is treated as a filled area)\n")

  row <- data.frame(cluster = cl, surface = "terra::distance",
                    n = length(d_terra), mean_km = mean(d_terra),
                    median_km = median(d_terra), max_km = max(d_terra))

  if (has_dt) {
    sd_ <- try(surface_dt(owf_u[[cl]]), silent = TRUE)
    if (inherits(sd_, "try-error")) {
      cat("\ndistanceto::distance_raster failed:\n  ",
          conditionMessage(attr(sd_, "condition")), "\n")
    } else {
      sd_r <- terra::rast(sd_)
      terra::crs(sd_r) <- main_crs$wkt
      d_dt <- as.numeric(stats::na.omit(terra::extract(sd_r, bnd)[[2]]))

      cat("\ndistanceto::distance_raster\n")
      cat("  grid       :", terra::ncol(sd_r), "x", terra::nrow(sd_r),
          " res", paste(round(terra::res(sd_r), 4), collapse = " x "), "\n")
      cat("  contour n  :", length(d_dt), "\n")
      cat("  mean       : ", sprintf("%.4f km\n", mean(d_dt)))
      cat("  inside OWF : ", sprintf("%.4f\n", mean(terra::extract(
            sd_r, terra::vect(st_geometry(owf_u[[cl]])))[[2]], na.rm = TRUE)))

      ## Same locations, both surfaces
      pts <- terra::vect(st_as_sf(
        as.data.frame(terra::crds(st, na.rm = TRUE)[
          sample(seq_len(sum(!is.na(terra::values(st)))),
                 min(5000, sum(!is.na(terra::values(st))))), , drop = FALSE]),
        coords = c("x", "y"), crs = main_crs))
      a <- terra::extract(st,   pts)[[2]]
      b <- terra::extract(sd_r, pts)[[2]]
      ok <- !is.na(a) & !is.na(b)
      dif <- b[ok] - a[ok]

      cat("\ncell-by-cell (", sum(ok), " random locations)\n", sep = "")
      cat("  distanceto - terra :\n")
      cat("    mean  ", sprintf("%+.4f km\n", mean(dif)))
      cat("    median", sprintf("%+.4f km\n", median(dif)))
      cat("    range ", sprintf("%+.4f to %+.4f km\n", min(dif), max(dif)))
      cat("    ratio (median of b/a, a>1km): ",
          sprintf("%.5f\n", median((b[ok] / a[ok])[a[ok] > 1])))
      cat("    correlation: ", sprintf("%.6f\n", cor(a[ok], b[ok])))

      row <- rbind(row, data.frame(
        cluster = cl, surface = "distanceto::distance_raster",
        n = length(d_dt), mean_km = mean(d_dt),
        median_km = median(d_dt), max_km = max(d_dt)))
    }
  }
  out[[cl]] <- row
}

res <- bind_rows(out)
cat("\n\n=== MEAN EFFECT DISTANCE OVER THE IDENTICAL 2022 ZONE ===\n")
print(as.data.frame(res), digits = 5, row.names = FALSE)
cat("\npublished: north 7.59 km, south 6.62 km\n\n")

write.csv(res, p_out("distance_surface_comparison.csv"), row.names = FALSE)
