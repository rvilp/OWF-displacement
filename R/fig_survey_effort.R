###############################################################################
#  Survey coverage by year
#
#  Small multiples in the style of Vilela et al. (2021), Fig. 2, showing where
#  the survey effort fell in each spring season.
#
#  NOTE ON WHAT IS PLOTTED. The archived dataset holds counts aggregated to
#  the nodes of the 5 km mesh, not the raw transect geometry, so the flight
#  lines themselves cannot be redrawn from it. Each point is therefore a mesh
#  node that was surveyed that year, sized by the effort recorded at it. This
#  is what the model actually sees, which for an appendix on the sensitivity
#  of the results to survey coverage is the more relevant depiction.
#
#  Only the eleven seasons that enter the analysis are shown: 2001-2005 and
#  2008 for the before period, 2017-2021 for the after period. The seasons
#  between 2009 and 2016 fall in the construction phase, belong to neither
#  period and are not plotted. Panel titles carry the period each year belongs
#  to, so that the gaps the sensitivity analysis probes are visible directly.
#
#      source("R/fig_survey_effort.R")
###############################################################################

library(sf); library(dplyr); library(ggplot2)

DATA_DIR <- "data"; OUT_DIR <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)

sf_use_s2(FALSE)
main_crs <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")

BEFORE_PHASE_MAX <- 8
AFTER_PHASE_MIN  <- 17
MAP_MARGIN_KM    <- 8

as_km <- function(x) { y <- if (inherits(x,"sf")) x else st_as_sf(x); st_crs(y) <- main_crs; y }
D <- new.env(parent = emptyenv()); load(p_data("diver_owf_data.RData"), envir = D)

countdata    <- as_km(D$counts_5km)
farms        <- as_km(D$owf_polygons)
pxl          <- as_km(D$prediction_pxl)
pred_outline <- as_km(D$prediction_mask)
conc_area    <- as_km(D$hd_mask)
spa_area     <- as_km(D$spa_mask)

owf <- farms %>% group_by(Name, cluster) %>% summarise(.groups = "drop")

## Survey technique collapses to two classes for display: the distinction
## that matters visually is visual against digital, not which digital system.
d <- countdata %>%
  filter(!is.na(area), area > 0,
         phase <= BEFORE_PHASE_MAX | phase >= AFTER_PHASE_MIN) %>%
  mutate(year = 2000L + as.integer(phase),
         type = ifelse(as.character(method) == "conventional",
                       "conventional", "digital"),
         period = ifelse(phase <= BEFORE_PHASE_MAX, "before", "after"))

## One panel per year, labelled with the period so the design is legible
## without cross-referencing the text.
lab <- d %>% distinct(year, period) %>% arrange(year) %>%
  mutate(panel = sprintf("%d  (%s)", year, period))
d <- d %>% left_join(lab, by = c("year", "period")) %>%
  mutate(panel = factor(panel, levels = lab$panel))

cat("years: ", nrow(lab), " | panels: ", paste(levels(d$panel), collapse = " | "),
    "\n", sep = "")

bb <- st_bbox(pxl)
xlim <- c(bb[["xmin"]] - MAP_MARGIN_KM, bb[["xmax"]] + MAP_MARGIN_KM)
ylim <- c(bb[["ymin"]] - MAP_MARGIN_KM, bb[["ymax"]] + MAP_MARGIN_KM)

## Coastline if download_emodnet_covariates.R has been run, otherwise omitted
## rather than substituted by something that is not a coastline.
## The coastline is an optional decoration and lives with the other
## environmental layers outside the repository, so its absence is not an
## error: the figure is simply drawn without it.
coast <- NULL
cf <- try(local({ source("R/covariate_paths.R", local = TRUE)
                  COV_FILES[["coastline"]] }), silent = TRUE)
if (!inherits(cf, "try-error") && file.exists(cf)) {
  coast <- st_read(cf, quiet = TRUE); cat("coastline included\n")
} else {
  cat("no coastline layer; drawing without it\n")
}

p <- ggplot() +
  { if (!is.null(coast))
      geom_sf(data = coast, colour = "grey45", linewidth = 0.25) } +
  geom_sf(data = pred_outline, fill = NA, colour = "grey55", linewidth = 0.25) +
  geom_sf(data = conc_area,    fill = NA, colour = "grey35", linewidth = 0.35) +
  geom_sf(data = spa_area,     fill = NA, colour = "#1B7837", linewidth = 0.4) +
  geom_sf(data = d, aes(colour = type, size = area), alpha = 0.55, stroke = 0) +
  geom_sf(data = owf, fill = NA, colour = "#B2182B", linewidth = 0.3) +
  facet_wrap(~ panel, ncol = 4) +
  scale_colour_manual(values = c(conventional = "#1B9E77", digital = "#3F007D"),
                      name = "Survey type") +
  scale_size_area(max_size = 2.2, name = expression(paste("Effort [km"^2, "]")),
                  breaks = c(25, 100, 250)) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, datum = st_crs(4326)) +
  labs(x = "Longitude", y = "Latitude") +
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_line(colour = "grey93", linewidth = 0.2),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8, margin = margin(2, 2, 2, 2)),
        axis.text = element_text(size = 6),
        legend.position = "bottom", legend.box = "horizontal")

ggsave(p_out("figS_survey_effort.png"), p, width = 8.5, height = 7.2, dpi = 300)
cat("written: outputs/figS_survey_effort.png\n")

## The same information without the map, for the years actually used: effort
## by year and latitude band. Kept as a compact companion panel.
band <- function(y) ifelse(y >= 6050, "north", ifelse(y < 5990, "south", "middle"))
d$band <- band(st_coordinates(d)[, 2])
eff <- d %>% st_drop_geometry() %>%
  group_by(year, period, band) %>% summarise(effort = sum(area), .groups = "drop") %>%
  mutate(band = factor(band, levels = c("north", "middle", "south")))

## The gap between 2008 and 2017 is real and is left visible rather than
## closed up, because the two periods are separated by the construction phase.
p2 <- ggplot(eff, aes(x = factor(year), y = effort, fill = band)) +
  geom_col(width = 0.75) +
  facet_wrap(~ band, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(north = "#2166AC", middle = "grey55",
                               south = "#B2182B"), guide = "none") +
  labs(x = NULL, y = expression(paste("Survey effort [km"^2, "]"))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA))

ggsave(p_out("figS_survey_effort_by_band.png"), p2, width = 7, height = 5.5, dpi = 300)
cat("written: outputs/figS_survey_effort_by_band.png\n")
