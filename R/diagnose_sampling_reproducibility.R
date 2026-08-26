###############################################################################
#  Is the posterior sampling actually reproducible?
#
#  Two runs of joint_likelihood_analysis.R against the SAME cached fit and the
#  SAME package versions produced different numbers:
#
#      16 Aug   north  n = 252   mean = 7.3115 km   area = 1137.67 km2
#      18 Aug   north  n = 255   mean = 7.4498 km   area = 1147.07 km2
#
#  The fits were not refitted (outputs/*_m3.rds untouched), so the difference
#  can only come from generate() / predict(), which draw from the posterior.
#
#  INLA's inla.posterior.sample() is only reproducible from `seed` when it
#  runs SINGLE-THREADED. inlabru is documented to force num.threads = "1:1:1"
#  when seed != 0, but that behaviour has not been verified for the version
#  actually in use here.
#
#  This script draws twice with the same seed inside one session and compares.
#  It costs a few minutes and refits nothing.
#
#      source("R/diagnose_sampling_reproducibility.R")
###############################################################################

library(INLA); library(inlabru); library(fmesher)
library(sf); library(terra); library(dplyr)

DATA_DIR <- "data"; OUT_DIR <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)

sf_use_s2(FALSE)
main_crs <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")

MODEL_TAG      <- "m3"
N_SAMPLES_DIST <- 1000
SEED           <- 20260811L

cat("inlabru ", as.character(packageVersion("inlabru")),
    " | INLA ", as.character(packageVersion("INLA")), "\n", sep = "")
cat("inla.getOption('num.threads'): ", INLA::inla.getOption("num.threads"), "\n\n")

as_km <- function(x) {
  y <- if (inherits(x, "sf")) x else st_as_sf(x); st_crs(y) <- main_crs; y
}
D <- new.env(parent = emptyenv())
load(p_data("diver_owf_data.RData"), envir = D)

pxl          <- as_km(D$prediction_pxl)
pred_outline <- as_km(D$prediction_mask)
if (!all(as.character(st_geometry_type(pxl)) == "POINT")) {
  pxl <- suppressWarnings(st_centroid(pxl))
}
pxl <- suppressWarnings(st_filter(pxl, st_union(st_geometry(pred_outline))))
cat("prediction grid: ", nrow(pxl), " cells\n\n", sep = "")

fit_file <- p_out(paste0("change_fit_", MODEL_TAG, ".rds"))
if (!file.exists(fit_file)) stop("no cached fit at ", fit_file, call. = FALSE)
change_fit <- readRDS(fit_file)

draw <- function(label, ...) {
  t0 <- Sys.time()
  s <- as.matrix(generate(change_fit, pxl, ~ log10(exp(spde_change)),
                          n.samples = N_SAMPLES_DIST, seed = SEED + 3L, ...))
  cat(sprintf("%-28s %5.1f s\n", label,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  s
}

compare <- function(a, b, label) {
  pa <- rowMeans(a < 0); pb <- rowMeans(b < 0)
  same_draws <- identical(a, b)
  n975 <- c(sum(pa >= 0.975), sum(pb >= 0.975))
  cat("\n--- ", label, " ---\n", sep = "")
  cat("  identical sample matrices : ", same_draws, "\n")
  cat("  max |P(dec) difference|   : ", sprintf("%.5f", max(abs(pa - pb))), "\n")
  cat("  cells with P(dec) >= 0.975: ", n975[1], " vs ", n975[2],
      "  (difference ", n975[1] - n975[2], ")\n", sep = "")
  cat("  cells changing side of the threshold: ",
      sum((pa >= 0.975) != (pb >= 0.975)), "\n")
  invisible(same_draws)
}

# ---------------------------------------------------------------- as the script does it
cat("A. exactly as joint_likelihood_analysis.R calls it\n")
s1 <- draw("  draw 1")
s2 <- draw("  draw 2")
ok_default <- compare(s1, s2, "same seed, default threading")

# ---------------------------------------------------------------- forced single thread
cat("\n\nB. with num.threads forced to a single thread\n")
s3 <- draw("  draw 3", num.threads = "1:1")
s4 <- draw("  draw 4", num.threads = "1:1")
ok_single <- compare(s3, s4, "same seed, num.threads = '1:1'")

# ---------------------------------------------------------------- verdict
cat("\n\n================ VERDICT ================\n")
if (ok_default) {
  cat("Sampling is already reproducible as called. The 16 vs 18 August\n")
  cat("difference has another cause - look at the prediction grid or at\n")
  cat("anything edited in Sections 1-2 between the two runs.\n")
} else if (ok_single) {
  cat("Sampling is NOT reproducible as currently called, but IS reproducible\n")
  cat("with num.threads = '1:1'. Add num.threads = '1:1' to every predict()\n")
  cat("and generate() call in joint_likelihood_analysis.R, re-run, and treat\n")
  cat("the resulting numbers as the ones to report.\n")
} else {
  cat("Sampling is not reproducible even single-threaded. Do not chase the\n")
  cat("seed: instead report the Monte Carlo uncertainty by drawing the whole\n")
  cat("measurement several times and quoting the spread, and raise\n")
  cat("N_SAMPLES_DIST until that spread is smaller than the last digit\n")
  cat("reported.\n")
}
cat("=========================================\n")
