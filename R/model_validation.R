###############################################################################
#  Predictive validation of the joint-likelihood model actually fitted here
#
#  The manuscript's validation section rests on Vilela et al. (2021), which
#  validated the same framework at annual resolution on the same survey
#  programme. That remains relevant, but it does not describe the behaviour
#  of the five-year joint-likelihood formulation used in this paper.
#
#  The fit already computes everything needed for that, and discards it:
#
#      control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE)
#
#  CPO is leave-one-out predictive: cpo_i = p(y_i | y_-i). PIT is the
#  corresponding probability integral transform, pit_i = P(Y_i <= y_i | y_-i).
#  Both are per-observation leave-one-out quantities, so this is a
#  cross-validation of THIS model on THESE data, not a reference to earlier
#  work.
#
#  DISCRETENESS. For a continuous response a well calibrated model gives
#  uniform PIT values. For counts it does not: 39.5% of the observations here
#  are zeros, and for an observed zero pit_i = P(Y <= 0) = P(Y = 0), which is
#  large. Zeros therefore pile up near one and the raw PIT is far from
#  uniform even when the model is perfectly adequate. Reporting the raw PIT
#  as a calibration check would misrepresent the fit.
#
#  The randomised PIT (Czado, Gneiting and Held, 2009) removes that artefact
#  and IS uniform under a calibrated model. It needs only what INLA already
#  returns, since P(Y < y) = pit - cpo:
#
#      u_i = (pit_i - cpo_i) + v_i * cpo_i,     v_i ~ Uniform(0, 1)
#
#  Both are reported below; the randomised version is the one to interpret.
#
#  Nothing is refitted. Reads the cached fit and takes seconds.
#
#      source("R/model_validation.R")
###############################################################################

library(INLA); library(inlabru); library(sf); library(dplyr); library(ggplot2)

DATA_DIR <- "data"; OUT_DIR <- "outputs"
p_data <- function(f) file.path(DATA_DIR, f)
p_out  <- function(f) file.path(OUT_DIR,  f)

MODEL_TAG        <- "m3"
BEFORE_PHASE_MAX <- 8
AFTER_PHASE_MIN  <- 17

fit_file <- p_out(paste0("change_fit_", MODEL_TAG, ".rds"))
if (!file.exists(fit_file))
  stop("cached fit not found: ", fit_file, call. = FALSE)
fit <- readRDS(fit_file)

cpo <- fit$cpo$cpo
pit <- fit$cpo$pit
fail <- fit$cpo$failure

## Which observations belong to which period. The joint fit stacks the before
## rows first, then the after rows, in the order the two bru_obs() were given.
sf::sf_use_s2(FALSE)
D <- new.env(parent = emptyenv()); load(p_data("diver_owf_data.RData"), envir = D)
dat <- sf::st_as_sf(D$counts_5km) %>%
  sf::st_drop_geometry() %>%
  filter(!is.na(area), !is.na(NHAT), area > 0)
n_before <- sum(dat$phase <= BEFORE_PHASE_MAX)
n_after  <- sum(dat$phase >= AFTER_PHASE_MIN)
n_obs    <- n_before + n_after

period <- rep(NA_character_, length(cpo))
if (length(cpo) >= n_obs) {
  period[seq_len(n_before)] <- "before"
  period[n_before + seq_len(n_after)] <- "after"
}

cat("\n=================== MODEL FIT ===================\n")
cat(sprintf("DIC   %10.1f   effective parameters %6.1f\n",
            fit$dic$dic,   fit$dic$p.eff))
cat(sprintf("WAIC  %10.1f   effective parameters %6.1f\n",
            fit$waic$waic, fit$waic$p.eff))
cat(sprintf("marginal log-likelihood %10.1f\n", fit$mlik[1]))

cat("\n============ LEAVE-ONE-OUT PREDICTIVE ============\n")
cat(sprintf("observations with CPO           %6d\n", sum(!is.na(cpo))))
cat(sprintf("CPO failures (> 0 = unreliable) %6d  (%.2f%%)\n",
            sum(fail > 0, na.rm = TRUE),
            100 * mean(fail > 0, na.rm = TRUE)))
cat(sprintf("mean log CPO                    %8.3f\n", mean(log(cpo), na.rm = TRUE)))
cat(sprintf("sum -log CPO (smaller is better)%10.1f\n", -sum(log(cpo), na.rm = TRUE)))

## Uniformity of PIT. For discrete responses PIT is not exactly uniform even
## under a correct model, so this is read as "no gross departure" rather than
## as a formal test; the mean and the tail masses are the informative parts.
ok <- !is.na(pit) & !is.na(cpo)
u <- pit[ok]

summarise_pit <- function(x, label) {
  ks <- suppressWarnings(ks.test(x, "punif"))
  cat("\n", label, "\n", sep = "")
  cat(sprintf("  n        %6d\n", length(x)))
  cat(sprintf("  mean     %6.3f   (0.500 expected)\n", mean(x)))
  cat(sprintf("  sd       %6.3f   (0.289 expected)\n", sd(x)))
  cat(sprintf("  P(<0.05) %6.3f   (0.050 expected)\n", mean(x < 0.05)))
  cat(sprintf("  P(>0.95) %6.3f   (0.050 expected)\n", mean(x > 0.95)))
  cat(sprintf("  KS       %6.4f\n", unname(ks$statistic)))
  c(mean = mean(x), sd = sd(x), lo = mean(x < 0.05),
    hi = mean(x > 0.95), ks = unname(ks$statistic))
}

raw <- summarise_pit(u, "Raw PIT - NOT a calibration check for count data")
cat(sprintf("  zeros in the response: %.1f%%, which is what drives the",
            100 * mean(dat$NHAT[dat$phase <= BEFORE_PHASE_MAX |
                                dat$phase >= AFTER_PHASE_MIN] < 0.5)))
cat(" mass near one\n")

## Randomised PIT, averaged over several draws of v so the result does not
## depend on one realisation of the randomisation.
N_RAND <- 20
set.seed(1)
rand <- replicate(N_RAND, {
  v <- runif(sum(ok))
  (pit[ok] - cpo[ok]) + v * cpo[ok]
})
ks_r <- apply(rand, 2, function(x) unname(suppressWarnings(ks.test(x, "punif"))$statistic))
ur <- as.vector(rand)
rpit <- summarise_pit(ur, "Randomised PIT - uniform under a calibrated model")
cat(sprintf("  KS across %d randomisations: %.4f to %.4f\n",
            N_RAND, min(ks_r), max(ks_r)))
ks <- suppressWarnings(ks.test(ur, "punif"))

by_period <- data.frame(period = period, cpo = cpo, pit = pit, fail = fail) %>%
  filter(!is.na(period), !is.na(cpo)) %>%
  group_by(period) %>%
  summarise(n = n(), mean_log_cpo = mean(log(cpo)),
            pit_mean = mean(pit, na.rm = TRUE),
            pit_lo = mean(pit < 0.05, na.rm = TRUE),
            pit_hi = mean(pit > 0.95, na.rm = TRUE),
            cpo_failures = sum(fail > 0), .groups = "drop")
cat("\nBy period\n")
print(as.data.frame(by_period), digits = 3, row.names = FALSE)

out <- data.frame(
  dic = fit$dic$dic, dic_p_eff = fit$dic$p.eff,
  waic = fit$waic$waic, waic_p_eff = fit$waic$p.eff,
  mlik = fit$mlik[1],
  n_cpo = sum(!is.na(cpo)), cpo_failures = sum(fail > 0, na.rm = TRUE),
  mean_log_cpo = mean(log(cpo), na.rm = TRUE),
  pit_raw_mean = raw[["mean"]], pit_raw_ks = raw[["ks"]],
  rpit_mean = rpit[["mean"]], rpit_sd = rpit[["sd"]],
  rpit_below_05 = rpit[["lo"]], rpit_above_95 = rpit[["hi"]],
  rpit_ks = rpit[["ks"]],
  rpit_ks_min = min(ks_r), rpit_ks_max = max(ks_r))
write.csv(out, p_out("model_validation.csv"), row.names = FALSE)
write.csv(by_period, p_out("model_validation_by_period.csv"), row.names = FALSE)

pd <- rbind(
  data.frame(pit = u,  kind = "Raw PIT (discreteness artefact)"),
  data.frame(pit = rand[, 1], kind = "Randomised PIT"))
p <- ggplot(pd, aes(x = pit)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#6A51A3",
                 colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = length(u) / 20, linetype = 2, colour = "grey30") +
  facet_wrap(~ kind) +
  labs(x = "PIT", y = "Observations",
       caption = paste("Dashed line: the count expected under a uniform",
                       "distribution. With 39.5% zeros the raw PIT cannot be",
                       "uniform;\nthe randomised version is the calibration",
                       "check.")) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA))
ggsave(p_out("figS_pit_histogram.png"), p, width = 8, height = 4, dpi = 300)

cat("\nwritten: outputs/model_validation.csv, model_validation_by_period.csv,",
    "figS_pit_histogram.png\n")
