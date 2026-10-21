# -------------------------------------------------------------------------
# 2017 Tehuantepec earthquake sequence: b-value estimation and plots
# -------------------------------------------------------------------------
# Run from the repository root directory. The catalog is expected in data/.
# The 400-event windows advance by 50 events. Uniform jittering is applied
# ONCE to the complete catalog, so overlapping windows share the same
# perturbed magnitude for each earthquake.


# -------------------------------------------------------------------------
# Working directory and input data
# -------------------------------------------------------------------------
# Run this script from the directory containing the R scripts.
# Place the file BDf2_tehuantepec_2017_100km.rds in the DATA/
# subdirectory before running the analysis.
# -------------------------------------------------------------------------
my_dir <- "."

source(file.path(my_dir, "extras_estimators.r"))
source(file.path(my_dir, "extras_plots.r"))

# Catalog in chronological order; identify the mainshock in the original data.
BDf2 <- readRDS(file.path(my_dir, "DATA/", "BDf2_tehuantepec_2017_100km.rds"))
BDf2 <- BDf2[order(BDf2$datetime_local), ]
rownames(BDf2) <- NULL
pm <- BDf2[which.max(BDf2$Magnitud), ]

# Estimation parameters.
a <- 0.05
bin <- 0.1
dc <- 0.2
ell <- 10

# Generate a single, reproducible jitter realization for the entire catalog.
set.seed(123)
BDf2$Magnitud_jit <- BDf2$Magnitud + bin*runif(nrow(BDf2))


# Overlapping event windows (400 events, moving forward by 50).
W400 <- make_event_windows(BDseq = BDf2, pm = pm, w = 400, step = 50)
k <- nrow(W400)

# Each result has columns: lower CI, point estimate, upper CI, effective n.
bGH <- array(NA_real_, c(k, 4))
bp  <- array(NA_real_, c(k, 4))
bmp <- array(NA_real_, c(k, 4))
bd  <- array(NA_real_, c(k, 4))
mc_vec <- rep(NA_real_, k)

colnames(bGH) <- c("lower", "bhat", "upper", "n_eff")
colnames(bp)  <- c("lower", "bhat", "upper", "n_eff")
colnames(bmp) <- c("lower", "bhat", "upper", "n_eff")
colnames(bd)  <- c("lower", "bhat", "upper", "n_eff")

for (j in seq_len(k)) {
  M0 <- W400$window[[j]]$Magnitud
  M00 <- M0

  # 1. Discrete Guttorp--Hopkins with window-specific maximum-curvature mc.
  mc <- mc_max_curvature(M00, bin = bin, correction = 0.2)
  mc_vec[j] <- mc
  mc <- mc
  M1 <- M00[M00 >= mc - 1e-10]
  nr <- length(M1)
  b_gut_hop_val <- b_gut_hop(nr, mean(M1), sd(M1), mc, bin, a = a)
  bGH[j, "bhat"]  <- b_gut_hop_val$bn
  bGH[j, "lower"] <- b_gut_hop_val$Le
  bGH[j, "upper"] <- b_gut_hop_val$Ue
  bGH[j, "n_eff"] <- nr

  # Use the SAME event-level jitter in every window containing an event.
  mags2 <- W400$window[[j]]$Magnitud_jit

  # 2. Continuous b-positive on jittered magnitudes.
  dfs <- diff(mags2)
  tdfs <- dfs[dfs > 0 & dfs >= dc - 1e-10]
  r <- length(tdfs)
  b_vde_val <- b_aki(r, mean(tdfs), sd(tdfs), dc, a = a)
  bp[j, "bhat"]  <- b_vde_val$bn
  bp[j, "lower"] <- b_vde_val$Le
  bp[j, "upper"] <- b_vde_val$Ue
  bp[j, "n_eff"] <- r

  # 3. Continuous b-more-positive (maximum forward lag ell), with jitter.
  res_bmp <- basic_stats_bmore_positive(mags2, ell, dc)
  b_bmp_val <- b_aki(res_bmp$rl, res_bmp$bbarl, res_bmp$bsdl, dc, a = a)
  bmp[j, "bhat"]  <- b_bmp_val$bn
  bmp[j, "lower"] <- b_bmp_val$Le
  bmp[j, "upper"] <- b_bmp_val$Ue
  bmp[j, "n_eff"] <- res_bmp$rl

  # 4. Discrete non-overlapping absolute differences, exact interval.
  npairs <- floor(length(M00) / 2)
  dfs2 <- M00[2 * seq_len(npairs)] - M00[2 * seq_len(npairs) - 1]
  dfs2 <- round(dfs2 / bin) * bin
  N <- length(dfs2)
  adbar <- mean(abs(dfs2))
  b_diff_disc_exact_val <- b_dif_disc_exact(N, adbar, bin, a = a)
  bd[j, "bhat"]  <- b_diff_disc_exact_val$bn
  bd[j, "lower"] <- b_diff_disc_exact_val$Le
  bd[j, "upper"] <- b_diff_disc_exact_val$Ue
  bd[j, "n_eff"] <- N
}

# Window estimates and effective sample sizes.
data.frame(time = W400$t_plot,
           bGH = round(bGH[, "bhat"], 1),
           bp = round(bp[, "bhat"], 1),
           bmp = round(bmp[, "bhat"], 1),
           bd = round(bd[, "bhat"], 1))

summary(mc_vec)
summary(bGH[, "n_eff"])
sd(bGH[, "n_eff"])
summary(bp[, "n_eff"])
sd(bp[, "n_eff"])
summary(bmp[, "n_eff"])
sd(bmp[, "n_eff"])
summary(bd[, "n_eff"])
sd(bd[, "n_eff"])

# Main figure: point estimates and confidence intervals.
dir.create(file.path(my_dir, "GRAF"), showWarnings = FALSE, recursive = TRUE)
pdf(file.path(my_dir, "GRAF", "b_graf.pdf"),
    width = 12, height = 4, pointsize = 12, paper = "special")
par(mfrow = c(1, 2), mai = c(0.85, 0.9, 0.05, 0.05))
plot_point_series(W400, bGH, bp, bmp, bd, pm,
                  ylim = c(1, 3.5), y_by = 0.5)
plot_interval_series(W400, bp, bd, pm,
                     ylim = c(1, 3.5), y_by = 0.5)
dev.off()

# Auxiliary plot: estimated completeness magnitude by window.
keep <- is.finite(bp[, 1]) & is.finite(bp[, 3]) &
        is.finite(bd[, 1]) & is.finite(bd[, 3])
tt <- W400$t_plot[keep]
mc <- mc_vec[keep]
ylim_mc <- range(mc, na.rm = TRUE)

plot(tt, mc,
     type = "b", pch = 16, lwd = 2, col = "gray20",
     xlab = "Months, 2017", ylab = expression(M[c]),
     xaxt = "n", yaxt = "n", ylim = ylim_mc, bty = "l")

y_ticks <- seq(floor(ylim_mc[1] * 10) / 10,
               ceiling(ylim_mc[2] * 10) / 10, by = 0.1)
axis(2, at = y_ticks, labels = y_ticks, las = 2,
     cex.axis = 0.85, col.axis = "gray25", col.ticks = "gray50")
abline(h = y_ticks, col = "gray92", lwd = 0.8)

# Month labels in English.
tz_use <- attr(tt, "tzone")
if (is.null(tz_use) || length(tz_use) == 0) {
  tz_use <- "America/Mexico_City"
}
x0 <- as.POSIXct(format(min(tt), "%Y-%m-01 00:00:00"), tz = tz_use)
x1 <- as.POSIXct(format(max(tt), "%Y-%m-01 00:00:00"), tz = tz_use)
at_month <- seq(x0, x1, by = "1 month")
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
month_labs <- month_names[as.integer(format(at_month, "%m"))]
axis.POSIXct(1, at = at_month, labels = month_labs, las = 2,
             cex.axis = 0.85, col.axis = "gray25", col.ticks = "gray50")
abline(v = at_month, col = "gray95", lwd = 0.8)
abline(v = pm$datetime_local[1], col = "red3", lty = 2, lwd = 1.3)
box(lwd = 2)
