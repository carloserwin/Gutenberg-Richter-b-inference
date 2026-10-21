#---------------------------------------------------
# Directory containing the R scripts
# Run this script from the repository's root directory.
my_dir <- "."
#---------------------------------------------------
source(file.path(my_dir, "extras_catalog_sim.r"))
source(file.path(my_dir, "extras_estimators.r"))

######################
# Catalog simulation #
######################
set.seed(232321)  # Reproducible catalog and jittering
b <- 1
mc <- 0
bin <- 0.1


mu <- 0.5
K <- 0.05
c0 <- 0.001
p0 <- 1.06
alpha <- 1.0
Mr <- mc
t0 <- 0
m_main <- 5
G <- 5.6
H <- 1
sigma <- 0.2
Mcat <- mc - 3*sigma

lev <- 1 # 1: complete catalog; 2: time-dependent incompleteness
n <- 2000
cat <- sim_catalog(lev, n, b, mc = mc, bin = bin,
                   mu = mu, K = K, c0 = c0,
                   p0 = p0, alpha = alpha, Mr = Mr,
                   t0 = t0, m_main = m_main,
                   G = G, H = H, Mcat = Mcat, sigma = sigma)
tt <- cat$time
# Binned magnitudes
mags <- cat$magnitude
# Jittering
mags_jit <- mags + bin*runif(length(mags))

# Calculating basic stats: magnitudes, abs-differences and b-positive methods
dc <- 0.2
basic_stats <- basic_magnitde_statistics(tt, mags, dc)
basic_stats_jit <- basic_magnitde_statistics(tt, mags_jit, dc)

# Check the statistics before using them in the estimators.
if (!basic_stats$valid) {
  stop("Invalid original-magnitude statistics: ", basic_stats$reason)
}
if (!basic_stats_jit$valid) {
  stop("Invalid jittered-magnitude statistics: ", basic_stats_jit$reason)
}

# Calculating basic stats: b-more-positive method
ell <- 10
basic_stats_bmp <- basic_stats_bmore_positive(basic_stats$magso, ell, dc)
basic_stats_bmp_jit <- basic_stats_bmore_positive(basic_stats_jit$magso, ell, dc)

if (!basic_stats_bmp$valid) {
  stop("Invalid b-more-positive statistics: ", basic_stats_bmp$reason)
}
if (!basic_stats_bmp_jit$valid) {
  stop("Invalid jittered b-more-positive statistics: ", basic_stats_bmp_jit$reason)
}


###############
# Basic stats #
###############
# magnitude-based
nr <- basic_stats$nr 
mbar <- basic_stats$mbar 
msd <- basic_stats$msd

nr_jit <- basic_stats_jit$nr
mbar_jit <- basic_stats_jit$mbar
msd_jit <- basic_stats_jit$msd

# magnitude-difference-based
N <- basic_stats$N 
adbar <- basic_stats$adbar 
s2_ad <- basic_stats$s2_ad

N_jit <- basic_stats_jit$N 
adbar_jit <- basic_stats_jit$adbar 
s2_ad_jit <- basic_stats_jit$s2_ad

# b-positive
r <- basic_stats$r 
bbar <- basic_stats$bbar 
bsd <- basic_stats$bsd

r_jit <- basic_stats_jit$r 
bbar_jit <- basic_stats_jit$bbar 
bsd_jit <- basic_stats_jit$bsd

# b-more-positive
rl <- basic_stats_bmp$rl
bbarl <- basic_stats_bmp$bbarl 
bsdl <- basic_stats_bmp$bsdl

rl_jit <- basic_stats_bmp_jit$rl
bbarl_jit <- basic_stats_bmp_jit$bbarl 
bsdl_jit <- basic_stats_bmp_jit$bsdl

###############
# Estimation  #
###############
a <- 0.05 # 95% confidence intervals
# Magnitude based estimators
b_aki_val <- b_aki(nr, mbar, msd, mc, a = a)
b_aki_val_jit <- b_aki(nr_jit, mbar_jit, msd_jit, mc, a = a)

b_bender_val <- b_bender_tgeom(nr, mbar, msd, mc, bin, mmax = 10, a = a, tol = 1e-8)
b_gut_hop_val <- b_gut_hop(nr, mbar, msd, mc, bin, a = a)
b_marz_val <- b_marzocchi(nr, mbar, msd, mc, bin, a = a)
# Difference based estimators
# Section 3.1.2
b_dif_val <- b_aki(N, adbar, sqrt(s2_ad), 0, a = a)
b_dif_val_jit <- b_aki(N_jit, adbar_jit, sqrt(s2_ad_jit), 0, a = a)
# Section 3.2.3
b_dif_disc_val <- b_dif_disc(N, adbar, s2_ad, bin, a)
b_dif_disc_exact_val <- b_dif_disc_exact(N, adbar, bin, a)


# b-positive
b_pos_cont_val <- b_aki(r, bbar, bsd, dc, a = a)
b_pos_cont_val_jit <- b_aki(r_jit, bbar_jit, bsd_jit, dc, a = a)
b_pos_disc_val <- b_gut_hop(r, bbar, bsd, dc, bin, a = a)
# b-more-positive
bm_pos_cont_val <- b_aki(rl, bbarl, bsdl, dc, a = a)
bm_pos_cont_val_jit <- b_aki(rl_jit, bbarl_jit, bsdl_jit, 
                             dc, a = a)


# -------------------------------------------------------------------------
# Results: all available confidence intervals for each estimator
# -------------------------------------------------------------------------
# C  = continuous estimator/model; D = discrete estimator/model.
# E  = exact under the stated reference model; MB = model-based asymptotic;
# NP = nonparametric asymptotic; jit = uniform jittering of binned magnitudes.
# The E label for b-positive and b-more-positive refers to the *working*
# independent exponential/geometric model, not the full dependent process.
# Marzocchi's intervals are plug-in intervals, NOT exact discrete intervals.

make_ci_rows <- function(method, fit, intervals) {
  data.frame(
    Method = paste(method, names(intervals), sep = "-"),
    bhat = rep(fit$bn, length(intervals)),
    L = vapply(intervals, function(fields) fit[[fields[1]]], numeric(1)),
    U = vapply(intervals, function(fields) fit[[fields[2]]], numeric(1)),
    row.names = NULL
  )
}

ci_all <- list(
  E  = c("Le", "Ue"),
  MB = c("La_model", "Ua_model"),
  NP = c("La_np", "Ua_np")
)
ci_asymptotic <- ci_all[c("MB", "NP")]
ci_marzocchi <- list(
  `Utsu-plugin` = c("L_utsu_plugin", "U_utsu_plugin")
)

results <- do.call(rbind, list(
  # Magnitude-based estimators
  make_ci_rows("Utsu-Aki-C",     b_aki_val,       ci_all),
  make_ci_rows("Utsu-Aki-C-jit", b_aki_val_jit,   ci_all),
  make_ci_rows("Bender-D",       b_bender_val,    ci_asymptotic),
  make_ci_rows("Guttorp-Hopkins-D", b_gut_hop_val, ci_all),
  make_ci_rows("Marzocchi-D",    b_marz_val,      ci_marzocchi),

  # Non-overlapping absolute-difference estimators
  make_ci_rows("AbsDiff-C",     b_dif_val,            ci_all),
  make_ci_rows("AbsDiff-C-jit", b_dif_val_jit,        ci_all),
  make_ci_rows("AbsDiff-D",     b_dif_disc_exact_val, ci_all["E"]),
  make_ci_rows("AbsDiff-D",     b_dif_disc_val,       ci_asymptotic),

  # b-positive estimators
  make_ci_rows("b-positive-C",     b_pos_cont_val,     ci_all),
  make_ci_rows("b-positive-C-jit", b_pos_cont_val_jit, ci_all),
  make_ci_rows("b-positive-D",     b_pos_disc_val,     ci_all),

  # b-more-positive estimators
  make_ci_rows("b-more-positive-C",     bm_pos_cont_val,     ci_all),
  make_ci_rows("b-more-positive-C-jit", bm_pos_cont_val_jit, ci_all)
))

stopifnot(nrow(results) == 36L)
print(results, row.names = FALSE, digits = 5)
