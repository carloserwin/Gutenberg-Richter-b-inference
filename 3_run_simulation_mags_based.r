#---------------------------------------------------
# Directory containing the R scripts
# Run this script from the repository's root directory.
my_dir <- "."
#---------------------------------------------------
source(file.path(my_dir, "extras_catalog_sim.r"))
source(file.path(my_dir, "extras_estimators.r"))
source(file.path(my_dir, "extras_simulation.r"))

# Directory to save the simulation results
dir.create(file.path(my_dir, "TAB"), showWarnings = FALSE)

set.seed(233434)  # Reproducibility


# -------------------------------------------------------------------------
# Simulation settings
# -------------------------------------------------------------------------
# Reproduces the Level 2 magnitude-based simulation study, with the
# completeness threshold estimated separately in each replication
# using the maximum-curvature method.
#
# m      : Number of Monte Carlo replications.
# n_grid : Complete catalogue sizes before applying the detection filter.
# -------------------------------------------------------------------------
m <- 10000
lev <- 2
n_grid <- c(100, 190, 400, 2120)

# -------------------------------------------------------------------------
# Model and estimation parameters
# -------------------------------------------------------------------------
# These values reproduce the simulation settings reported in the paper.
# See extras_catalog_sim.r for the catalogue generation procedures
# and extras_estimators.r for the estimation methods.
# -------------------------------------------------------------------------

# Significance level for confidence intervals  
a <- 0.05
# True Gutenberg--Richter parameters
b <- 1
mc <- 0
bin <- 0.1

# ETAS model parameters
mu <- 0.5
K <- 0.05
c0 <- 0.001
p0 <- 1.06
alpha <- 1.0
Mr <- mc

# Difference-based estimation parameters
dc <- 0.2   # Minimum retained magnitude difference
ell <- 10   # Maximum forward-search lag for b-more-positive

# -------------------------------------------------------------------------
# Run magnitude-based Monte Carlo simulations
# -------------------------------------------------------------------------
# For each catalogue size, record:
# - Execution time (minutes).
# - Average sample size after filtering by the estimated completeness level.
# - Average estimated completeness threshold.
#
# Simulation performance results are saved in the TAB/ directory.
# -------------------------------------------------------------------------
tobs <- rep(NA, 4)
nobsbar <- tobs
mcestbar <- tobs
for (j in 1:4) {
  n <- n_grid[j]
  t0 <- system.time(
    res <- run_simulation_magnitude_based(
      m, a, lev, n, b, ell,
      mc = mc, dc = dc, bin = bin,
      mu = mu, K = K, c0 = c0,
      p0 = p0, alpha = alpha, Mr = Mr,
      m_main = 5, Mcat = -0.6
    )
  )
  # Record simulation summaries
  tobs[j] <- t0[3] / 60
  nobsbar[j] <- mean(res$n0)
  mcestbar[j] <- mean(res$mc_est)
  # Save simulation results
  A <- round(rbind(res$p, res$pae, res$md), 1)
  rownames(A) <- c("TC", "PMAB", "bar_d")
  write.csv(t(A),
    file.path(my_dir, "TAB",
              paste0("RES_MAG_", lev, "_", n, ".csv")))
  # Display simulation progress
  cat("Simulation completed:\n")
  cat("Level:", lev, "\n")
  cat("Catalogue size:", n, "\n")
  cat("Mean retained sample size:", nobsbar[j], "\n")
  cat("Mean estimated Mc:", mcestbar[j], "\n")
  cat("Execution time (minutes):", tobs[j], "\n\n")
}