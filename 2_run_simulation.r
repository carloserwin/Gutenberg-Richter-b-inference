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
# m      : Number of Monte Carlo replications.
# jitter : 0 = use the original discretised magnitudes;
#          1 = apply uniform jittering to continuous estimators.
# lev    : 1 = Level 1 (complete ETAS--GR catalogues);
#          2 = Level 2 (transiently incomplete ETAS--GR catalogues).
#
# For Level 1, n_grid specifies the catalogue sizes.
# For Level 2, n_grid specifies the sizes of the complete catalogues
# before applying the time-dependent detection filter.
# -------------------------------------------------------------------------
m <- 10000
jitter <- 1
lev <- 2
if(lev == 1){
  n_grid <- c(20, 100, 300, 2000)
}else{
  n_grid <- c(100, 190, 400, 2120)
}
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
# Run Monte Carlo simulations
# -------------------------------------------------------------------------
# For each catalogue size, run the selected simulation scenario and
# record the execution time and average observed catalogue size.
#
# Results are saved as CSV files in the TAB/ directory.
# Each file contains TC, PMAB, and mean confidence interval length
# for the corresponding estimation methods.
# -------------------------------------------------------------------------
my_times <- rep(0, 4)
nn <- rep(0, 4)
for (j in 1:4) {
  n <- n_grid[j]
  if(jitter == 0){
    t0 <- system.time(res <- run_simulation(m, a, lev, n, b, ell, mc = mc, dc = dc, bin = bin,
                                            mu = mu, K = K, c0 = c0, p0 = p0, alpha = alpha, Mr = Mr, m_main = 5, Mcat = -0.6))
    simtypo <- "ALL"
  }else{
    t0 <- system.time(res <- run_simulation_cont_jitter(m, a, lev, n, b, ell, mc = mc, dc = dc, bin = bin,
                                                        mu = mu, K = K, c0 = c0, p0 = p0, alpha = alpha, Mr = Mr, m_main = 5, Mcat = -0.6))
    simtypo <- "JITTER"
  }
  # Record execution time (minutes)
  my_times[j] <- t0[3]/60
  # Display simulation progress
  cat("Simulation completed:\n")
  cat("Level:", lev, "\n")
  cat("Catalogue size:", n, "\n")
  cat("Type:", simtypo, "\n")
  cat("Execution time (minutes):", my_times[j], "\n")
  # Average number of observed events
  nn[j] <- mean(res$n0)
  # Save simulation results
  A <- round(rbind(res$p, res$pae, res$md), 1)
  rownames(A) <- c("TC", "PMAB", "bar_d")
  write.csv(t(A),
            file.path(my_dir, "TAB",
                      paste0("RES_", lev, "_", n, "_", simtypo, ".csv")))
}
my_times
nn
