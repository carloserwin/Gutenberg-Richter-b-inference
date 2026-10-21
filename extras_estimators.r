# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Inference over b - all the strategies
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 1- Aki--Utsu continuous exponential estimator
# -------------------------------------------------------------------------
# Computes Aki's continuous-exponential MLE for the b-value.
#
# Inputs are summary statistics from magnitudes above the threshold mc:
# sample size, sample mean, and sample standard deviation.
#
# Returns:
# - point estimate
# - Utsu exact chi-square confidence interval
# - model-based asymptotic confidence interval
# - Shi--Bolt nonparametric asymptotic confidence interval
#
# Input:
# n    - sample size
# mbar - sample mean of the retained magnitudes
# msd  - sample standard deviation of the retained magnitudes
# mc   - lower magnitude threshold
# a    - significance level; the confidence level is 1-a
b_aki <- function(n, mbar, msd, mc, a = 0.05) {
  # Aki's continuous exponential MLE
  bn <- 1 / (log(10) * (mbar - mc))
  # Utsu exact confidence interval under the exponential model
  Le <- (bn / (2 * n)) * qchisq(a / 2, df = 2 * n)
  Ue <- (bn / (2 * n)) * qchisq(1 - a / 2, df = 2 * n)
  # Model-based asymptotic standard error
  se_mb <- bn / sqrt(n)
  # Shi--Bolt nonparametric asymptotic standard error
  # Here msd should be the usual sample standard deviation, as returned by sd().
  se_np <- log(10) * bn^2 * msd / sqrt(n)
  # Wald confidence intervals
  z <- qnorm(1 - a / 2)
  La_model <- bn - z * se_mb
  Ua_model <- bn + z * se_mb
  La_np <- bn - z * se_np
  Ua_np <- bn + z * se_np
  list(
    bn = bn,
    Le = Le,
    Ue = Ue,
    La_model = La_model,
    Ua_model = Ua_model,
    La_np = La_np,
    Ua_np = Ua_np
  )
}
# -------------------------------------------------------------------------
# 2. Computes the Bender estimator for the b-value under a truncated geometric
# -------------------------------------------------------------------------
# The estimator is obtained by solving the
# likelihood equation for theta on the finite support {0, ..., J - 1}, where
# J is determined by mc, mmax, and the bin width. 
# The function returns the point estimate and two delta-method 
# Wald intervals: a model-based interval
# using the truncated-geometric variance, and a Shi-Bolt-type interval using
# the empirical magnitude standard deviation.
#
# Input:
# n     - sample size
# mbar  - sample mean of the retained discretized magnitudes
# msd   - sample standard deviation of the retained discretized magnitudes
# mc    - first retained magnitude on the discretized scale
# bin   - magnitude bin width Delta
# mmax  - fixed upper endpoint of the retained discretized support
# a - significance level for the confidence intervals
# tol   - numerical tolerance for grid alignment and root solving
b_bender_tgeom <- function(n, mbar, msd, mc, bin, mmax = 10,
                           a = 0.05, tol = 1e-8) {
  J_raw <- (mmax - mc) / bin + 1
  J <- round(J_raw)
  if (abs(J_raw - J) > tol) {
    stop("mmax and mc are not aligned with the discretized grid.")
  }
  if (mbar <= mc) {
    stop("The sample mean must be larger than mc for a finite estimate.")
  }
  theta_hat <- mle_theta_tgeom(
    mbar = mbar,
    mc = mc,
    bin = bin,
    mmax = mmax,
    tol = tol
  )
  # Boundary cases: the usual interior asymptotic approximation is not valid.
  if (theta_hat <= 0) {
    return(list(
      bn = Inf,
      beta_hat = Inf,
      theta_hat = theta_hat,
      se_model = NA_real_,
      La_model = NA_real_,
      Ua_model = NA_real_,
      se_np = NA_real_,
      La_np = NA_real_,
      Ua_np = NA_real_,
      n = n,
      J = J,
      mbar = mbar,
      mc = mc,
      mmax = mmax,
      bin = bin
    ))
  }
  
  if (theta_hat >= 1) {
    return(list(
      bn = 0,
      beta_hat = 0,
      theta_hat = theta_hat,
      se_model = NA_real_,
      La_model = NA_real_,
      Ua_model = NA_real_,
      se_np = NA_real_,
      La_np = NA_real_,
      Ua_np = NA_real_,
      n = n,
      J = J,
      mbar = mbar,
      mc = mc,
      mmax = mmax,
      bin = bin
    ))
  }
  
  beta_hat <- -log(theta_hat) / bin
  bn <- beta_hat / log(10)
  sigma2_hat <- var_tgeom(theta_hat, J)
  # Model-based delta-method standard error
  se_model <- 1 / (bin * log(10) * sqrt(n * sigma2_hat))
  # Shi-Bolt-type standard error:
  # estimate Var(X') nonparametrically through S_M.
  se_np <- msd / (bin^2 * log(10) * sigma2_hat * sqrt(n))
  z <- qnorm(1 - a / 2)
  list(
    bn = bn,
    beta_hat = beta_hat,
    theta_hat = theta_hat,
    se_model = se_model,
    La_model = bn - z * se_model,
    Ua_model = bn + z * se_model,
    se_np = se_np,
    La_np = bn - z * se_np,
    Ua_np = bn + z * se_np,
    n = n,
    J = J,
    mbar = mbar,
    mc = mc,
    mmax = mmax,
    bin = bin
  )
}
# Mean and variance of truncated geometric on {0, ..., J - 1}
mu_tgeom <- function(theta, J) {
  if (theta <= 0) return(0)
  if (theta >= 1) return((J - 1) / 2)
  
  theta / (1 - theta) -
    J * theta^J / (1 - theta^J)
}
# Variance of the truncated geometric distribution on {0, ..., J - 1}.
var_tgeom <- function(theta, J) {
  if (theta <= 0) return(0)
  if (theta >= 1) return((J^2 - 1) / 12)
  
  theta / (1 - theta)^2 -
    J^2 * theta^J / (1 - theta^J)^2
}
# Solves the Bender estimating equation
# (mbar - mc) / bin = mu_J(theta)
# to obtain the MLE of theta.
mle_theta_tgeom <- function(mbar, mc, bin, mmax, tol = 1e-10) {
  
  J <- round((mmax - mc) / bin + 1)
  xbar <- (mbar - mc) / bin
  upper_mean <- (J - 1) / 2
  if (xbar <= tol) return(0)
  if (xbar >= upper_mean - tol) return(1)
  uniroot(
    f = function(theta) mu_tgeom(theta, J) - xbar,
    interval = c(.Machine$double.eps, 1 - .Machine$double.eps^0.5),
    tol = tol
  )$root
}
# -------------------------------------------------------------------------
# 3. Guttorp--Hopkins estimator 
# -------------------------------------------------------------------------
# Computes the Guttorp--Hopkins MLE for the b-value under the untruncated
# discretized exponential model.
#
# Inputs are summary statistics from magnitudes above the threshold mc:
# sample size, sample mean, sample standard deviation, threshold and bin width.
#
# Returns:
# - point estimate
# - exact Tinti--Mulargia confidence interval via beta quantiles
# - model-based asymptotic delta-method confidence interval
# - nonparametric asymptotic delta-method confidence interval
#
# Input:
# n     - sample size
# mbar  - sample mean of the retained discretized magnitudes
# msd   - sample standard deviation of the retained discretized magnitudes
# mc    - first retained magnitude on the discretized scale
# bin   - magnitude bin width Delta
# a     - significance level; the confidence level is 1-a
b_gut_hop <- function(n, mbar, msd, mc, bin, a = 0.05) {
  d <- mbar - mc
  if (d <= 0) {
    stop("The sample mean must be larger than mc for a finite estimate.")
  }
  # Guttorp--Hopkins estimator
  bn <- 1 / (bin * log(10)) * log(1 + bin / d)
  # MLE of theta = exp(-beta Delta)
  theta_hat <- d / (d + bin)
  # Observed sufficient statistic:
  # t_obs = sum X_i' = n (mbar - mc) / Delta
  t_obs <- round(n * d / bin)
  # Exact Tinti--Mulargia interval via beta quantiles
  theta_L <- if (t_obs == 0) {
    0
  } else {
    qbeta(a / 2, shape1 = t_obs, shape2 = n)
  }
  theta_U <- qbeta(1 - a / 2, shape1 = t_obs + 1, shape2 = n)
  Le <- -log(theta_U) / (bin * log(10))
  Ue <- if (theta_L == 0) {
    Inf
  } else {
    -log(theta_L) / (bin * log(10))
  }
  # Asymptotic delta-method standard errors
  se_model <- 1 / (log(10) * sqrt(n * d * (d + bin)))
  se_np <- msd / (log(10) * d * (d + bin) * sqrt(n))
  z <- qnorm(1 - a / 2)
  La_model <- bn - z * se_model
  Ua_model <- bn + z * se_model
  La_np <- bn - z * se_np
  Ua_np <- bn + z * se_np
  list(
    bn = bn,
    Le = Le,
    Ue = Ue,
    La_model = La_model,
    Ua_model = Ua_model,
    La_np = La_np,
    Ua_np = Ua_np
  )
}
# -------------------------------------------------------------------------
# Marzocchi et al. binning-corrected estimator
# -------------------------------------------------------------------------
# Computes the Marzocchi et al. binning-corrected b-value estimator.
# The point estimate is algebraically equivalent to the Guttorp--Hopkins
# discretized-exponential estimator.
#
# Uncertainty is attached pragmatically by plugging the corrected estimate
# into the Utsu chi-square interval form and into a Shi--Bolt-type
# nonparametric asymptotic Wald interval.
#
# Input:
# n     - sample size
# mbar  - sample mean of the retained magnitudes
# msd   - sample standard deviation of the retained magnitudes
# mc    - first retained magnitude on the discretized scale
# bin   - magnitude bin width Delta
# a     - significance level; the confidence level is 1-a
b_marzocchi <- function(n, mbar, msd, mc, bin, a = 0.05) {
  d <- mbar - mc
  if (d <= 0) {
    stop("The sample mean must be larger than mc for a finite estimate.")
  }
  # Marzocchi et al. corrected estimator
  # Algebraically equivalent to the Guttorp--Hopkins estimator
  bn <- 1 / (bin * log(10)) * log(1 + bin / d)
  # Utsu-form plug-in interval using the corrected estimate
  L_utsu_plugin <- (bn / (2 * n)) * qchisq(a / 2, df = 2 * n)
  U_utsu_plugin <- (bn / (2 * n)) * qchisq(1 - a / 2, df = 2 * n)
  # Shi--Bolt-type plug-in nonparametric asymptotic standard error
  se_np <- log(10) * bn^2 * msd / sqrt(n)
  z <- qnorm(1 - a / 2)
  L_np <- bn - z * se_np
  U_np <- bn + z * se_np
  
  list(
    bn = bn,
    L_utsu_plugin = L_utsu_plugin,
    U_utsu_plugin = U_utsu_plugin,
    L_np = L_np,
    U_np = U_np,
    se_np = se_np
  )
}
# -------------------------------------------------------------------------
# 4.1 Discrete absolute-difference estimator
# model based and non-parametric intervals
# -------------------------------------------------------------------------
# Computes the b-value estimator based on absolute binned magnitude
# differences under the folded discrete Laplace model.
#
# Inputs:
# N      - number of absolute differences
# adbar  - sample mean of the absolute differences, on the magnitude scale
# s2_ad  - sample variance of the absolute differences, on the magnitude scale
# bin    - magnitude bin width Delta
# a      - significance level; the confidence level is 1-a
b_dif_disc <- function(N, adbar, s2_ad, bin = 0.1, a = 0.05) {
  
  ybar <- adbar / bin
  # Stable version of theta_hat = (sqrt(1 + ybar^2) - 1) / ybar
  theta_hat <- ybar / (sqrt(1 + ybar^2) + 1)
  # Point estimate
  bn <- -log(theta_hat) / (bin * log(10))
  z <- qnorm(1 - a / 2)
  # |g'(adbar)| = 1 / {log(10) adbar sqrt(adbar^2 + bin^2)}
  delta_factor <- 1 / (log(10) * adbar * sqrt(adbar^2 + bin^2))
  # Model-based interval
  varA_model <- adbar * sqrt(adbar^2 + bin^2)
  se_model <- delta_factor * sqrt(varA_model / N)
  La_model <- bn - z * se_model
  Ua_model <- bn + z * se_model
  # Nonparametric interval
  se_np <- delta_factor * sqrt(s2_ad / N)
  La_np <- bn - z * se_np
  Ua_np <- bn + z * se_np
  list(bn = bn,
       La_model = La_model,
       Ua_model = Ua_model,
       La_np = La_np,
       Ua_np = Ua_np)
}
# -------------------------------------------------------------------------
# 4.2 Discrete absolute-difference estimator
# exact confidence interval
# -------------------------------------------------------------------------
# The estimator is based on absolute binned magnitude differences
#
#   A_i' = |D_i'|,
#
# and on the statistic
#
#   T_N = sum_i A_i' / Delta,
#
# where Delta is the magnitude bin width. Under the folded discrete Laplace
# model, T_N has an exact finite-sample distribution depending on theta.
# The interval is obtained by pivoting its CDF, in the same spirit as the
# exact interval of Tinti and Mulargia for binned magnitudes.
#
# Arguments:
#   N        : number of absolute differences.
#   adbar    : sample mean of A_i', on the magnitude scale.
#   bin      : magnitude bin width Delta.
#   a        : significance level; default gives a 95% confidence interval.
#   tol_grid : tolerance used to check that N * adbar / bin is an integer.
#   tol_root : numerical tolerance used in one-dimensional root finding.
#
# Returns:
#   A list with the point estimate, the exact confidence interval for b,
#   the corresponding theta limits, and auxiliary quantities.
#
# Numerical note:
# For large N, pnbinom() may produce internal pbeta() underflow
# warnings when evaluating extremely small tail probabilities.
# These warnings typically arise in mixture components with
# negligible binomial weights and do not necessarily affect
# the resulting CDF.
# The numerical accuracy of the interval endpoints should
# be checked by evaluating the CDF at the computed roots.
# with(b_dif_disc_exact_val, {
#  c(lower = pT_absdisc(t_obs - 1, N, theta_L),
#    upper = pT_absdisc(t_obs, N, theta_U)
#    )})
b_dif_disc_exact <- function(N, adbar, bin = 0.1, a = 0.05,
                              tol_grid = 1e-8,
                              tol_root = 1e-10) {
  # Observed value of
  #
  #   T_N = sum_i A_i' / bin.
  #
  # Since adbar = mean(A_i'), this is N * adbar / bin.
  t_raw <- N * adbar / bin
  t_obs <- round(t_raw)
  # In exact arithmetic, t_raw should be an integer. The tolerance only
  # protects against floating-point representation errors.
  if (abs(t_raw - t_obs) > tol_grid) {
    warning("N * adbar / bin is not close to an integer; using round(N * adbar / bin).")
  }
  t_obs <- as.integer(t_obs)
  # Point estimator.
  #
  # If ybar = adbar / bin, the MLE of theta is
  #
  #   theta_hat = (sqrt(1 + ybar^2) - 1) / ybar.
  #
  # We use the algebraically equivalent but numerically stable form below.
  ybar <- adbar / bin
  if (ybar == 0) {
    theta_hat <- 0
    beta_hat <- Inf
    bn <- Inf
  } else {
    theta_hat <- ybar / (sqrt(1 + ybar^2) + 1)
    beta_hat <- -log(theta_hat) / bin
    bn <- beta_hat / log(10)
  }
  # Small numerical buffer for root finding on (0, 1)
  eps <- 1e-12
  # Lower theta limit.
  #
  # The CDF F_T(t | theta) is decreasing in theta for fixed t.
  # The lower equal-tail limit solves
  #
  #   F_T(t_obs - 1 | theta_L) = 1 - a / 2.
  #
  # If t_obs = 0, the convention is theta_L = 0.
  if (t_obs == 0L) {
    theta_L <- 0
  } else {
    f_lower <- function(theta) {
      pT_absdisc(t_obs - 1, N, theta) - (1 - a / 2)
    }
    theta_L <- uniroot(
      f_lower,
      interval = c(eps, 1 - eps),
      tol = tol_root
    )$root
  }
  # Upper theta limit.
  #
  # The upper equal-tail limit solves
  #
  #   F_T(t_obs | theta_U) = a / 2.
  f_upper <- function(theta) {
    pT_absdisc(t_obs, N, theta) - a / 2
  }
  theta_U <- uniroot(
    f_upper,
    interval = c(eps, 1 - eps),
    tol = tol_root)$root
  # Transform the theta interval into a b interval.
  #
  # Since
  #
  #   b = -log(theta) / (bin log(10))
  #
  # is decreasing in theta, the endpoints are reversed.
  L_exact <- -log(theta_U) / (bin * log(10))
  
  if (theta_L == 0) {
    U_exact <- Inf
  } else {
    U_exact <- -log(theta_L) / (bin * log(10))
  }
  
  list(
    bn = bn,
    beta_hat = beta_hat,
    theta_hat = theta_hat,
    
    # Exact confidence interval for b
    Le = L_exact,
    Ue = U_exact,
    
    # Equal-tail limits in theta
    theta_L = theta_L,
    theta_U = theta_U,
    
    # Inputs and auxiliary quantities
    N = N,
    adbar = adbar,
    ybar = ybar,
    t_obs = t_obs,
    bin = bin,
    alpha = a
  )
}

# Numerically stable computation of log(sum(exp(x))).
# Used to add probabilities computed on the log scale.
log_sum_exp <- function(x) {
  mmm <- max(x)
  if (!is.finite(mmm)) {
    return(mmm)
  }
  mmm + log(sum(exp(x - mmm)))
}

# CDF of T_N = sum_i A_i' / Delta under the folded discrete Laplace model.
#
# The computation uses the binomial--negative-binomial mixture representation:
# first condition on the number of non-zero absolute differences.
pT_absdisc <- function(tt, N, theta, lower.tail = TRUE, log.p = FALSE) {
  
  if (length(tt) != 1L || !is.finite(tt)) {
    stop("tt must be a finite number.")
  }
  
  if (length(N) != 1L || !is.finite(N) || N <= 0 || N != floor(N)) {
    stop("N must be a positive integer.")
  }
  
  if (length(theta) != 1L || !is.finite(theta) || theta < 0 || theta >= 1) {
    stop("theta must be in [0, 1).")
  }
  
  tt <- floor(tt)
  
  if (tt < 0) {
    p <- 0
    if (!lower.tail) p <- 1
    if (log.p) p <- log(p)
    return(p)
  }
  
  if (theta == 0) {
    p <- 1
    if (!lower.tail) p <- 0
    if (log.p) p <- log(p)
    return(p)
  }
  
  rmax <- min(N, tt)
  r <- 0:rmax
  
  # Probability that an absolute difference is non-zero
  p_pos <- 2 * theta / (1 + theta)
  
  # R_N = number of non-zero absolute differences
  log_terms <- dbinom(r, size = N, prob = p_pos, log = TRUE)
  
  # Conditional on R_N = r > 0,
  # T_N - r has a negative-binomial distribution with parameters
  # size = r and prob = 1 - theta.
  log_cond <- rep(0, length(r))
  
  id <- r > 0
  
  if (any(id)) {
    log_cond[id] <- pnbinom(
      q = tt - r[id],
      size = r[id],
      prob = 1 - theta,
      log.p = TRUE
    )
  }
  
  log_cdf <- log_sum_exp(log_terms + log_cond)
  cdf <- exp(log_cdf)
  
  # Guard against tiny numerical errors
  cdf <- min(max(cdf, 0), 1)
  
  if (!lower.tail) {
    cdf <- 1 - cdf
  }
  
  if (log.p) {
    cdf <- log(cdf)
  }
  
  cdf
}
# -------------------------------------------------------------------------
# Basic summary statistics for magnitude-based and difference-based methods
# -------------------------------------------------------------------------
# Computes the summary statistics needed by the inference functions up to,
# but not including, the b-more-positive pre-processing.
#
# The function assumes that:
# - tt contains the occurrence times of the events.
# - mags contains the corresponding magnitudes.
# - dc is the minimum retained positive-difference threshold.
#
# The magnitudes are first summarized directly, giving the sample size,
# sample mean and sample standard deviation used by magnitude-based
# estimators.
#
# Then the catalogue is ordered by occurrence time. From the ordered
# magnitudes, the function computes:
#
# 1. Successive positive differences:
#      D_i = M_{i+1} - M_i,
#    retaining only differences satisfying D_i > 0 and D_i >= dc.
#    These summaries are used for the b-positive estimator.
#
# 2. Non-overlapping absolute differences:
#      D_i = M_{2i} - M_{2i-1},
#      A_i = |D_i|,
#    for i = 1,...,floor(n/2).
#    These summaries are used for the continuous and discrete
#    absolute-difference estimators.
#
# The function returns valid = TRUE when all required quantities are finite
# and usable. If any required check fails, it returns valid = FALSE together
# with a short reason. This allows the main simulation loop to stop or handle
# the failed replicate before calling the inference functions.
#
# Output when valid = TRUE:
#   valid - logical indicator.
#   reason - NULL.
#   nr - number of magnitudes.
#   mbar - mean magnitude.
#   msd - sample standard deviation of magnitudes.
#   r - number of retained successive positive differences.
#   bbar - mean retained positive difference.
#   bsd - sample standard deviation of retained positive differences.
#   N - number of non-overlapping differences.
#   adbar - mean absolute non-overlapping difference.
#   s2_ad - sample variance of absolute non-overlapping differences.
#
# Output when valid = FALSE:
#   valid - FALSE.
#   reason - character string describing the failed check.
# -------------------------------------------------------------------------
basic_magnitde_statistics <- function(tt, mags, dc){
  valid <- TRUE
  reason <- NULL
  invalid <- function(msg){
    valid <- FALSE
    list(valid = valid, reason = msg)
  }
  
  nr <- length(mags)
  if (nr < 8L) return(invalid("nr < 8"))
  if (!all(is.finite(tt)) || !all(is.finite(mags))) {
    return(invalid("non-finite tt or mags"))
  }
  
  mbar <- mean(mags)
  msd <- sd(mags)
  if (!is.finite(mbar) || !is.finite(msd) || msd <= 0) {
    return(invalid("invalid mbar or msd"))
  }
  o <- order(tt)
  mags2 <- mags[o]
  # Successive positive differences (b-positive)
  # Allow for floating-point errors at the threshold dc.
  dfs <- diff(mags2)
  tdfs <- dfs[dfs > 0 & dfs >= dc - 1e-10]
  r <- length(tdfs)
  if (r < 2L) return(invalid("r < 2"))
  if (!all(is.finite(tdfs))) return(invalid("non-finite tdfs"))
  
  bbar <- mean(tdfs)
  bsd <- sd(tdfs)
  # Check that the retained excesses yield a finite estimate.
  if (!is.finite(bbar) || !is.finite(bsd) ||
      bbar <= dc + 1e-10 || bsd <= 1e-10) {
    return(invalid("invalid b-positive statistics"))
  }
  
  # Non-overlapping differences
  npairs <- floor(length(mags2) / 2)
  if (npairs < 2L) return(invalid("npairs < 2"))
  dfs2 <- mags2[2 * seq_len(npairs)] - mags2[2 * seq_len(npairs) - 1]
  N <- length(dfs2)
  if (N < 2L) return(invalid("N < 2"))
  if (!all(is.finite(dfs2))) return(invalid("non-finite dfs2"))
  
  ad <- abs(dfs2)
  if (length(ad) < 2L || any(!is.finite(ad))) {
    return(invalid("invalid ad"))
  }
  
  adbar <- mean(ad)
  s2_ad <- var(ad)
  
  if (!is.finite(adbar) || adbar <= 0) {
    return(invalid("invalid adbar"))
  }
  if (!is.finite(s2_ad) || s2_ad < 0) {
    return(invalid("invalid s2_ad"))
  }
  
  list(valid = valid,
       reason = reason,
       magso = mags2,
       nr = nr, mbar = mbar, msd = msd,
       r = r, bbar = bbar, bsd = bsd,
       N = N, adbar = adbar, s2_ad = s2_ad)
}
# -------------------------------------------------------------------------
# # Pre-processing for the continuous maximum-lag b-more-positive estimator
# -------------------------------------------------------------------------
# For each event i, search forward up to ell events and select the first
# event j such that M_j > M_i. The difference D = M_j - M_i is retained
# only if D >= dc.
#
# Inputs:
# M   - vector of magnitudes, already ordered by occurrence time
# ell - maximum forward lag
# dc  - minimum positive difference threshold
#
# Output:
# vector of retained b-more-positive differences
get_b_more_positive_diffs <- function(M, ell = 5, dc = 0) {
  n <- length(M)
  if (n < 2) {
    return(numeric(0))
  }
  out <- numeric(0)
  for (i in 1:(n - 1)) {
    jmax <- min(n, i + ell)
    found <- FALSE
    j <- i + 1
    while (!found && j <= jmax) {
      if (M[j] > M[i]) {
        d <- M[j] - M[i]
        # Allow for floating-point errors at the threshold dc.
        if (d >= dc - 1e-10) {
          out <- c(out, d)
        }
        found <- TRUE
      }
      j <- j + 1
    }
  }
  out
}
# -------------------------------------------------------------------------
# Summary statistics for the b-more-positive estimator
# -------------------------------------------------------------------------
# Computes the summary statistics of the positive magnitude differences
# selected by the b-more-positive rule with maximum forward lag ell
# and minimum difference threshold dc.
#
# The input magnitudes must already be ordered by occurrence time.
# The selected differences are obtained using get_b_more_positive_diffs().
#
# Input:
# magso - magnitudes ordered by occurrence time
# ell   - maximum number of subsequent events examined
# dc    - minimum retained positive-difference threshold
#
# Output:
# valid - TRUE if the retained differences yield usable statistics
# reason - explanation when valid = FALSE; otherwise NULL
# rl    - number of retained b-more-positive differences
# bbarl - sample mean of the retained differences
# bsdl  - sample standard deviation of the retained differences
#
# At least two retained differences with positive sample standard
# deviation are required. Otherwise, valid = FALSE is returned.
# -------------------------------------------------------------------------
basic_stats_bmore_positive <- function(magso, ell, dc){
  valid <- TRUE
  reason <- NULL
  D_more <- get_b_more_positive_diffs(magso, ell = ell, dc = dc)
  rl <- length(D_more)
  bbarl <- NA_real_
  bsdl <- NA_real_
  if (rl < 2L) {
    valid <- FALSE
    reason <- "rl < 2"
  } else {
    bbarl <- mean(D_more)
    bsdl <- sd(D_more)
    # Check that the retained excesses yield a finite estimate.
    if (!is.finite(bbarl) || !is.finite(bsdl) ||
        bbarl <= dc + 1e-10 || bsdl <= 1e-10) {
      valid <- FALSE
      reason <- "invalid b-more-positive statistics"
    }
  }
  
  list(
    valid = valid,
    reason = reason,
    rl = rl,
    bbarl = bbarl,
    bsdl = bsdl
  )
}
# -------------------------------------------------------------------------
# Completeness magnitude estimation: maximum curvature method
# -------------------------------------------------------------------------
# Estimates the magnitude of completeness mc using the maximum curvature
# method applied to the frequency-magnitude distribution.
#
# Magnitudes are first rounded to the specified bin width.
# The magnitude class with the highest observed frequency is selected,
# and a fixed correction is added to obtain the estimated completeness
# magnitude. The resulting estimate is rounded to the magnitude grid.
#
# Input:
# M          - vector of observed magnitudes
# bin        - magnitude bin width
# correction - additive correction to the modal magnitude
#
# Output:
# Estimated completeness magnitude mc on the discretized scale.
# -------------------------------------------------------------------------
mc_max_curvature <- function(M, bin = 0.1, correction = 0.2){
  Mb <- round(M / bin) * bin
  vals <- sort(unique(Mb))
  tab <- tabulate(match(Mb, vals))
  mc0 <- vals[which.max(tab)]
  mc <- mc0 + correction
  round(mc / bin) * bin
}