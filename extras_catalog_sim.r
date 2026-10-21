# -------------------------------------------------------------------------
# 0) Simulate a seismic catalog
# -------------------------------------------------------------------------
#
# Wrapper function to simulate either a complete ETAS-GR catalog
# or an incomplete ETAS-GR catalog with time-dependent detection.
#
# Level 1:
#   Simulates a complete fixed-n ETAS-GR catalog.
#   Magnitudes are generated from a binned Gutenberg-Richter distribution,
#   and occurrence times are generated from an ETAS-type sequential model.
#
# Level 2:
#   First simulates a complete fixed-n ETAS-GR catalog as in Level 1.
#   Then applies a time-dependent incompleteness filter, retaining each
#   event with probability
#
#     p_det(m_i, t_i) = Phi{ (m_i - Mc(t_i)) / sigma }.
#
# Input:
# lev    - simulation level: 1 for complete catalogs, 2 for incomplete catalogs
# n      - number of events in the complete catalog
# b      - Gutenberg-Richter b-value used to generate magnitudes
# mc     - first retained magnitude on the discretized scale
# bin    - magnitude bin width
# mu     - background seismicity rate in the ETAS model
# K      - productivity parameter of the ETAS triggering kernel
# c0     - Omori-Utsu offset parameter; avoids singularity at zero
# p0     - Omori-Utsu temporal decay parameter
# alpha  - magnitude productivity parameter
# Mr     - reference magnitude for ETAS productivity
# t0     - initial time
# m_main - optional magnitude of an initial mainshock at time t0
# G      - magnitude offset in the time-dependent completeness model
# H      - temporal decay parameter in the completeness model
# Mcat   - background/network completeness magnitude
# sigma  - smoothness parameter of the detection probability
#
# Output:
# A data frame with the simulated catalog.
# For lev = 1, the output is the complete catalog.
# For lev = 2, the output is the observed catalog after filtering.
sim_catalog <- function(lev, n, b, mc = 0, bin = 0.1,
                        mu = 0.5, K = 0.05, c0 = 0.01,
                        p0 = 1.1, alpha = 1.0, Mr = mc,
                        t0 = 0, m_main = NULL,
                        G = 5.6, H = 1,
                        Mcat = 0.95, sigma = 0.2) {
  
  if (!lev %in% c(1, 2)) {
    stop("lev must be either 1 or 2.")
  }
  
  if (lev == 1) {
    
    catalog <- simulate_catalog_level1(
      n = n,
      b = b,
      mc = mc,
      bin = bin,
      mu = mu,
      K = K,
      c0 = c0,
      p0 = p0,
      alpha = alpha,
      Mr = Mr,
      t0 = t0,
      m_main = m_main
    )
  } else {
    catalog <- simulate_catalog_level2(
      n = n,
      b = b,
      mc = mc,
      bin = bin,
      mu = mu,
      K = K,
      c0 = c0,
      p0 = p0,
      alpha = alpha,
      Mr = Mr,
      t0 = t0,
      m_main = m_main,
      G = G,
      H = H,
      Mcat = Mcat,
      sigma = sigma,
      log_base = 10,
      return_complete = FALSE
    )
  }
  
  catalog
}
# ------------------------------------------------------------
# Level 1: Simulate a fixed-n ETAS-GR catalog, Ogata 1988
# ------------------------------------------------------------
#
# Simulates a synthetic earthquake catalog with exactly n events.
# Magnitudes are independent GR marks. Each magnitude is independent
# of its own occurrence time conditional on the preceding history,
# but affects the occurrence rates of subsequent events.
#
# The ETAS conditional rate is
#
# lambda(t) = mu + sum_{t_j < t}
#             K*exp{alpha*(m_j - Mr)} / (t - t_j + c0)^p0.
#
# Thus, each past event contributes to the future occurrence rate,
# with productivity depending on its magnitude and temporal decay
# following an Omori-Utsu type kernel.
#
# Input:
# n       - number of events to simulate
# b       - Gutenberg-Richter b-value used to generate magnitudes
# mc      - first retained magnitude on the discretized scale
# bin     - magnitude bin width
# mu      - background seismicity rate
# K       - productivity parameter of the triggering kernel
# c0      - Omori-Utsu offset parameter; avoids singularity at zero
# p0      - Omori-Utsu temporal decay parameter
# alpha   - magnitude productivity parameter
# Mr      - reference magnitude for productivity
# t0      - initial time
# m_main  - optional magnitude of an initial mainshock at time t0
#
# If m_main is NULL, the catalog starts with an empty history.
# If m_main is supplied, the simulation is conditional on an initial
# mainshock at time t0, but this mainshock is not included in the
# returned catalog.
#
# Output:
# A data frame with columns:
# id        - event index
# time      - simulated occurrence time
# magnitude - simulated binned magnitude
simulate_catalog_level1 <- function(n, b, mc, bin = 0.1,
                                    mu = 0.5, K = 0.05, c0 = 0.01,
                                    p0 = 1.1, alpha = 1.0, Mr = mc,
                                    t0 = 0, m_main = NULL) {
  if (mu <= 0) {
    stop("This fixed-n sequential simulator assumes mu > 0.")
  }
  
  if (c0 <= 0) {
    stop("Parameter c must be positive.")
  }
  # Generate exactly n magnitudes first.
  m <- simulate_magnitudes_gr(n, b, mc = mc, bin = bin)
  t <- numeric(n)
  
  # History used by the ETAS intensity.
  # If m_main is supplied, we condition on a mainshock at t0.
  if (is.null(m_main)) {
    t_hist <- numeric(0)
    m_hist <- numeric(0)
  } else {
    t_hist <- t0
    m_hist <- m_main
  }
  
  t_prev <- t0
  for (i in seq_len(n)) {
    W_i <- simulate_etas_waiting_time(
      t_prev = t_prev,
      t_hist = t_hist,
      m_hist = m_hist,
      mu = mu,
      K = K,
      c0 = c0,
      p0 = p0,
      alpha = alpha,
      Mr = Mr
    )
    
    t[i] <- t_prev + W_i
    
    # Update history: the new event affects future events.
    t_hist <- c(t_hist, t[i])
    m_hist <- c(m_hist, m[i])
    
    t_prev <- t[i]
  }
  
  data.frame(
    id = seq_len(n),
    time = t,
    magnitude = m
  )
}
# Simulates binned magnitudes from a Gutenberg-Richter exponential tail.
#
# M' = mc + bin*X',  P(X'=x) = (1-theta)*theta^x,
# theta = exp(-beta*bin), beta = b*log(10).
#
# Input:
# n   - sample size
# b   - Gutenberg-Richter b-value
# mc  - first retained magnitude on the discretized scale
# bin - magnitude bin width
#
# mc should be on the same lattice as the binned magnitudes,
# i.e. mc/bin should be an integer.
simulate_magnitudes_gr <- function(n, b, mc = 0, bin = 0.1) {
  if (abs(mc/bin - round(mc/bin)) > 1e-8) {
    warning("mc is not aligned with the magnitude bin width.")
  }
  bbeta <- b*log(10)
  theta <- exp(-bbeta*bin)
  xs <- rgeom(n, prob = 1 - theta)
  mags <- mc + bin*xs
  mags
}
#  Integrated conditional rate H_i(w)
etas_H <- function(w, t_prev, t_hist, m_hist,
                   mu, K, c0, p0, alpha, Mr) {
  
  H <- mu * w
  if (length(t_hist) == 0) {
    return(H)
  }
  d <- t_prev - t_hist + c0
  A <- K*exp(alpha*(m_hist - Mr))
  if (abs(p0 - 1) < 1e-10) {
    trigger_int <- log((d + w) / d)
  } else {
    trigger_int <- ((d + w)^(1 - p0) - d^(1 - p0)) / (1 - p0)
  }
  H + sum(A * trigger_int)
}
# Simulate one ETAS waiting time
simulate_etas_waiting_time <- function(t_prev, t_hist, m_hist,
                                       mu, K, c0, p0, alpha, Mr) {
  
  E <- rexp(1, rate = 1)
  f_root <- function(w) {
    etas_H(w, t_prev, t_hist, m_hist,
           mu = mu, K = K, c0 = c0, p0 = p0,
           alpha = alpha, Mr = Mr) - E
  }
  # Initial upper bound.
  # If there were no triggering, W would be approximately E / mu.
  upper <- max(E / mu, 1e-8)
  while (f_root(upper) < 0) {
    upper <- 2 * upper
  }
  uniroot(f_root, lower = 0, upper = upper, tol = 1e-12)$root
}
# ------------------------------------------------------------
# Level 2: ETAS-GR complete catalog + time-dependent incompleteness
# ------------------------------------------------------------
#
# Simulates a complete fixed-n ETAS-GR catalog using simulate_catalog_level1()
# and then applies a time-dependent detection filter to mimic short-term
# catalog incompleteness during active earthquake sequences.
#
# The complete catalog is first generated from an ETAS-type temporal model
# with Gutenberg-Richter magnitudes. Then, for each event time t_i, an
# instantaneous completeness magnitude Mc(t_i) is computed from the previous
# events. The detection probability of event i is modeled as
#
# p_det(m_i, t_i) = Phi{ (m_i - Mc(t_i)) / sigma },
#
# where Phi is the standard normal distribution function. Events are retained
# independently with probability p_det(m_i, t_i).
#
# The time-dependent completeness magnitude is based on the rule
#
# Mc(t) = max{ Mcat, max_{t_j < t} [m_j - G - H*log(t - t_j)] },
#
# with the logarithm taken in base log_base. This construction makes
# completeness worse immediately after large events and allows it to improve
# gradually with time.
#
# Input:
# n               - number of events in the complete catalog
# b               - Gutenberg-Richter b-value used to generate magnitudes
# mc              - first retained magnitude on the discretized scale
# bin             - magnitude bin width
# mu              - background seismicity rate
# K               - productivity parameter of the ETAS triggering kernel
# c0              - Omori-Utsu offset parameter; avoids singularity at zero
# p0              - Omori-Utsu temporal decay parameter
# alpha           - magnitude productivity parameter
# Mr              - reference magnitude for productivity
# t0              - initial time
# m_main          - optional magnitude of an initial mainshock at time t0
# G               - magnitude offset in the completeness model
# H               - temporal decay parameter in the completeness model
# Mcat            - background/network completeness magnitude
# sigma           - smoothness parameter of the detection probability
# log_base        - base of the logarithm used in Mc(t)
# return_complete - if TRUE, returns both complete and observed catalogs
#
# If m_main is NULL, the catalog starts with an empty history.
# If m_main is supplied, the simulation is conditional on an initial
# mainshock at time t0. This mainshock is used in the ETAS history and in
# the completeness calculation, but it is not included in the returned
# observed catalog unless it is explicitly added elsewhere.
#
# Output:
# If return_complete = FALSE:
#   A data frame containing only the observed events after filtering.
#
# If return_complete = TRUE:
#   A list with:
#   complete    - complete simulated catalog with detection information
#   observed    - filtered observed catalog
#   n_complete  - number of events in the complete catalog
#   n_observed  - number of retained observed events
simulate_catalog_level2 <- function(n, b, mc, bin = 0.1,
                                    mu = 0.5, K = 0.05, c0 = 0.01,
                                    p0 = 1.1, alpha = 1.0, Mr = mc,
                                    t0 = 0, m_main = NULL,
                                    G = 5.6, H = 1,
                                    Mcat = 0.95, sigma = 0.2,
                                    log_base = 10,
                                    return_complete = FALSE) {
  
  cat_complete <- simulate_catalog_level1(
    n = n,
    b = b,
    mc = mc,
    bin = bin,
    mu = mu,
    K = K,
    c0 = c0,
    p0 = p0,
    alpha = alpha,
    Mr = Mr,
    t0 = t0,
    m_main = m_main
  )
  
  mc_t <- compute_mc_time(
    cat = cat_complete,
    G = G,
    H = H,
    Mcat = Mcat,
    t0 = t0,
    m_main = m_main,
    log_base = log_base
  )
  
  p_det <- pnorm((cat_complete$magnitude - mc_t) / sigma)
  observed <- runif(nrow(cat_complete)) <= p_det
  
  cat_complete$Mc_time <- mc_t
  cat_complete$p_det <- p_det
  cat_complete$observed <- observed
  
  cat_observed <- cat_complete[observed, ]
  rownames(cat_observed) <- NULL
  cat_observed$id_obs <- seq_len(nrow(cat_observed))
  if (return_complete) {
    return(list(
      complete = cat_complete,
      observed = cat_observed,
      n_complete = nrow(cat_complete),
      n_observed = nrow(cat_observed)
    ))
  }
  cat_observed
}
# The detection threshold is computed from the complete event history,
# including events that will subsequently be removed by thinning.
# An optional initial mainshock is also included in this history.
compute_mc_time <- function(cat, G = 5.6, H = 1,
                            Mcat = 0.95,
                            t0 = 0, m_main = NULL,
                            log_base = 10,
                            eps = 1e-12) {
  
  n <- nrow(cat)
  mc_t <- numeric(n)
  
  for (i in seq_len(n)) {
    ti <- cat$time[i]
    # Previous complete events
    if (i == 1) {
      t_prev <- numeric(0)
      m_prev <- numeric(0)
    } else {
      t_prev <- cat$time[seq_len(i - 1)]
      m_prev <- cat$magnitude[seq_len(i - 1)]
    }
    # Optional initial mainshock
    if (!is.null(m_main)) {
      t_prev <- c(t0, t_prev)
      m_prev <- c(m_main, m_prev)
    }
    if (length(t_prev) == 0) {
      mc_t[i] <- Mcat
    } else {
      dt <- pmax(ti - t_prev, eps)
      mc_candidates <- m_prev - G - H * log(dt, base = log_base)
      # Optional upper bound: an event cannot obscure events larger than itself.
      mc_candidates <- pmin(mc_candidates, m_prev)
      mc_t[i] <- max(Mcat, mc_candidates)
    }
  }
  mc_t
}
# Apply detection filter
apply_detection_filter <- function(cat, G = 5.6, H = 1,
                                   Mcat = 0.95, sigma = 0.2,
                                   t0 = 0, m_main = NULL,
                                   log_base = 10) {
  
  mc_t <- compute_mc_time(
    cat = cat,
    G = G,
    H = H,
    Mcat = Mcat,
    t0 = t0,
    m_main = m_main,
    log_base = log_base
  )
  
  p_det <- pnorm((cat$magnitude - mc_t) / sigma)
  observed <- runif(nrow(cat)) <= p_det
  cat$Mc_time <- mc_t
  cat$p_det <- p_det
  cat$observed <- observed
  cat_obs <- cat[observed, ]
  rownames(cat_obs) <- NULL
  cat_obs$id_obs <- seq_len(nrow(cat_obs))
  cat_obs
}