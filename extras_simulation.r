# -------------------------------------------------------------------------
# Monte Carlo evaluation of b-value estimators without jittering
# -------------------------------------------------------------------------
# Evaluates 24 estimation and confidence-interval strategies using
# magnitudes generated from a discretized Gutenberg-Richter model.
#
# The procedures include continuous and discrete magnitude-based
# estimators, absolute-difference estimators, b-positive, and
# b-more-positive methods. Confidence intervals are evaluated using
# exact, model-based, and nonparametric constructions, as applicable.
#
# Each replicate generates a seismic catalog using sim_catalog().
# Level 1 produces a complete catalog, whereas Level 2 applies
# time-dependent incompleteness, resulting in a random observed size.
#
# Invalid catalogs are regenerated, up to a maximum of 20 attempts.
# All methods are evaluated on the same accepted catalog.
#
# Input:
# m     - number of Monte Carlo replicates
# a     - significance level for confidence intervals
# lev   - catalog simulation level (1 or 2)
# n     - number of events in the complete simulated catalog
# b     - true Gutenberg-Richter b-value
# ell   - maximum forward lag for b-more-positive
# mc    - lower magnitude threshold
# dc    - minimum retained positive-difference threshold
# bin   - magnitude bin width
#
# Additional arguments specify the ETAS and detection-model parameters.
#
# Output:
# p   - empirical confidence-interval coverage probabilities (%)
# pae - mean absolute relative estimation errors (%)
# md  - mean confidence-interval widths
# n0  - observed catalog sizes in the accepted replicates
# -------------------------------------------------------------------------
run_simulation <- function(m, a, lev, n, b, ell, mc = 0, dc = 0.1, bin = 0.1,
                           mu = 0.5, K = 0.05, c0 = 0.01,
                           p0 = 1.1, alpha = 1.0, Mr = mc,
                           t0 = 0, m_main = NULL,
                           G = 5.6, H = 1,
                           Mcat = 0.95, sigma = 0.2) {
  max_attempts <- 20
  cont <- rep(0, 24)
  ae <- rep(0, 24)
  d <- rep(0, 24)
  n0 <- rep(0, m)
  attempts_used <- rep(0, m)
  for (j in 1:m) {
    if (j == 1 || j %% 500 == 0) {
      cat("Iteración, j =", j, "\n")
    }
    valid <- FALSE
    cont2 <- 0
    # Initialize in case all attempts fail
    nr <- NA_integer_
    reason <- NULL
    while (!valid && cont2 < max_attempts) {
      cont2 <- cont2 + 1
      # simulating a seismic catalog (t_i, m_i)
      rsim <- sim_catalog(lev, n, b, mc = mc, bin = bin,
                          mu = mu, K = K, c0 = c0,
                          p0 = p0, alpha = alpha, Mr = Mr,
                          t0 = t0, m_main = m_main,
                          G = G, H = H,
                          Mcat = Mcat, sigma = sigma)
      tt <- rsim$time
      mags <- rsim$magnitude
      nr <- length(mags)
      ######################################
      # Calculating basic stats internally #
      ######################################
      # Basic stats: mags base, abs-dif, b-pos
      basic_stats <- basic_magnitde_statistics(tt, mags, dc)
      if (!basic_stats$valid) {
        valid <- FALSE
        reason <- basic_stats$reason
        
      } else {
        # Basic stats: b-more-positive
        basic_stats_bmp <- basic_stats_bmore_positive(basic_stats$magso, ell, dc)
        if (!basic_stats_bmp$valid) {
          valid <- FALSE
          reason <- basic_stats_bmp$reason
        } else {
          valid <- TRUE
          reason <- NULL
        }
      }
    }
    
    if (!valid) {
      stop(paste0(
        "Iteration ", j,
        " failed after ", max_attempts,
        " attempts. Last observed sample size was nr = ", nr, ". ",
        "Invalid basic statistics: ", reason
      ))
    }
    ###############
    # Basic stats #
    ###############
    # magnitudes based
    nr <- basic_stats$nr 
    mbar <- basic_stats$mbar 
    msd <- basic_stats$msd
    # magnitudes difference based
    N <- basic_stats$N 
    adbar <- basic_stats$adbar 
    s2_ad <- basic_stats$s2_ad
    # b-positive
    r <- basic_stats$r 
    bbar <- basic_stats$bbar 
    bsd <- basic_stats$bsd
    # b-more-positive
    rl <- basic_stats_bmp$rl 
    bbarl <- basic_stats_bmp$bbarl 
    bsdl <- basic_stats_bmp$bsdl
    #
    n0[j] <- nr
    attempts_used[j] <- cont2
    ###############
    # Estimation  #
    ###############
    # Magnitude based estimators
    b_aki_val <- b_aki(nr, mbar, msd, mc, a = a)
    b_bender_val <- b_bender_tgeom(nr, mbar, msd, mc, bin, mmax = 10, a = a, tol = 1e-8)
    b_gut_hop_val <- b_gut_hop(nr, mbar, msd, mc, bin, a = a)
    b_marz_val <- b_marzocchi(nr, mbar, msd, mc, bin, a = a)
    # Difference magnitudes based estimators
    b_dif_val <- b_aki(N, adbar, sqrt(s2_ad), 0, a = a)
    b_dif_disc_val <- b_dif_disc(N, adbar, s2_ad, bin, a)
    b_dif_disc_exact_val <- b_dif_disc_exact(N, adbar, bin, a)
    # b-positive
    b_pos_cont_val <- b_aki(r, bbar, bsd, dc, a = a)
    b_pos_disc_val <- b_gut_hop(r, bbar, bsd, dc, bin, a = a)
    # b-more-positive
    bm_pos_cont_val <- b_aki(rl, bbarl, bsdl, dc, a = a)
    ##############
    # Evaluation #
    ##############
    # Aku-Utsu asymptotic mb
    res1 <- eval_strategy(b, b_aki_val$bn , b_aki_val$La_model, b_aki_val$Ua_model)
    cont[1] <- res1$c + cont[1]
    ae[1] <- res1$s + ae[1]
    d[1] <- res1$d + d[1]
    # Aku-Utsu asymptotic np -shi and bolt
    res2 <- eval_strategy(b, b_aki_val$bn , b_aki_val$La_np, b_aki_val$Ua_np)
    cont[2] <- res2$c + cont[2]
    ae[2] <- res2$s + ae[2]
    d[2] <- res2$d + d[2]
    # Aku-Utsu exact
    res3 <- eval_strategy(b, b_aki_val$bn , b_aki_val$Le, b_aki_val$Ue)
    cont[3] <- res3$c + cont[3]
    ae[3] <- res3$s + ae[3]
    d[3] <- res3$d + d[3]
    # Bender model
    res4 <- eval_strategy(b, b_bender_val$bn , b_bender_val$La_model, b_bender_val$Ua_model)
    cont[4] <- res4$c + cont[4]
    ae[4] <- res4$s + ae[4]
    d[4] <- res4$d + d[4]
    # Bender np
    res5 <- eval_strategy(b, b_bender_val$bn , b_bender_val$La_np, b_bender_val$Ua_np)
    cont[5] <- res5$c + cont[5]
    ae[5] <- res5$s + ae[5]
    d[5] <- res5$d + d[5]
    # Guttorp-Hopkins asymptotic model
    res6 <- eval_strategy(b, b_gut_hop_val$bn , b_gut_hop_val$La_model, b_gut_hop_val$Ua_model)
    cont[6] <- res6$c + cont[6]
    ae[6] <- res6$s + ae[6]
    d[6] <- res6$d + d[6]
    # Guttorp-Hopkins asymptotic np
    res7 <- eval_strategy(b, b_gut_hop_val$bn , b_gut_hop_val$La_np, b_gut_hop_val$Ua_np)
    cont[7] <- res7$c + cont[7]
    ae[7] <- res7$s + ae[7]
    d[7] <- res7$d + d[7]
    # Guttorp-Hopkins exact -> Tinti and Mulargia
    res8 <- eval_strategy(b, b_gut_hop_val$bn , b_gut_hop_val$Le, b_gut_hop_val$Ue)
    cont[8] <- res8$c + cont[8]
    ae[8] <- res8$s + ae[8]
    d[8] <- res8$d + d[8]
    # Marzocchi plug-in
    res9 <- eval_strategy(b, b_marz_val$bn, b_marz_val$L_utsu_plugin, b_marz_val$U_utsu_plugin)
    cont[9] <- res9$c + cont[9]
    ae[9] <- res9$s + ae[9]
    d[9] <- res9$d + d[9]
    # b diff cont asymptotic model
    res10 <- eval_strategy(b, b_dif_val$bn, b_dif_val$La_model, b_dif_val$Ua_model)
    cont[10] <- res10$c + cont[10]
    ae[10] <- res10$s + ae[10]
    d[10] <- res10$d + d[10]
    # b diff cont asymptotic np
    res11 <- eval_strategy(b, b_dif_val$bn, b_dif_val$La_np, b_dif_val$Ua_np)
    cont[11] <- res11$c + cont[11]
    ae[11] <- res11$s + ae[11]
    d[11] <- res11$d + d[11]
    # b diff cont exact
    res12 <- eval_strategy(b, b_dif_val$bn, b_dif_val$Le, b_dif_val$Ue)
    cont[12] <- res12$c + cont[12]
    ae[12] <- res12$s + ae[12]
    d[12] <- res12$d + d[12]
    # b diff disc asymptotic model
    res13 <- eval_strategy(b, b_dif_disc_val$bn, b_dif_disc_val$La_model, b_dif_disc_val$Ua_model)
    cont[13] <- res13$c + cont[13]
    ae[13] <- res13$s + ae[13]
    d[13] <- res13$d + d[13]
    # b diff disc asymptotic np
    res14 <- eval_strategy(b, b_dif_disc_val$bn, b_dif_disc_val$La_np, b_dif_disc_val$Ua_np)
    cont[14] <- res14$c + cont[14]
    ae[14] <- res14$s + ae[14]
    d[14] <- res14$d + d[14]
    # b diff disc exact
    res15 <- eval_strategy(b, b_dif_disc_exact_val$bn, b_dif_disc_exact_val$Le, b_dif_disc_exact_val$Ue)
    cont[15] <- res15$c + cont[15]
    ae[15] <- res15$s + ae[15]
    d[15] <- res15$d + d[15]
    # b+ cont asymptotic mb
    res16 <- eval_strategy(b, b_pos_cont_val$bn, b_pos_cont_val$La_model, b_pos_cont_val$Ua_model)
    cont[16] <- res16$c + cont[16]
    ae[16] <- res16$s + ae[16]
    d[16] <- res16$d + d[16]
    # b+ cont asymptotic np
    res17 <- eval_strategy(b, b_pos_cont_val$bn, b_pos_cont_val$La_np, b_pos_cont_val$Ua_np)
    cont[17] <- res17$c + cont[17]
    ae[17] <- res17$s + ae[17]
    d[17] <- res17$d + d[17]
    # b+ cont exact
    res18 <- eval_strategy(b, b_pos_cont_val$bn, b_pos_cont_val$Le, b_pos_cont_val$Ue)
    cont[18] <- res18$c + cont[18]
    ae[18] <- res18$s + ae[18]
    d[18] <- res18$d + d[18]
    # b+ disc asymptotic mb
    res19 <- eval_strategy(b, b_pos_disc_val$bn, b_pos_disc_val$La_model, b_pos_disc_val$Ua_model)
    cont[19] <- res19$c + cont[19]
    ae[19] <- res19$s + ae[19]
    d[19] <- res19$d + d[19]
    # b+ disc asymptotic np
    res20 <- eval_strategy(b, b_pos_disc_val$bn, b_pos_disc_val$La_np, b_pos_disc_val$Ua_np)
    cont[20] <- res20$c + cont[20]
    ae[20] <- res20$s + ae[20]
    d[20] <- res20$d + d[20]
    # b+ disc exact
    res21 <- eval_strategy(b, b_pos_disc_val$bn, b_pos_disc_val$Le, b_pos_disc_val$Ue)
    cont[21] <- res21$c + cont[21]
    ae[21] <- res21$s + ae[21]
    d[21] <- res21$d + d[21]
    # b-more-positive cont asymptotic mb
    res22 <- eval_strategy(b, bm_pos_cont_val$bn, bm_pos_cont_val$La_model, bm_pos_cont_val$Ua_model)
    cont[22] <- res22$c + cont[22]
    ae[22] <- res22$s + ae[22]
    d[22] <- res22$d + d[22]
    # b-more-positive cont asymptotic np
    res23 <- eval_strategy(b, bm_pos_cont_val$bn, bm_pos_cont_val$La_np, bm_pos_cont_val$Ua_np)
    cont[23] <- res23$c + cont[23]
    ae[23] <- res23$s + ae[23]
    d[23] <- res23$d + d[23]
    # b-more-positive cont exact
    res24 <- eval_strategy(b, bm_pos_cont_val$bn, bm_pos_cont_val$Le, bm_pos_cont_val$Ue)
    cont[24] <- res24$c + cont[24]
    ae[24] <- res24$s + ae[24]
    d[24] <- res24$d + d[24]
  }
  p <- 100*cont/m 
  pae <- 100*ae/m
  md <- d/m
  list(p = p, pae = pae, md = md, n0 = n0)
}
# -------------------------------------------------------------------------
# Monte Carlo evaluation of continuous estimators with uniform jittering
# -------------------------------------------------------------------------
# Evaluates 12 confidence-interval strategies based on continuous
# b-value estimators applied to jittered magnitudes.
#
# Each catalog is generated from the discretized Gutenberg-Richter
# model and, for Level 2, subjected to time-dependent incompleteness.
# Uniform jittering is then applied:
#
#   M_i^U = M_i' + bin * U_i,  U_i ~ Uniform(0,1).
#
# A single jittered catalog is shared by all continuous estimators
# within each replicate.
#
# The evaluated procedures are Aki-Utsu, continuous absolute
# differences, b-positive, and b-more-positive, each with exact,
# model-based, and nonparametric interval constructions.
#
# Invalid catalogs are regenerated, up to a maximum of 20 attempts.
#
# Input:
# m, a, lev, n, b, ell, mc, dc, bin - as in run_simulation()
# Additional arguments specify the ETAS and detection-model parameters.
#
# Output:
# p   - empirical confidence-interval coverage probabilities (%)
# pae - mean absolute relative estimation errors (%)
# md  - mean confidence-interval widths
# n0  - observed catalog sizes in the accepted replicates
# -------------------------------------------------------------------------
run_simulation_cont_jitter <- function(m, a, lev, n, b, ell, mc = 0, dc = 0.1, bin = 0.1,
                                mu = 0.5, K = 0.05, c0 = 0.01,
                                p0 = 1.1, alpha = 1.0, Mr = mc,
                                t0 = 0, m_main = NULL,
                                G = 5.6, H = 1, Mcat = 0.95, sigma = 0.2) {
  max_attempts <- 20
  cont <- rep(0, 12)
  ae <- rep(0, 12)
  d <- rep(0, 12)
  n0 <- rep(0, m)
  attempts_used <- rep(0, m)
  for (j in 1:m) {
    if (j == 1 || j %% 500 == 0) {
      cat("Iteración, j =", j, "\n")
    }
    valid <- FALSE
    cont2 <- 0
    # Initialize in case all attempts fail
    nr <- NA_integer_
    reason <- NULL
    while (!valid && cont2 < max_attempts) {
      cont2 <- cont2 + 1
      # simulating a seismic catalog (t_i, m_i)
      rsim <- sim_catalog(lev, n, b, mc = mc, bin = bin,
                          mu = mu, K = K, c0 = c0,
                          p0 = p0, alpha = alpha, Mr = Mr,
                          t0 = t0, m_main = m_main,
                          G = G, H = H,
                          Mcat = Mcat, sigma = sigma)
      tt <- rsim$time
      nr <- length(rsim$magnitude)
      u <- runif(nr)
      mags <- rsim$magnitude + bin*u
      ######################################
      # Calculating basic stats internally #
      ######################################
      # Basic stats: mags base, abs-dif, b-pos
      basic_stats <- basic_magnitde_statistics(tt, mags, dc)
      if (!basic_stats$valid) {
        valid <- FALSE
        reason <- basic_stats$reason
      } else {
        # Basic stats: b-more-positive
        basic_stats_bmp <- basic_stats_bmore_positive(basic_stats$magso, ell, dc)
        if (!basic_stats_bmp$valid) {
          valid <- FALSE
          reason <- basic_stats_bmp$reason
        } else {
          valid <- TRUE
          reason <- NULL
        }
      }
    }
    
    if (!valid) {
      stop(paste0(
        "Iteration ", j,
        " failed after ", max_attempts,
        " attempts. Last observed sample size was nr = ", nr, ". ",
        "Invalid basic statistics: ", reason
      ))
    }
    ###############
    # Basic stats #
    ###############
    # magnitudes based
    nr <- basic_stats$nr 
    mbar <- basic_stats$mbar 
    msd <- basic_stats$msd
    # magnitudes difference based
    N <- basic_stats$N 
    adbar <- basic_stats$adbar 
    s2_ad <- basic_stats$s2_ad
    # b-positive
    r <- basic_stats$r 
    bbar <- basic_stats$bbar 
    bsd <- basic_stats$bsd
    # b-more-positive
    rl <- basic_stats_bmp$rl 
    bbarl <- basic_stats_bmp$bbarl 
    bsdl <- basic_stats_bmp$bsdl
    #
    n0[j] <- nr
    attempts_used[j] <- cont2
    ###############
    # Estimation  #
    ###############
    # Magnitude based estimators
    b_aki_val <- b_aki(nr, mbar, msd, mc, a = a)
    # Difference magnitudes based estimators
    b_dif_val <- b_aki(N, adbar, sqrt(s2_ad), 0, a = a)
    # b-positive
    b_pos_cont_val <- b_aki(r, bbar, bsd, dc, a = a)
    # b-more-positive
    bm_pos_cont_val <- b_aki(rl, bbarl, bsdl, dc, a = a)
    ##############
    # Evaluation #
    ##############
    # Aku-Utsu asymptotic mb
    res1 <- eval_strategy(b, b_aki_val$bn , b_aki_val$La_model, b_aki_val$Ua_model)
    cont[1] <- res1$c + cont[1]
    ae[1] <- res1$s + ae[1]
    d[1] <- res1$d + d[1]
    # Aku-Utsu asymptotic np -shi and bolt
    res2 <- eval_strategy(b, b_aki_val$bn , b_aki_val$La_np, b_aki_val$Ua_np)
    cont[2] <- res2$c + cont[2]
    ae[2] <- res2$s + ae[2]
    d[2] <- res2$d + d[2]
    # Aku-Utsu exact
    res3 <- eval_strategy(b, b_aki_val$bn , b_aki_val$Le, b_aki_val$Ue)
    cont[3] <- res3$c + cont[3]
    ae[3] <- res3$s + ae[3]
    d[3] <- res3$d + d[3]
    # b diff cont asymptotic model
    res10 <- eval_strategy(b, b_dif_val$bn, b_dif_val$La_model, b_dif_val$Ua_model)
    cont[4] <- res10$c + cont[4]
    ae[4] <- res10$s + ae[4]
    d[4] <- res10$d + d[4]
    # b diff cont asymptotic np
    res11 <- eval_strategy(b, b_dif_val$bn, b_dif_val$La_np, b_dif_val$Ua_np)
    cont[5] <- res11$c + cont[5]
    ae[5] <- res11$s + ae[5]
    d[5] <- res11$d + d[5]
    # b diff cont exact
    res12 <- eval_strategy(b, b_dif_val$bn, b_dif_val$Le, b_dif_val$Ue)
    cont[6] <- res12$c + cont[6]
    ae[6] <- res12$s + ae[6]
    d[6] <- res12$d + d[6]
    # b+ cont asymptotic mb
    res16 <- eval_strategy(b, b_pos_cont_val$bn, b_pos_cont_val$La_model, b_pos_cont_val$Ua_model)
    cont[7] <- res16$c + cont[7]
    ae[7] <- res16$s + ae[7]
    d[7] <- res16$d + d[7]
    # b+ cont asymptotic np
    res17 <- eval_strategy(b, b_pos_cont_val$bn, b_pos_cont_val$La_np, b_pos_cont_val$Ua_np)
    cont[8] <- res17$c + cont[8]
    ae[8] <- res17$s + ae[8]
    d[8] <- res17$d + d[8]
    # b+ cont exact
    res18 <- eval_strategy(b, b_pos_cont_val$bn, b_pos_cont_val$Le, b_pos_cont_val$Ue)
    cont[9] <- res18$c + cont[9]
    ae[9] <- res18$s + ae[9]
    d[9] <- res18$d + d[9]
    # b-more-positive cont asymptotic mb
    res22 <- eval_strategy(b, bm_pos_cont_val$bn, bm_pos_cont_val$La_model, bm_pos_cont_val$Ua_model)
    cont[10] <- res22$c + cont[10]
    ae[10] <- res22$s + ae[10]
    d[10] <- res22$d + d[10]
    # b-more-positive cont asymptotic np
    res23 <- eval_strategy(b, bm_pos_cont_val$bn, bm_pos_cont_val$La_np, bm_pos_cont_val$Ua_np)
    cont[11] <- res23$c + cont[11]
    ae[11] <- res23$s + ae[11]
    d[11] <- res23$d + d[11]
    # b-more-positive cont exact
    res24 <- eval_strategy(b, bm_pos_cont_val$bn, bm_pos_cont_val$Le, bm_pos_cont_val$Ue)
    cont[12] <- res24$c + cont[12]
    ae[12] <- res24$s + ae[12]
    d[12] <- res24$d + d[12]
  }
  p <- 100*cont/m 
  pae <- 100*ae/m
  md <- d/m
  list(p = p, pae = pae, md = md, n0 = n0)
}
# -------------------------------------------------------------------------
# Monte Carlo evaluation of magnitude-based estimators with estimated mc
# -------------------------------------------------------------------------
# Evaluates nine confidence-interval strategies for magnitude-based
# b-value estimators when the completeness threshold is estimated
# separately for each simulated catalog.
#
# The completeness magnitude is estimated from the original binned
# magnitudes using the maximum curvature method with a correction
# of 0.2 magnitude units.
#
# Events below the estimated threshold are excluded.
# The discrete estimators are applied to the original binned
# magnitudes, whereas the continuous Aki-Utsu estimator is applied
# after uniform jittering.
#
# The evaluated estimators are Aki-Utsu, Bender, Guttorp-Hopkins,
# and Marzocchi, with their corresponding confidence intervals.
#
# Invalid catalogs are regenerated, up to a maximum of 20 attempts.
#
# Input:
# m, a, lev, n, b, ell, mc, dc, bin - as in run_simulation()
# Additional arguments specify the ETAS and detection-model parameters.
#
# Output:
# p      - empirical confidence-interval coverage probabilities (%)
# pae    - mean absolute relative estimation errors (%)
# md     - mean confidence-interval widths
# n0     - sample sizes after filtering at the estimated mc
# mc_est - estimated completeness magnitudes across replicates
# -------------------------------------------------------------------------
run_simulation_magnitude_based <- function(m, a, lev, n, b, ell, mc = 0, dc = 0.1, bin = 0.1,
                           mu = 0.5, K = 0.05, c0 = 0.01,
                           p0 = 1.1, alpha = 1.0, Mr = mc,
                           t0 = 0, m_main = NULL,
                           G = 5.6, H = 1,
                           Mcat = 0.95, sigma = 0.2) {
  max_attempts <- 20
  cont <- rep(0, 9)
  ae <- rep(0, 9)
  d <- rep(0, 9)
  n0 <- rep(0, m)
  mc_est <- rep(0, m)
  attempts_used <- rep(0, m)
  for (j in 1:m) {
    if (j == 1 || j %% 500 == 0) {
      cat("Iteración, j =", j, "\n")
    }
    valid <- FALSE
    cont2 <- 0
    # Initialize in case all attempts fail
    nr <- NA_integer_
    reason <- NULL
    while (!valid && cont2 < max_attempts) {
      cont2 <- cont2 + 1
      # simulating a seismic catalog (t_i, m_i)
      rsim <- sim_catalog(lev, n, b, mc = mc, bin = bin,
                          mu = mu, K = K, c0 = c0,
                          p0 = p0, alpha = alpha, Mr = Mr,
                          t0 = t0, m_main = m_main,
                          G = G, H = H,
                          Mcat = Mcat, sigma = sigma)
      mc_hat <- mc_max_curvature(rsim$magnitude, bin = bin, correction = 0.2)
      id <- rsim$magnitude >= mc_hat
      tt <- rsim$time[id]
      mags <- rsim$magnitude[id]
      nn <- length(rsim$magnitude)
      ######################################
      # Calculating basic stats internally #
      ######################################
      # Basic stats: mags base, abs-dif, b-pos
      # Validate magnitude-based statistics only
      nr <- length(mags)
      mbar <- if (nr > 0L) mean(mags) else NA_real_
      msd <- if (nr > 1L) sd(mags) else NA_real_
      
      valid <- nr >= 8L &&
        all(is.finite(mags)) &&
        is.finite(mbar) &&
        is.finite(msd) &&
        msd > 0 &&
        mbar > mc_hat
      
      reason <- if (valid) NULL else "invalid magnitude statistics"
    }
    if (!valid) {
      stop(paste0(
        "Iteration ", j,
        " failed after ", max_attempts,
        " attempts. Last observed sample size was nr = ", nn, ". ",
        "Invalid basic statistics: ", reason
      ))
    }
    ###############
    # Basic stats #
    ###############
    n0[j] <- nr
    mc_est[j] <- mc_hat
    attempts_used[j] <- cont2
    ###############
    # Estimation  #
    ###############
    u <- runif(length(rsim$magnitude))
    mags2 <- rsim$magnitude + bin*u
    id2 <- mags2 >= mc_hat
    mags2 <- mags2[id2]
    nr2 <- length(mags2)
    # Magnitude based estimators
    b_aki_val <- b_aki(nr2, mean(mags2), sd(mags2), mc_hat, a = a)
    b_bender_val <- b_bender_tgeom(nr, mbar, msd, mc_hat, bin, mmax = 10, a = a, tol = 1e-8)
    b_gut_hop_val <- b_gut_hop(nr, mbar, msd, mc_hat, bin, a = a)
    b_marz_val <- b_marzocchi(nr, mbar, msd, mc_hat, bin, a = a)
    ##############
    # Evaluation #
    ##############
    # Aku-Utsu asymptotic mb
    res1 <- eval_strategy(b, b_aki_val$bn , b_aki_val$La_model, b_aki_val$Ua_model)
    cont[1] <- res1$c + cont[1]
    ae[1] <- res1$s + ae[1]
    d[1] <- res1$d + d[1]
    # Aku-Utsu asymptotic np -shi and bolt
    res2 <- eval_strategy(b, b_aki_val$bn , b_aki_val$La_np, b_aki_val$Ua_np)
    cont[2] <- res2$c + cont[2]
    ae[2] <- res2$s + ae[2]
    d[2] <- res2$d + d[2]
    # Aku-Utsu exact
    res3 <- eval_strategy(b, b_aki_val$bn , b_aki_val$Le, b_aki_val$Ue)
    cont[3] <- res3$c + cont[3]
    ae[3] <- res3$s + ae[3]
    d[3] <- res3$d + d[3]
    # Bender model
    res4 <- eval_strategy(b, b_bender_val$bn , b_bender_val$La_model, b_bender_val$Ua_model)
    cont[4] <- res4$c + cont[4]
    ae[4] <- res4$s + ae[4]
    d[4] <- res4$d + d[4]
    # Bender np
    res5 <- eval_strategy(b, b_bender_val$bn , b_bender_val$La_np, b_bender_val$Ua_np)
    cont[5] <- res5$c + cont[5]
    ae[5] <- res5$s + ae[5]
    d[5] <- res5$d + d[5]
    # Guttorp-Hopkins asymptotic model
    res6 <- eval_strategy(b, b_gut_hop_val$bn , b_gut_hop_val$La_model, b_gut_hop_val$Ua_model)
    cont[6] <- res6$c + cont[6]
    ae[6] <- res6$s + ae[6]
    d[6] <- res6$d + d[6]
    # Guttorp-Hopkins asymptotic np
    res7 <- eval_strategy(b, b_gut_hop_val$bn , b_gut_hop_val$La_np, b_gut_hop_val$Ua_np)
    cont[7] <- res7$c + cont[7]
    ae[7] <- res7$s + ae[7]
    d[7] <- res7$d + d[7]
    # Guttorp-Hopkins exact -> Tinti and Mulargia
    res8 <- eval_strategy(b, b_gut_hop_val$bn , b_gut_hop_val$Le, b_gut_hop_val$Ue)
    cont[8] <- res8$c + cont[8]
    ae[8] <- res8$s + ae[8]
    d[8] <- res8$d + d[8]
    # Marzocchi plug-in
    res9 <- eval_strategy(b, b_marz_val$bn, b_marz_val$L_utsu_plugin, b_marz_val$U_utsu_plugin)
    cont[9] <- res9$c + cont[9]
    ae[9] <- res9$s + ae[9]
    d[9] <- res9$d + d[9]
  }
  p <- 100*cont/m 
  pae <- 100*ae/m
  md <- d/m
  list(p = p, pae = pae, md = md, 
       n0 = n0, mc_est = mc_est)
}
# -------------------------------------------------------------------------
# Performance measures for a single estimation strategy
# -------------------------------------------------------------------------
# Evaluates a point estimator and its confidence interval against
# the true Gutenberg-Richter b-value.
#
# Input:
# b    - true b-value
# bhat - estimated b-value
# L    - lower confidence-interval endpoint
# U    - upper confidence-interval endpoint
#
# Output:
# c - coverage indicator: 1 if L <= b <= U, 0 otherwise
# s - absolute relative estimation error: |b - bhat| / b
# d - confidence-interval width: U - L
#
# The Monte Carlo routines average these measures over replicates.
# -------------------------------------------------------------------------
eval_strategy <- function(b, bhat, L, U){
  cc <- 0
  if(b >= L && U >= b){
    cc <- 1
  }
  ss <- abs(b - bhat)/b
  dd <- U - L
  list(c = cc, s = ss, d = dd)
}