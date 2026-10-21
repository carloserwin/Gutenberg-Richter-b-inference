# -------------------------------------------------------------------------
# Construction of overlapping event-based windows
# -------------------------------------------------------------------------
# Divides a seismic catalog into overlapping windows containing a fixed
# number of consecutive events, ordered by occurrence time.
#
# Each window is represented by the occurrence time of its central event.
# The function also records its time relative to the mainshock and
# whether the window contains the mainshock.
#
# Input:
# BDseq - seismic catalog with a datetime_local column
# pm    - data frame containing the mainshock information
# w     - number of events per window
# step  - number of events between consecutive window starting positions
#
# Output:
# A data frame containing window indices, dates, relative times,
# a mainshock indicator, and a list-column with the corresponding
# event subsets.
# -------------------------------------------------------------------------
make_event_windows <- function(BDseq, pm, w = 400, step = 50){
  
  BDseq <- BDseq[order(BDseq$datetime_local), ]
  rownames(BDseq) <- NULL
  
  n <- nrow(BDseq)
  
  if(n < w){
    stop("The data set has fewer events than the window size.")
  }
  
  starts <- seq(1, n - w + 1, by = step)
  ends <- starts + w - 1
  mids <- floor((starts + ends) / 2)
  
  out <- data.frame(
    window_id = seq_along(starts),
    i_start = starts,
    i_mid = mids,
    i_end = ends,
    n = w,
    t_start = BDseq$datetime_local[starts],
    t_mid = BDseq$datetime_local[mids],
    t_end = BDseq$datetime_local[ends]
  )
  
  # Tiempo usado para graficar:
  # ahora es el tiempo del evento central de la ventana
  out$t_plot <- out$t_mid
  
  # Por si quieres conservar también el criterio anterior
  out$t_plot_end <- out$t_end
  
  out$days_from_mainshock <- as.numeric(
    difftime(out$t_plot, pm$datetime_local[1], units = "days")
  )
  
  out$contains_mainshock <- out$t_start <= pm$datetime_local[1] &
    out$t_end >= pm$datetime_local[1]
  
  windows <- vector("list", length(starts))
  
  for(j in seq_along(starts)){
    windows[[j]] <- BDseq[starts[j]:ends[j], ]
  }
  
  out$window <- I(windows)
  
  out
}
# -------------------------------------------------------------------------
# Temporal evolution of b-value confidence intervals
# -------------------------------------------------------------------------
# Plots the lower and upper confidence limits for two estimators:
# continuous b-positive with jittering and discrete absolute differences.
#
# Confidence intervals are displayed against the central event time
# of each window. The mainshock occurrence time is indicated by
# a vertical dashed line.
#
# Input:
# W400         - data frame containing the event windows
# bp           - b-positive results: lower, bhat, upper, n_eff
# bd           - discrete absolute-difference results, same structure
# pm           - mainshock information
# ylim         - optional vertical axis limits
# y_by         - spacing between vertical axis ticks
# xlab, ylab   - axis labels
# main         - plot title
# month_labels - whether to display monthly labels on the time axis
#
# Output:
# A plot produced using base R graphics.
# -------------------------------------------------------------------------
plot_interval_series <- function(W400, bp, bd, pm,
                                 ylim = NULL,
                                 y_by = 0.2,
                                 xlab = "Months, 2017",
                                 ylab = "95% confidence intervals",
                                 main = "",
                                 month_labels = TRUE) {
  
  tt <- W400$t_plot
  bpL <- bp[, 1]
  bpU <- bp[, 3]
  bdL <- bd[, 1]
  bdU <- bd[, 3]
  if(is.null(ylim)){
    ylim <- range(c(bpL, bpU, bdL, bdU), na.rm = TRUE)
  }
  y_ticks <- seq(ylim[1], ylim[2], by = y_by)
  plot(tt, bpL,
       type = "n",
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       bty = "l")
  
  abline(h = y_ticks,
         col = "gray92",
         lwd = 0.8)
  
  # Eje x con meses en inglés
  if(month_labels) {
    
    tz_use <- attr(tt, "tzone")
    if(is.null(tz_use) || length(tz_use) == 0) {
      tz_use <- "America/Mexico_City"
    }
    
    x0 <- as.POSIXct(format(min(tt), "%Y-%m-01 00:00:00"),
                     tz = tz_use)
    
    x1 <- as.POSIXct(format(max(tt), "%Y-%m-01 00:00:00"),
                     tz = tz_use)
    
    at_month <- seq(x0, x1, by = "1 month")
    
    month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    
    month_labs <- month_names[as.integer(format(at_month, "%m"))]
    
    axis.POSIXct(1,
                 at = at_month,
                 labels = month_labs,
                 las = 2,
                 cex.axis = 0.85,
                 col.axis = "gray25",
                 col.ticks = "gray50")
    
  } else {
    
    axis.POSIXct(1,
                 at = tt,
                 format = "%b",
                 las = 2,
                 cex.axis = 0.85,
                 col.axis = "gray25",
                 col.ticks = "gray50")
  }
  
  axis(2,
       at = y_ticks,
       labels = y_ticks,
       las = 2,
       cex.axis = 0.85,
       col.axis = "gray25",
       col.ticks = "gray50")
  
  # Intervalos
  lines(tt, bpL, lwd = 2.2, col = "darkgray")
  lines(tt, bpU, lwd = 2.2, col = "darkgray")
  
  lines(tt, bdL, lwd = 2.2, col = "black")
  lines(tt, bdU, lwd = 2.2, col = "black")
  
  # Mainshock
  abline(v = pm$datetime_local[1],
         col = "red",
         lty = 2,
         lwd = 1.5)
  
  legend("topleft",
         legend = c("b+, cont + jitter",
                    "Discrete absolute-difference interval",
                    "Mainshock, M 8.2"),
         col = c("darkgray", "black", "red"),
         lty = c(1, 1, 2),
         lwd = c(2.2, 2.2, 1.5),
         bty = "n",
         cex = 0.85)
  
  box(lwd = 2)
}
# -------------------------------------------------------------------------
# Temporal evolution of b-value point estimates
# -------------------------------------------------------------------------
# Plots the b-value estimates obtained from four methods:
# Guttorp-Hopkins, continuous b-positive with jittering,
# continuous b-more-positive with jittering, and discrete
# absolute differences.
#
# Estimates are displayed against the central event time of each
# window. The mainshock occurrence time is indicated by a
# vertical dashed line.
#
# Input:
# W400         - data frame containing the event windows
# bGH          - Guttorp-Hopkins results
# bp           - b-positive results
# bmp          - b-more-positive results
# bd           - discrete absolute-difference results
# pm           - mainshock information
# ylim         - optional vertical axis limits
# y_by         - spacing between vertical axis ticks
# xlab, ylab   - axis labels
# main         - plot title
# month_labels - whether to display monthly labels on the time axis
#
# The result matrices have columns: lower, bhat, upper, n_eff.
#
# Output:
# A plot produced using base R graphics.
# -------------------------------------------------------------------------
plot_point_series <- function(W400, bGH, bp, bmp, bd, pm,
                              ylim = NULL,
                              y_by = 0.2,
                              xlab = "Months, 2017",
                              ylab = "Point estimates for b",
                              main = "",
                              month_labels = TRUE) {
  
  tt <- W400$t_plot
  pbGH <- bGH[,2]
  pbp <- bp[,2]
  pbmp <- bmp[,2]
  pbd <- bd[,2]
  if(is.null(ylim)){
    ylim <- range(c(pbGH, pbp, pbmp, pbd), na.rm = TRUE)
  }
  
  y_ticks <- seq(ylim[1], ylim[2], by = y_by)
  
  plot(tt, pbGH,
       type = "n",
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       bty = "l")
  abline(h = y_ticks, col = "gray92", lwd = 0.8)
  # eje x con meses en inglés
  if(month_labels) {
    tz_use <- attr(tt, "tzone")
    if(is.null(tz_use) || length(tz_use) == 0) {
      tz_use <- "America/Mexico_City"
    }
    
    x0 <- as.POSIXct(format(min(tt), "%Y-%m-01 00:00:00"),
                     tz = tz_use)
    
    x1 <- as.POSIXct(format(max(tt), "%Y-%m-01 00:00:00"),
                     tz = tz_use)
    
    at_month <- seq(x0, x1, by = "1 month")
    
    month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    
    month_labs <- month_names[as.integer(format(at_month, "%m"))]
    
    axis.POSIXct(1,
                 at = at_month,
                 labels = month_labs,
                 las = 2,
                 cex.axis = 0.85,
                 col.axis = "gray25",
                 col.ticks = "gray50")
    
    abline(v = at_month,
           col = "gray95",
           lwd = 0.8)
    
  } else {
    
    axis.POSIXct(1,
                 at = tt,
                 format = "%b",
                 las = 2,
                 cex.axis = 0.85,
                 col.axis = "gray25",
                 col.ticks = "gray50")
  }
  
  axis(2,
       at = y_ticks,
       labels = y_ticks,
       las = 2,
       cex.axis = 0.85,
       col.axis = "gray25",
       col.ticks = "gray50")
  
  # series puntuales
  lines(tt, pbGH,
        lwd = 2.2,
        col = "darkgreen")
  lines(tt, pbp,
        lwd = 2.2,
        col = "darkgray")
  lines(tt, pbmp,
        lwd = 2.2,
        col = "blue")
  lines(tt, pbd,
        lwd = 2.2,
        col = "black")
  
  # mainshock
  abline(v = pm$datetime_local[1],
         col = "red",
         lty = 2,
         lwd = 1.5)
  legend("topleft",
         legend = c("GH",
                    "b+, cont + jitter",
                    "b++, cont + jitter",
                    "Discrete absolute-difference",
                    "Mainshock, M 8.2"),
         col = c("darkgreen", "darkgray", "blue", "black", "red"),
         lty = c(1, 1, 1, 1, 2),
         lwd = c(2.2, 2.2, 2.2, 2.2, 1.5),
         bty = "n",
         cex = 0.85)
  box(lwd = 2)
}