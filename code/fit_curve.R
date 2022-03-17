# Fit a curve to the d18O values
fit_curve <- function(isodata) {
  
  # due to a bug in nls, d must be in the global environment
  d <<- data.frame(X = -isodata$measure, Y = isodata$d18O)
  
  # # Exploration of plausible parameter ranges
  # source("code/simple_iso_plot.R")
  # FD1 <- function(x, A, x_0, z, M) {
  #   A * cos(2 * pi * ((x - x_0) / z)) + M
  # }
  # x <- seq(50,0,0.5)
  # # A: Amplitude
  # # x_0: phase
  # # z: period
  # # M: y-axis offset
  # 
  # FD <- function(x) { FD1(
  #   x,
  #   #1,0,10,-15 # low_1
  #   #5,10,45,-5 # start_1
  #   #5,30,45,-5 # start_2
  #   #3,22,16,-8 # start3
  #   #15,45,90,10 #up_1
  #   fit$m$getPars()[["A"]], fit$m$getPars()[["x_0"]], fit$m$getPars()[["z"]], fit$m$getPars()[["M"]]
  # )}
  # simple_iso_plot(d$X, d$Y, isodata$d13C, FD)
  
  # cos curve fitting starting parameters and search range
  low     <-  list(A = 1,  x_0 = 0,  z = 10, M = -15)
  starts <- list(
    start_1 = list(A = 3,  x_0 = 22, z = 16, M = -8 ),
    start_2 = list(A = 5,  x_0 = 10, z = 45, M = -5 ),
    start_3 = list(A = 5,  x_0 = 30, z = 45, M = -5 )
  )
  up      <-  list(A = 15, x_0 = 45, z = 90, M =  10)
  
  # fitting attempts
  potential_fits <- purrr::map(starts, function(start) {
    try(
      stats::nls(
        Y ~ A * cos(2 * pi * ((X - x_0) / z)) + M,
        data = d,
        lower = low, start = start, upper = up,
        algorithm = "port",
        control = nls.control(maxiter = 100000)
      ),
      silent = T
    )
  })
  # remove failed fits
  working_fits <- purrr::discard(potential_fits, function(x) {class(x) == "try-error"})
  # select best fit
  fit <- working_fits[[which.min(purrr::map_dbl(working_fits, function(x) {x$m$deviance()}))]]

  # get error bar for the fitted curve
  # adapted from https://stackoverflow.com/questions/32613119/plot-the-median-confidence-interval-of-a-bootstrap-output-in-ggplot2
  curveBoot <- nlstools::nlsBoot(fit, niter = 10000)
  theta_mat <- curveBoot$coefboot # Matrix with the bootstrapped parameter estimates
  
  # Points where to evaluate the model
  x_eval <- seq(min(d$X), max(d$X), length.out = 1000)
  
  # Matrix with the predictions
  pred_mat <- apply(theta_mat, 1, function(theta) { 
    theta["A"] * cos(2 * pi * ((x_eval - theta["x_0"]) / theta["z"])) + theta["M"]
  })
  
  # Pack the estimates in an output data.frame
  estim_mat <- dplyr::bind_cols(
    x = x_eval,
    tibble::as_tibble(t(apply(pred_mat, 1, function(y_est) {
      c(
        median_est = median(y_est),
        ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE),
        ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
      )
    })))
  )
  
  # prepare output list
  return(list(
    specimen = isodata$specimen[1],
    fit = fit,
    estim_mat = estim_mat,
    theta_mat = tibble::as_tibble(theta_mat) # plausible fittings from nlstools::nlsBoot
  ))

}
