# Fit a curve to the d18O values
fit_curve <- function(isodata) {
  
  # due to a bug in nls, d must be in the global environment
  d <<- data.frame(X = -isodata$measure, Y = isodata$d18O)
  
  # # Exploration of plausible parameter ranges
  # FD1 <- function(x, A, x_0, z, M) {
  #   A * cos(2 * pi * ((x - x_0) / z)) + M
  # }
  # 
  # x <- -50:0
  # 
  # # A: Amplitude
  # # x_0: phase
  # # z: period
  # # M: y-axis offset
  # plot(x, FD1(x, 5, 10, 45, -5))
  # points(d, col = "red")
  # plot(x, FD1(x, 5, -10, 45, -5))
  
  # cos curve fitting starting parameters and search range
  low_1   <- low2 <- list(A = 0,  x_0 = 0,  z = 20, M = -20)
  start_1 <-         list(A = 5,  x_0 = 10, z = 45, M = -5 )
  start_2 <-         list(A = 5,  x_0 = 30, z = 45, M = -5 )
  up_1    <- up2  <- list(A = 15, x_0 = 45, z = 90, M =  10)
  
  # first fitting attempt
  fit1 <- try(
    stats::nls(
      Y ~ A * cos(2 * pi * ((X - x_0) / z)) + M,
      data = d,
      lower = low_1, start = start_1, upper = up_1,
      algorithm = "port",
      control = nls.control(maxiter = 100000)
    ),
    silent = T
  )
  
  # second fitting attempt after moving the starting curve
  fit2 <- try(
      stats::nls(
      Y ~ A * cos(2 * pi * ((X - x_0) / z)) + M,
      data = d,
      lower = low_1, start = start_2, upper = up_1,
      algorithm = "port",
      control = nls.control(maxiter = 100000)
    ),
    silent = T
  )

  fit <- if (class(fit1) == "try-error" & class(fit2) == "try-error") {
    stop("No fitting modell found")
  } else if (class(fit1) != "try-error" & class(fit2) != "try-error") {
    # return the model with the smaller residual sum-of-squares
    if (fit1$m$deviance() <= fit2$m$deviance()) { fit1 } else { fit2 }
  } else if (class(fit1) != "try-error") {
    fit1
  } else {
    fit2
  }

  # get error bar for the fitted curve
  # adapted from https://stackoverflow.com/questions/32613119/plot-the-median-confidence-interval-of-a-bootstrap-output-in-ggplot2
  curveBoot <- nlstools::nlsBoot(fit, niter = 1000)
  theta_mat <- curveBoot$coefboot # Matrix with the bootstrapped parameter estimates
  
  fun <- function(x, theta) {
    theta["A"] * cos(2 * pi * ((x - theta["x_0"]) / theta["z"])) + theta["M"]
  }
  
  # Points where to evaluate the model
  x_eval <- seq(min(d$X), max(d$X), length.out = 100)
  
  # Matrix with the predictions
  pred_mat <- apply(theta_mat, 1, function(theta) { fun(x_eval, theta) })
  
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
    fit = fit,
    estim_mat = estim_mat
  ))

}
