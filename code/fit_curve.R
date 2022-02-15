# Fit a curve to the d18O values
fit_curve <- function(d) {
  
  d <<- d # due to a bug in nls, d must be in the global environment
  
  # fit cos curve model
  fit <- stats::nls(
    Y ~ A * cos(2 * pi * ((X - x_0) / z)) + M,
    data = d,
    lower = list(
      A = -50,
      M = -50,
      x_0 = -50,
      z = -50
    ),
    start = list(
      A = 5,
      M = 5,
      x_0 = -20,
      z = 25
    ),
    upper = list(
      A = 50,
      M = 50,
      x_0 = 50,
      z = 50
    ),
    algorithm = "port",
    control = nls.control(maxiter = 100000),
  ) # If curve is offset from data adjust x_0 variable first then try A

  # fit <- minpack.lm::nlsLM(
  #   formula = Y ~ A * cos(2 * pi * ((X - x_0) / z)) + M,
  #   data = d,
  #   start = list(
  #     A = 5,
  #     M = 5,
  #     x_0 = -20,
  #     z = 25
  #   ),
  #   control = minpack.lm::nls.lm.control(maxiter = 1024)
  # )
    
  FD1 <- function(x) { stats::predict(fit, data.frame(X = x))}
  
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
  
  return(estim_mat)

}
