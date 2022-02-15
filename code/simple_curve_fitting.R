source("code/simple_iso_plot.R")

# Fit a curve to the d18O values
d <- data.frame(X = xdata, Y = ydata)
e <- 2.71

# Oringal: start=list(A=2,M=-5,x_0=1,z=30
fit <- stats::nls(
  Y ~ A * cos(2 * pi * ((X - x_0) / z)) + M,
  data = d,
  start = list(
    A = 5,
    M = 5,
    x_0 = -20,
    z = 25
  ),
  nls.control(maxiter = 100000)
) # If curve is offset from data adjust x_0 variable first then try A

fitdata <- as.data.frame(summary(fit)$coefficients)
fitdata
estimate <- fitdata$Estimate
estimate
A <- estimate[1]
M <- estimate[2]
x_0 <- estimate[3]
z <- estimate[4]

FD1 <- function(x) {
  A * cos(2 * pi * ((x - x_0) / z)) + M
}

# simple exploration plot
simple_iso_plot(xdata, ydata, cdata)

# if data$increment_rej feeds xdata:
curve(
  FD1,
  from = min(xdata - 10, na.rm = T), to = max(xdata + 15, na.rm = T),
  n = 36,
  add = T,
  lty = "dashed",
  lwd = 3
)


# Adapted from https://stackoverflow.com/questions/32613119/plot-the-median-confidence-interval-of-a-bootstrap-output-in-ggplot2
# nlsBoot
curveBoot <- nlstools::nlsBoot(fit, niter = 1000)
# plot(curveBoot)

# Matrix with the bootstrapped parameter estimates
Theta_mat <- curveBoot$coefboot
fun <- function(x, theta) {
  theta["A"] * cos(2 * pi * ((x - theta["x_0"]) / theta["z"])) + theta["M"]
}

# Points where to evaluate the model
x_eval <- seq(min(d$X), max(d$X), length.out = 100)

# Matrix with the predictions
Pred_mat <- apply(Theta_mat, 1, function(theta) fun(x_eval, theta))

# Pack the estimates for plotting
Estims_plot <- cbind(
  x = x_eval,
  as.data.frame(t(apply(Pred_mat, 1, function(y_est) {
    c(
      median_est = median(y_est),
      ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE),
      ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
    )
  })))
)

cols <- c("c1" = "dark blue", "c2" = "#488957", "c3" = "#3399ff")
shapes <- c("s1" = 16, "s2" = 16)

p1 <- ggplot(data = Estims_plot, aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est)) +
  geom_ribbon(alpha = 0.4, fill = "blue") +
  geom_point(data = d, aes(x = X, y = Y), size = rel(2), color = "blue", alpha = 0.8, inherit.aes = FALSE) +
  geom_line(data = d, aes(x = X, y = cdata), size = rel(.5), colour = "black", linetype = "dashed", inherit.aes = FALSE) +
  geom_point(data = d, aes(x = X, y = cdata), size = rel(3.3), colour = "green4", inherit.aes = FALSE) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 18)
  ) +
  ggtitle(paste(isodata$site[1], " - ", isodata$period[1], " - ", isodata$chronology[1]),
          subtitle = paste(isodata$specimen[1], " - ", isodata$taxon[1], " - ", isodata$element[1])
  ) +
  scale_y_continuous(
    expand = c(0, 0), limits = c(-14, 2),
    breaks = seq(-14, 2, by = 2), name = "‰ (VPDB)",
    minor_breaks = NULL
  ) +
  scale_x_continuous(
    expand = c(0, 0), limits = c(-40, 0), labels = c("40", "30", "20", "10", "0"),
    name = c("Distance from root-enamel junction (mm)")
  ) +
  scale_color_manual(
    name = "",
    breaks = c("c1", "c2"),
    values = "blue",
    labels = c("δ18O", "δ13C")
  ) +
  scale_shape_manual(
    name = "",
    breaks = c("s1", "s2"),
    values = shapes,
    labels = c("δ18O", "δ13C")
  )
print(p1)

ggsave(
  file.path("plots", paste("tooth_seq_FINAL_", isodata$specimen[1], ".pdf", sep = "")),
  p1,
  width = 55, height = 40, units = c("cm"), scale = .35, useDingbats = FALSE
)