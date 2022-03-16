library(ggplot2)

plot_single_curve_with_fitted_curve <- function(isodata, estim_mat) {

  cols <- c("c1" = "dark blue", "c2" = "#488957", "c3" = "#3399ff")
  shapes <- c("s1" = 20, "s2" = 20)
  
  ggplot() +
    geom_ribbon(
      data = estim_mat,
      aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est),
      alpha = 0.4, fill = "blue"
    ) +
    geom_point(
      data = isodata,
      aes(x = -measure, y = d18O), size = rel(3.6),
      color = "blue", alpha = 0.8, inherit.aes = FALSE
    ) +
    geom_line(
      data = isodata,
      aes(x = -measure, y = d13C), size = rel(.5),
      colour = "black", linetype = "dashed", inherit.aes = FALSE
    ) +
    geom_point(
      data = isodata,
      aes(x = -measure, y = d13C), size = rel(3.6),
      colour = "green4", inherit.aes = FALSE
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 18)
    ) +
    ggtitle(
      paste(isodata$site[1], " - ", isodata$period[1]),
      subtitle = paste(isodata$specimen[1], " - ", isodata$taxon[1], " - ", isodata$element[1])
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = seq(-14, 2, by = 2), name = "‰ (VPDB)",
      minor_breaks = NULL
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      labels = c("40", "30", "20", "10", "0"),
      name = c("Distance from root-enamel junction (mm)")
    ) +
    coord_cartesian(
      xlim = c(-40, 0),
      ylim = c(-15, 2)
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

}
