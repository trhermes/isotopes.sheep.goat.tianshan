# A function to build a simple base-r plot to show isodata
# and visualize adjustments to cosine fitting parameters

simple_iso_plot <- function(xdata, ydata, cdata, FD1){

  par(las = 1)
  par(mar = c(4.5, 4.5, 4.5, 3.5))
  xlimit <- c(-45, 1) # To be used if isodata$increment_rej feeds xdata
  ylimit <- c(-17, 3)
  
  # Set the plot title
  num_samples <- c(length(xdata))
  title <- paste(
    isodata$specimen[1],
    isodata$taxon[1],
    isodata$element[1],
    paste(num_samples, "samples"),
    sep = "  -  "
  )
  subtitle <- paste(
    isodata$site[1],
    isodata$phase[1],
    isodata$period[1],
    isodata$chronology[1],
    sep = "  -  "
  )
  
  # Plot the data
  plot(
    xdata, ydata,
    xlim = xlimit, ylim = ylimit,
    pch = 19,
    col = "dark blue",
    cex = 1.8,
    xlab = "Occlusal surface       -->       REJ",
    ylab = "‰",
    cex.lab = 1.5,
    main = title,
    cex.main = 2.2,
    axes = FALSE,
    cex.lab = 2.2,
    mar = c(5, 2, 4, 2) + 0.2
  )
  par(new = TRUE)
  plot(
    xdata, cdata,
    cex = 1.8,
    axes = FALSE,
    ylab = "",
    xlab = "",
    xlim = xlimit, ylim = ylimit,
    pch = 19,
    col = "dark green"
  )
  mtext(
    subtitle,
    outer = T,
    side = 3,
    line = -5,
    cex = 1.9
  )
  
  # Make the axes
  # axis(side = 1, at = seq(-1,35, by = 1), lwd = 2, cex.axis = 1.6, labels = F)
  axis(
    side = 1,
    at = seq(-45, 1, by = 1),
    lwd = 2,
    cex.axis = 1.6,
    labels = F
  ) # To be used if isodata$increment_rej feeds xdata
  axis(
    side = 2,
    at = seq(-16, 3, by = 1),
    lwd = 2,
    cex.axis = 2.2
  )
  # axis(side = 3, at = seq(-14,-3, by = 1), lwd = 2, cex.axis = 1.4)
  
  # Make the legend
  labels <- c(
    expression(delta^{{ 18 }} * "O"), 
    expression(delta^{{ 13 }} * "C")
  )
  # legend(16,6,legend = labels, pch = c(19,19), col = c("dark blue","dark green"),bty = "n", horiz = T, cex = 2.2, x.intersp = .3, text.width = 1.8, xjust = 0.5)
  
  # To be used if isodata$increment_rej feeds xdata:
  legend(
    -24,
    10,
    legend = labels,
    pch = c(19, 19),
    col = c("dark blue", "dark green"),
    bty = "n",
    horiz = T,
    cex = 2.2,
    x.intersp = .3,
    text.width = 1.8,
    xjust = 0.5
  )
  
  # if data$increment_rej feeds xdata:
  curve(
    FD1,
    from = min(xdata - 10, na.rm = T), to = max(xdata + 15, na.rm = T),
    n = 36,
    add = T,
    lty = "dashed",
    lwd = 3
  )

}
