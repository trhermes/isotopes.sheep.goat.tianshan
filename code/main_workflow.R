library(bayesboot)
library(nlstools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)
library(broom)
library(corrr)

# Load data
file <- "data/T032_C_O_meas.csv"
isodata <- readr::read_csv(file)
ydata <- isodata$d18O
xdata <- -isodata$measure # isodata$Increment_ is a vector of the the increments offset against 32 max increments of tooth TRH-14

# Correct for Seuss effect if sample is modern
if (isodata$chronology[1] == "modern") {
  isodata$d13C <- isodata$d13C + 1.5
}
cdata <- data$d13C


# Build simple base-r plot to show isodata and visualize adjustments to cosine fitting parameters
par(las = 1)
par(mar = c(4.5, 4.5, 4.5, 3.5))
xlimit <- c(-45, 1) # To be used if isodata$increment_rej feeds xdata
ylimit <- c(-17, 3)

# Set the plot title
num_samples <- c(length(xdata))
title <- paste(
  isodata$specimen[1], " - ", isodata$taxon[1], " - ", isodata$element[1], " - ", num_samples, "samples"
)
subtitle <- paste(isodata$site[1], " - ", isodata$phase[1], " - ", isodata$period[1], " - ", isodata$chronology[1])

# Plot the data
plot(xdata, ydata, xlim = xlimit, ylim = ylimit, pch = 19, col = "dark blue", cex = 1.8, xlab = "Occlusal surface       -->       REJ", ylab = "‰", cex.lab = 1.5, main = title, cex.main = 2.2, axes = FALSE, cex.lab = 2.2, mar = c(5, 2, 4, 2) + 0.2)
par(new = TRUE)
plot(xdata, cdata, cex = 1.8, axes = FALSE, ylab = "", xlab = "", xlim = xlimit, ylim = ylimit, pch = 19, col = "dark green")
mtext(subtitle, outer = T, side = 3, line = -5, cex = 1.9)

# Make the axes
# axis(side = 1, at = seq(-1,35, by = 1), lwd = 2, cex.axis = 1.6, labels = F)
axis(side = 1, at = seq(-45, 1, by = 1), lwd = 2, cex.axis = 1.6, labels = F) # To be used if isodata$increment_rej feeds xdata
axis(side = 2, at = seq(-16, 3, by = 1), lwd = 2, cex.axis = 2.2)
# axis(side = 3, at = seq(-14,-3, by = 1), lwd = 2, cex.axis = 1.4)

# Make the legend
labels <- c(expression(delta^{{ 18 }} * "O"), expression(delta^{{ 13 }} * "C"))
# legend(16,6,legend = labels, pch = c(19,19), col = c("dark blue","dark green"),bty = "n", horiz = T, cex = 2.2, x.intersp = .3, text.width = 1.8, xjust = 0.5)

# To be used if isodata$increment_rej feeds xdata:
legend(-24, 10, legend = labels, pch = c(19, 19), col = c("dark blue", "dark green"), bty = "n", horiz = T, cex = 2.2, x.intersp = .3, text.width = 1.8, xjust = 0.5)

d <<- data.frame(X = xdata, Y = ydata)

################################
# Fit a curve to the d18O values
e <- 2.71
# Oringal: start=list(A=2,M=-5,x_0=1,z=30
fit <- nls(Y ~ A * cos(2 * pi * ((X - x_0) / z)) + M, data = d, start = list(A = 5, M = 5, x_0 = -20, z = 25), nls.control(maxiter = 100000)) # If curve is offset from data adjust x_0 variable first then try A
fitdata <- as.data.frame(summary(fit)$coefficients)
fitdata
estimate <- fitdata$Estimate
estimate
A <- estimate[1]
M <- estimate[2]
x_0 <- estimate[3]
z <- estimate[4]

FD1 <- function(x, y) {
  A * cos(2 * pi * ((x - x_0) / z)) + M
}
# FD1=function(x,y){(A*((X-x_b)/x_a))*cos(2*pi*((X-x_0)/z))+M+(p*X)}
# curve(FD1, from=-25,to=max(xdata,na.rm=T), n=36, add=T, lty = "dashed", lwd = 3)

# if data$increment_rej feeds xdata:
curve(FD1, from = min(xdata - 10, na.rm = T), to = max(xdata + 15, na.rm = T), n = 36, add = T, lty = "dashed", lwd = 3)

###############################################
###############################################
###############################################

# Adapted from https://stackoverflow.com/questions/32613119/plot-the-median-confidence-interval-of-a-bootstrap-output-in-ggplot2
# nlsBoot
curveBoot <- nlsBoot(fit, niter = 1000)
# plot(curveBoot)

# Matrix with the bootstrapped parameter estimates
Theta_mat <- curveBoot$coefboot
fun <- function(x, theta) theta["A"] * cos(2 * pi * ((x - theta["x_0"]) / theta["z"])) + theta["M"]

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

ggsave(paste("tooth_seq_FINAL_", isodata$specimen[1], ".pdf", sep = ""), p1,
  width = 55, height = 40, units = c("cm"), scale = .35, useDingbats = FALSE
)



######
######
######
# Calculate how much tooth distance is represented in the period of the fitted curve using data pointa
# If data point do not include min and max values, they must be obtained by extrapolation
######
######
######
# Find the min and max (period) of the fitted curve within the range of xdata first then outside it

if (optimize(FD1, interval = c(min(xdata), max(xdata)))$minimum != optimize(FD1, interval = c(-40, 0))$minimum) {
  curvemin <- optimize(FD1, interval = c(-60, 0))
} else {
  curvemin <- optimize(FD1, interval = c(min(xdata), max(xdata)))
}

if (optimize(FD1, interval = c(min(xdata), max(xdata)), maximum = TRUE)$maximum != optimize(FD1, interval = c(-40, 0), maximum = TRUE)$maximum) {
  curvemax <- optimize(FD1, interval = c(-60, 0), maximum = TRUE)
} else {
  curvemax <- optimize(FD1, interval = c(min(xdata), max(xdata)), maximum = TRUE)
}

curvemax$maximum
curvemin$minimum

# Optimize returns two values. The min or max, and the objective. The objective
# is the value of curve, whereas the min or max are the x values of where they are
# curveminmax=c(curvemin$minimum,curvemax$maximum)
# This depends on whether the max is after the min, either add or substract accordingly
curveperiod <- abs(curvemax$maximum - curvemin$minimum)
curveperiod
isodata$curveperiod <- curveperiod

# Convert xdata to Julian days and set min temp day to Jan 15th
julian <- ((xdata - curvemin$minimum) / (curveperiod / 180)) + 15

birth <- abs(curvemax$maximum / (2 * curveperiod))

# Because the min and max may fall outside of the range of xdata, birth may be greater than 1
# if so, then simply subtract 1 to get a value between 0 and 1
if (birth > 1) {
  birth <- birth - 1
}

for (i in 1:length(julian)) {
  if (julian[i] <= 1) {
    julian[i] <- julian[i] + 365
  }
  if (julian[i] >= 366) {
    julian[i] <- julian[i] - 365
  }
}
isodata$julian <- julian
isodata$birth <- birth

write.csv(isodata, file = paste("julian/", isodata$specimen[1], "-julian.csv", sep = ""), row.names = FALSE)


#######################################
#######################################
#######################################

# Combine julian CSVs into one file for Chap and comparative data
all_data_comp <- list.files(path = "~/Dropbox (MPI SHH)/Margins or Nodes/Chap tooth SIA/01-data/julian/", full.names = TRUE) %>%
  lapply(read_csv, col_types = cols(specimen = "c", phase = "c")) %>%
  bind_rows()

# Write table (Table 1)
write.csv(all_data_comp, "../Figures/Table S3 - all sites isotope data.csv", row.names = F, quote = F)

# Make df with only Chap data
all_data <- all_data_comp %>% filter(site == "Chap" | site == "Jeti-Oguz")
# Generate summary stats figures and tables
#
#
#

# d18O boxplot for new data
d18O_box <- ggplot(all_data, aes(x = d18O, y = taxon)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(shape = 16, height = 0.1) +
  theme_bw()

# d13C boxplot for new data
d13C_box <- ggplot(all_data, aes(x = d13C, y = taxon)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(shape = 16, height = 0.1) +
  theme_bw()

# Summary stats table for new data
sum_stats <- all_data_comp %>%
  select(d13C, d18O, site) %>%
  group_by(site) %>%
  summarize(across(starts_with("d"), list(mean = mean, sd = sd, min = min, max = max)))

# Count unique teeth per site
counts <- all_data_comp %>%
  select(specimen, site, increment) %>%
  group_by(site) %>%
  summarize(n_teeth = n_distinct(specimen), n_increment = n())

# Combine sum_stats and counts
sum_stats_final <- cbind(sum_stats, counts) %>%
  select(-10) %>% # Remove second instance of "site" by its index
  arrange(match(site, c("Chap", "Jeti-Oguz", "Bayan-Zherek", "Begash", "Dali", "Kent", "Turgen"))) %>%
  mutate_if(is.numeric, round, digits = 2)

# Write table (Table 1)
write.csv(sum_stats_final, "../Figures/Chap_tooth_summary_stats.txt", row.names = F, quote = F)

# Average isotopic change per tooth
iso_change <- all_data_comp %>%
  select(d13C, d18O, site, specimen) %>%
  group_by(specimen, site) %>%
  summarize(across(starts_with("d"), list(mean = mean, sd = sd, min = min, max = max))) %>%
  mutate(d13c_range = d13C_max - d13C_min)

# Min and max isotope range per site
site_iso_change <- iso_change %>%
  select(d13c_range, site) %>%
  group_by(site) %>%
  summarize(across(starts_with("d"), list(min = min, max = max)))
write.csv(site_iso_change, "../Figures/Table S4 - individual min and max d13C ranges per site.csv", row.names = F, quote = F)

# Pearson's r correlation test between d13C and d18O values for each tooth
# Informed from https://dominicroye.github.io/en/2019/tidy-correlation-tests-in-r/
cor_fun <- function(df) cor.test(df$d13C, df$d18O, method = "pearson") %>% tidy()
corr_tests <- all_data_comp %>%
  select(d13C, d18O, specimen, site) %>%
  group_nest(specimen, site) %>%
  mutate(model = map(data, cor_fun)) %>%
  select(-data) %>%
  unnest(cols = c(model)) %>%
  mutate(significant = ifelse(p.value <= .05, "pass", "fail")) %>%
  select(specimen, site, Pearson_r = estimate, conf.low, conf.high, test_stat = statistic, df = parameter, p.value, significant, -method, -alternative) %>%
  mutate(across(where(is.numeric), round, 4))
write.csv(corr_tests, "../Figures/Table S5 - pearsons r between d13C and d18O values per tooth.csv", row.names = F, quote = F)

# Count pass and fail
corr_tests %>% count(significant)

# Plot the isotopic data per tooth x~y
corr_tests_plot <- all_data_comp %>%
  select(d13C, d18O, specimen, site) %>%
  group_nest(specimen, site) %>%
  mutate(., model = map(data, cor_fun)) %>%
  select(-model) %>%
  unnest() %>%
  ggplot(aes(d13C, d18O)) +
  geom_point() +
  xlim(-15, 0) +
  theme_bw() +
  coord_fixed() +
  facet_wrap(~ specimen + site)
ggsave("../Figures/Fig S1 - scatter plots of d13C and d18O values by tooth.pdf", corr_tests_plot,
  dpi = 300,
  width = 20, height = 30, units = c("cm"), scale = .85
)

# Birth seasonality chart for all sites -- all_data_comp
birth <- all_data_comp %>%
  filter(element == "M/2", increment == 1) %>%
  distinct(increment, specimen, .keep_all = T)
birth$site <- factor(birth$site, levels = c("Chap", "Jeti-Oguz", "Bayan-Zherek", "Begash", "Dali", "Kent", "Turgen"))

birth_plot <- ggplot(birth, aes(site, birth)) +
  geom_point(size = 2) +
  geom_label_repel(aes(label = specimen), force = 10, nudge_x = 0.3, size = 1.8) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous(
    limits = c(0, 1), expand = c(0, 0),
    breaks = seq(0, 1, by = .1), name = "Birth Seasonality",
    minor_breaks = NULL
  ) +
  scale_x_discrete(limits = rev(levels(birth$site)), name = "")
print(birth_plot)

ggsave("../Figures/birth_seasonality_plot_caprines.png", birth_plot,
  dpi = 300,
  width = 40, height = 30, units = c("cm"), scale = .35
)
