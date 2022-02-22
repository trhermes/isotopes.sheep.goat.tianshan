# Find data
isodata_path <- "data/input/isodata"
isodata_files_paths <- list.files(isodata_path, full.names = T)

# Read data
isodata_list <- purrr::map(
  isodata_files_paths,
  readr::read_csv
)

# Correct for Seuss effect if sample is modern
isodata_list_seuss_corrected <- purrr::map(
  isodata_list, 
  function(isodata) {
    if (!is.na(isodata$chronology[1]) & isodata$chronology[1] == "modern") {
      isodata$d13C <- isodata$d13C + 1.5
    }
    return(isodata)
  }
)

# Fit curve for every isodata file
source("code/fit_curve.R")
fitted_curves <- purrr::map2(
  isodata_files_paths,
  isodata_list_seuss_corrected,
  function(filepath, isodata) {
    message("Trying to fit: ", filepath)
    fit_curve(isodata)
  }
)

# plot for debugging
# source("code/simple_iso_plot.R")
# simple_iso_plot(d$X, d$Y, isodata$d13C, FD1)

# plot fitted curves
source("code/plot_single_tooth_with_fitted_curve.R")
purrr::walk2(
  isodata_list_seuss_corrected,
  fitted_curves,
  function(isodata, fitted_curve) {
    p <- plot_single_curve_with_fitted_curve(isodata, fitted_curve$estim_mat)
    ggsave(
      file.path("plots", paste("tooth_seq_FINAL_", isodata$specimen[1], ".pdf", sep = "")),
      p, width = 55, height = 40, units = c("cm"), scale = .35, useDingbats = FALSE
    )
  }
)

# Calculate how much tooth distance is represented in the period of the fitted curve
# If data point do not include min and max values, they must be obtained by extrapolation

curve_periods <- purrr::map2(
  isodata_list_seuss_corrected,
  fitted_curves,
  function(isodata, fitted_curve) {
    
    FD1 = function(x) { stats::predict(
      fitted_curve$fit,
      data.frame(X = x)
    ) }
    
    # Find the min and max (period) of the fitted curve within the range of xdata first then outside it
    covered_range <- range(-isodata$measure)
    total_range   <- c(-40, 0)
    super_range   <- c(-60, 0) # TODO: why are total_range and super_range different?
    tolerance     <- 0.1
    
    # Optimize returns two values. The min or max, and the objective. The objective
    # is the value of curve, whereas the min or max are the x values of where they are
    # curveminmax=c(curvemin$minimum,curvemax$maximum)
    min_in_covered <- optimize(FD1, interval = covered_range)$minimum
    min_in_total   <- optimize(FD1, interval = total_range)$minimum
    max_in_covered <- optimize(FD1, interval = covered_range, maximum = TRUE)$maximum
    max_in_total   <- optimize(FD1, interval = total_range, maximum = TRUE)$maximum
    
    curvemin <- if (abs(min_in_covered - min_in_total) > tolerance) {
      optimize(FD1, interval = super_range)$minimum
    } else {
      min_in_covered
    }
    
    curvemax <- if (abs(max_in_covered - max_in_total) > tolerance) {
      optimize(FD1, interval = super_range, maximum = TRUE)$maximum
    } else {
      max_in_covered
    }
    
    # This depends on whether the max is after the min, either add or substract accordingly
    curveperiod <- abs(curvemax - curvemin)
    
    return(curveperiod)
  }
)

curve_periods <- purrr::map2(
  isodata_list_seuss_corrected,
  fitted_curves,
  function(isodata, fitted_curve) {
    # curve parameters
    period <- fitted_curve$fit$m$getPars()[["z"]]
    phase_shift <- fitted_curve$fit$m$getPars()[["x_0"]]
    # positions
    oldest_meas_point <- min(-isodata$measure)
    youngest_meas_point <- max(-isodata$measure)
    # ranges
    tooth_length <- abs(oldest_meas_point - youngest_meas_point)
    length_of_year_in_mm <- period
    tooth_life_in_years <- tooth_length/length_of_year_in_mm
    # derived values
    fifteenth_jan_pos_mm <- -period/2 + phase_shift
    mm_per_day <- tooth_life_in_years/365
    
    # plot(-50:10, FD1(-50:10, fitted_curve$fit$m$getPars()[["A"]], fitted_curve$fit$m$getPars()[["x_0"]], fitted_curve$fit$m$getPars()[["z"]], fitted_curve$fit$m$getPars()[["M"]]))
  }
)

mm_pos_to_julian <- function()

#curve_periods[[1]] -> curveperiod

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
