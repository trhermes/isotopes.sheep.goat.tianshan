# Find data
isodata_path <- "data/input/isodata"
isodata_files_paths <- list.files(isodata_path, full.names = T)

# Read data
isodata_list <- purrr::map(
  isodata_files_paths,
  readr::read_csv,
  col_types = readr::cols()
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

# derive some core output parameters:
# the julian calendar equivalent of each sample position
# and the birth season proxy
time_and_birth <- purrr::map2(
  isodata_list_seuss_corrected,
  fitted_curves,
  function(isodata, fitted_curve) {
    # basic parameters
    period <- fitted_curve$fit$m$getPars()[["z"]] # identical to the length of one year in mm
    phase_shift <- fitted_curve$fit$m$getPars()[["x_0"]]
    oldest_meas_point <- min(-isodata$measure)
    # adjust curve to year by equating the minimum value of the curve with the 15th of January
    one_min_pos <- -period / 2 + phase_shift
    multi_min_pos <- one_min_pos + seq(-5,3,1) * period
    fift_jan_before_birth <- multi_min_pos[
      tail(which(multi_min_pos < oldest_meas_point), n = 1)
    ]
    # derive sampling days in julian calender format
    distance_to_fift_jan <- abs(fift_jan_before_birth - (-isodata$measure))
    day_sampled_in_julian_calender <- -365 + 15 + distance_to_fift_jan/period * 365
    julian <- day_sampled_in_julian_calender %% 365 # recycle year
    # calculate birth value (proxy for birth season)
    one_max_pos <- phase_shift
    multi_max_pos <- one_max_pos + seq(-5,3,1) * period
    last_max <- multi_max_pos[
      tail(which(multi_max_pos < 0), n = 1)
    ]
    birth <- abs(last_max) / period
    # translate birth value to a season
    birth_season <- dplyr::case_when(
      birth <  0.13 | birth >= 0.87 ~ "winter",
      birth >= 0.13 & birth <  0.38 ~ "spring",
      birth >= 0.38 & birth <  0.63 ~ "summer",
      birth >= 0.63 & birth <  0.87 ~ "fall",
    )
    # compile output
    list(
      julian = julian,
      birth = birth,
      birth_season = birth_season
    )
  }
)

###

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
