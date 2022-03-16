#### Part I ####
# The first part of this script focusses on fitting the theoretical seasonality
# model to the empirical data and deriving julian calender days and birth season
# estimates

# Fetch and prepare isotope data to read in via CSVs
source("code/fetch_data.R")

# Produce map
source("code/map.R")

# Read data
# specimen overview table CSV
specimen_overview <- readr::read_csv(
  "data/input/specimen.csv",
  col_types = readr::cols()
)
# measurements per tooth table CSV
isodata_list <- readr::read_csv(
  "data/input/isodata/all_data.csv") %>% 
  group_split(specimen)
  
#   purrr::map(
#   list.files("data/input/isodata", full.names = T),
#   readr::read_csv,
#   col_types = readr::cols(specimen = "c")
# )

# Correct for Seuss effect if sample is modern
isodata_list_seuss_corrected <- purrr::map(
  isodata_list, 
  function(isodata) {
    chron <- dplyr::filter(
      specimen_overview,
      specimen == isodata$specimen[1]
    )$period[1]
    if (!is.na(chron) & chron == "modern") {
      isodata$d13C <- isodata$d13C + 1.5
    }
    return(isodata)
  }
)

# Fit curve for every isodata file
source("code/fit_curve.R")
fitted_curves <- purrr::map(
  isodata_list_seuss_corrected,
  function(isodata) {
    message("Trying to fit specimen: ", isodata$specimen[1]) 
    fit_curve(isodata)
  }
)

# Plot fitted curves
source("code/plot_single_tooth_with_fitted_curve.R")
plot_list <- purrr::map2(
  isodata_list_seuss_corrected,
  fitted_curves,
  function(isodata, fitted_curve) {
    isodata_merged <- dplyr::left_join(
      isodata,
      specimen_overview,
      by = "specimen"
    )
    p <- plot_single_curve_with_fitted_curve(isodata_merged, fitted_curve$estim_mat)
    ggsave(
      file.path(
        "plots", "isodata_specimen",
        paste("isotope_tooth_seq_", isodata_merged$specimen[1], ".pdf", sep = "")
      ),
      p, width = 55, height = 40, units = c("cm"), scale = .35, 
      device = grDevices::cairo_pdf
    )
  return(p)}
)

# Create plot matrix for new data plots (Chap & Jeti-Orguz)

# Select and order the desired plots
new_data_plot_nums <- specimen_overview %>%
  dplyr::group_by(specimen) %>% 
  dplyr::mutate(position = cur_group_id()) %>% 
  dplyr::select(specimen, position, site) %>% 
  dplyr::filter(site %in% c("Chap", "Jeti-Oguz")) %>% 
  dplyr::pull(position)

plot_list_selection <- plot_list[min(new_data_plot_nums):max(new_data_plot_nums)]

# Create a hacky "legend" plot
legend_plot <- tibble::tibble(
  x = c(1,1,1),
  y = 3:1,
  label = c("Isotope values", "δ13C", "δ18O")
) %>% 
  ggplot() +
  geom_point(aes(x - 1, y, color = label), size = 8) +
  geom_text(aes(x + 1, y, label = label), size = 8) +
  scale_color_manual(
    values = c("δ18O" = "blue", "δ13C" = "green4"),
    na.translate = FALSE
  ) +
  coord_cartesian(xlim = c(-4,8), ylim = c(-3,5)) +
  theme_nothing()

# Merge plot components and write to a file
grid_plot <- cowplot::plot_grid(
  plotlist = append(plot_list_selection, list(legend_plot)),
  ncol = 3
)

ggsave(
  "plots/Figure2.pdf", grid_plot, 
  scale = 4, 
  width = 16, 
  height = 19, 
  units = "cm", 
  bg = "white",
  device = grDevices::cairo_pdf
)

# Derive the julian calendar equivalent of each sampling position on each tooth
julian_for_each_specimen <- purrr::map2(
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
    return(julian)
  }
)

# Merge julian into isodata tables
isodata_julian <- purrr::map2(
  isodata_list_seuss_corrected,
  julian_for_each_specimen,
  function(isodata, julian) {
    dplyr::mutate(
      isodata,
      julian = julian
    )
  }
)

# Helper function to calculate birth value (proxy for birth season)
derive_birth <- function(period, phase_shift) {
  one_max_pos <- phase_shift
  multi_max_pos <- one_max_pos + seq(-5,3,1) * period
  last_max <- multi_max_pos[
    tail(which(multi_max_pos < 0), n = 1)
  ]
  birth <- abs(last_max) / period
  return(birth)
}

# Helper function to get the birth season from a birth value
determine_birth_season <- function(birth) {
  dplyr::case_when(
    birth <  0.13 | birth >= 0.87 ~ "winter",
    birth >= 0.13 & birth <  0.38 ~ "spring",
    birth >= 0.38 & birth <  0.63 ~ "summer",
    birth >= 0.63 & birth <  0.87 ~ "fall",
  )
}

# Derive the birth season proxy
birth_proxy_for_each_specimen <- purrr::map_df(
  fitted_curves,
  function(fitted_curve) {
    # Calculate mean birth estimate based on simple curve fit
    period <- fitted_curve$fit$m$getPars()[["z"]] # identical to the length of one year in mm
    phase_shift <- fitted_curve$fit$m$getPars()[["x_0"]]
    birth_simple_fit <- derive_birth(period, phase_shift)
    # Calculate birth estimate based on fits obtained via bootstrap resampling
    periods <- fitted_curve$theta_mat$z
    phase_shifts <- fitted_curve$theta_mat$x_0
    births_resampling <- purrr::map2_dbl(periods, phase_shifts, derive_birth)
    birth_resampling_mean <- mean(births_resampling)
    birth_resampling_sd <- sd(births_resampling)
    # compile output
    tibble::tibble(
      specimen = fitted_curve$specimen,
      birth_simple_fit = birth_simple_fit,
      birth_season_simple_fit = determine_birth_season(birth_simple_fit),
      birth_resampling_mean = birth_resampling_mean,
      birth_resampling_sd = birth_resampling_sd,
      birth_season_resampling_mean = determine_birth_season(birth_resampling_mean)
    )
  }
)

# Merge birth info into specimen overview table 
specimen_overview_birth <- dplyr::left_join(
  specimen_overview,
  birth_proxy_for_each_specimen,
  by = "specimen"
)


#### Part II ####
# In the second part of this script we derive meaningful summary statistics and
# create some plots from and for the data prepared in part I

library(magrittr)
library(ggplot2)

# Summarize (mean) duplicate entries per increment
isodata_julian_summarized <- purrr::map(
  isodata_julian,
  function(isodata) {
    isodata %>%
    dplyr::group_by(increment) %>%
    dplyr::summarise(
      increment = dplyr::first(increment),
      specimen = dplyr::first(specimen),
      d13C = mean(d13C),
      d18O = mean(d18O),
      measure = mean(measure),
      julian = mean(julian)
    )
  }
)

# Create a super-dataset from the prepared data
all_data_comp <- isodata_julian_summarized %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(
    specimen_overview_birth,
    by = "specimen"
  )

# Write Table S3
write.csv(
  all_data_comp %>%
    dplyr::mutate(
      dplyr::across(
        tidyselect:::where(is.numeric),
        round, digits = 2
      )
    ),
  "tables/Table_S3.csv",
  row.names = F, quote = F
)

# Make df with only Chap data
all_data <- all_data_comp %>% dplyr::filter(site %in% c("Chap", "Jeti-Oguz"))

# Generate summary stats figures and tables

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
  dplyr::select(d13C, d18O, site) %>%
  dplyr::group_by(site) %>%
  dplyr::summarize(
    dplyr::across(
      tidyselect::starts_with("d"),
      list(mean = mean, sd = stats::sd, min = min, max = max)
    )
  )

# Count unique teeth per site
counts <- all_data_comp %>%
  dplyr::select(specimen, site, increment, individual) %>%
  dplyr::group_by(site) %>%
  dplyr::summarize(
    n_teeth = dplyr::n_distinct(specimen),
    n_increment = dplyr::n(),
    n_mni = dplyr::n_distinct(individual)
  )

# Combine sum_stats and counts
sum_stats_final <- sum_stats %>%
  dplyr::left_join(counts, by = "site") %>%
  dplyr::arrange(match(site, c("Chap", "Jeti-Oguz", "Bayan-Zherek", "Begash", "Dali", "Kent", "Turgen"))) %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::where(is.numeric),
      round, digits = 2
    )
  )

# Write Chap summary stats
write.csv(
  sum_stats_final,
  "tables/Chap_tooth_summary_stats.csv",
  row.names = F, quote = F
)

# Average isotopic change per tooth
iso_change <- all_data_comp %>%
  dplyr::select(d13C, d18O, site, specimen) %>%
  dplyr::group_by(specimen, site) %>%
  dplyr::summarize(
    dplyr::across(
      tidyselect::starts_with("d"), 
      list(mean = mean, sd = sd, min = min, max = max)
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(d13c_range = d13C_max - d13C_min)

# Min and max isotope range per site
site_iso_change <- iso_change %>%
  dplyr::select(d13c_range, site) %>%
  dplyr::group_by(site) %>%
  dplyr::summarize(
    dplyr::across(
      tidyselect::starts_with("d"),
      list(min = min, max = max)
    )
  )

# Write Table S4
write.csv(
  site_iso_change,
  "tables/Table_S4.csv",
  row.names = F, quote = F
)

# Pearson's r correlation test between d13C and d18O values for each tooth
# Informed from https://dominicroye.github.io/en/2019/tidy-correlation-tests-in-r/
cor_fun <- function(df) { cor.test(df$d13C, df$d18O, method = "pearson") %>% broom::tidy() }
corr_tests <- all_data_comp %>%
  dplyr::select(d13C, d18O, specimen, site) %>%
  dplyr::group_nest(specimen, site) %>%
  dplyr::mutate(model = purrr::map(data, cor_fun)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest(cols = c(model)) %>%
  dplyr::mutate(
    significant = ifelse(p.value <= .05, "pass", "fail")
  ) %>%
  dplyr::transmute(
    specimen,
    site,
    Pearson_r = estimate,
    conf.low, conf.high,
    test_stat = statistic,
    df = parameter,
    p.value, significant
  ) %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::where(is.numeric),
      round, digits = 4
    )
  )

# Write Table S5
write.csv(
  corr_tests,
  "tables/Table_S5.csv",
  row.names = F,
  quote = F
)

# Count pass and fail
# corr_tests %>% dplyr::count(significant)

# Plot the isotopic data per tooth x~y
corr_tests_plot <- all_data_comp %>%
  dplyr::select(d13C, d18O, specimen, site) %>%
  ggplot() +
    geom_point(mapping = aes(d13C, d18O)) +
    xlim(-15, 0) +
    theme_bw() +
    coord_fixed() +
    facet_wrap(~ specimen + site)

ggsave(
  "plots/Fig_S1.pdf",
  corr_tests_plot,
  dpi = 300, width = 20, height = 30, units = c("cm"), scale = .85
)

# Birth seasonality chart for all sites
birth <- specimen_overview_birth %>%
  dplyr::filter(element == "M/2") %>%
  dplyr::mutate(
    site = factor(
      site,
      levels = c("Chap", "Jeti-Oguz", "Bayan-Zherek", "Begash", "Dali", "Kent", "Turgen")
    ),
    specimen = forcats::fct_reorder(specimen, birth_resampling_mean)
  )

birth_plot <- birth %>%
  ggplot() +
  facet_grid(
    rows = dplyr::vars(site),
    space = "free",
    scales = "free_y"
  ) +
  geom_errorbarh(
    mapping = aes(
      y = specimen,
      xmin = birth_resampling_mean - birth_resampling_sd,    # 2-sigma, for 1-sigma insert /2
      xmax = birth_resampling_mean + birth_resampling_sd     # 2-sigma, for 1-sigma insert /2
    ),
    size = 0.3, height = 0.3
  ) +
  geom_point(
    aes(birth_resampling_mean, specimen),
    size = 1.5, shape = 16
  ) +
  # geom_point(
  #   aes(birth_simple_fit, specimen),
  #   size = 2,
  #   colour = "red"
  # ) +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0)
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(0, 1, by = .1), name = "Birth Seasonality",
    minor_breaks = NULL
  ) +
  coord_cartesian(
    xlim = c(0, 1)
  )

print(birth_plot)

ggsave(
  "plots/Figure3.png",
  birth_plot,
  dpi = 300,
  width = 40, height = 30, units = c("cm"), scale = .4
)
