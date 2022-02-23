#### Part I ####
# The first part of this script focusses on fitting the theoretical seasonality
# model to the empirical data and deriving julian calender days and birth season
# estimates

# Find data
isodata_path <- "data/input/isodata"
isodata_files_paths <- list.files(isodata_path, full.names = T)

# Read data
specimen_overview <- readr::read_csv(
  "data/input/specimen.csv",
  col_types = readr::cols()
)
isodata_list <- purrr::map(
  isodata_files_paths,
  readr::read_csv,
  col_types = readr::cols(specimen = "c")
)

# Correct for Seuss effect if sample is modern
isodata_list_seuss_corrected <- purrr::map(
  isodata_list, 
  function(isodata) {
    chron <- dplyr::filter(
      specimen_overview,
      specimen == isodata$specimen[1]
    )$chronology[1]
    if (!is.na(chron) & chron == "modern") {
      isodata$d13C <- isodata$d13C + 1.5
    }
    return(isodata)
  }
)

# Fit curve for every isodata file
# (parallelized with furrr)
source("code/fit_curve.R")
future::plan(future::multisession, workers = 6) # workers is the number of CPU cores
fitted_curves <- furrr::future_map(
  isodata_list_seuss_corrected,
  function(isodata) {
    # for debugging with a sequential purrr::map
    #message("Trying to fit specimen: ", isodata$specimen[1]) 
    fit_curve(isodata)
  },
  .options = furrr::furrr_options(seed = TRUE)
)

# Plot fitted curves
source("code/plot_single_tooth_with_fitted_curve.R")
purrr::walk2(
  isodata_list_seuss_corrected,
  fitted_curves,
  function(isodata, fitted_curve) {
    isodata_merged <- dplyr::left_join(
      isodata,
      specimen_overview_table,
      by = "specimen"
    )
    p <- plot_single_curve_with_fitted_curve(isodata_merged, fitted_curve$estim_mat)
    ggsave(
      file.path(
        "plots", "isodata_specimen",
        paste("tooth_seq_FINAL_", isodata_merged$specimen[1], ".pdf", sep = "")
      ),
      p, width = 55, height = 40, units = c("cm"), scale = .35, useDingbats = FALSE
    )
  }
)

# Derive some core output parameters:
# The julian calendar equivalent of each sample position
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
      specimen = isodata$specimen[1],
      julian = julian,
      birth = birth,
      birth_season = birth_season
    )
  }
)

# merge intermediate results
specimen_overview_birth <- dplyr::left_join(
  specimen_overview,
  purrr::map_df(
    time_and_birth, function(x) { 
    tibble::tibble(specimen = x$specimen, birth = x$birth, birth_season = x$birth_season)
  }),
  by = "specimen"
)
isodata_julian <- purrr::map2(
  isodata_list_seuss_corrected,
  time_and_birth,
  function(isodata, time_and_birth) {
    dplyr::mutate(
      isodata,
      julian = time_and_birth$julian
    )
  }
)

# write.csv(
#   isodata,
#   file = paste("julian/", isodata$specimen[1], "-julian.csv", sep = ""), 
#   row.names = FALSE
# )

#### Part II ####
# In the second part of this script we derive meaningful summary statistics and
# create some plots from and for the data prepared in part I

library(magrittr)
library(ggplot2)

# Create a super-dataset from the prepared data
all_data_comp <- isodata_julian %>%
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
  dplyr::select(specimen, site, increment) %>%
  dplyr::group_by(site) %>%
  dplyr::summarize(
    n_teeth = dplyr::n_distinct(specimen),
    n_increment = dplyr::n()
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
    )
  )

birth_plot <- ggplot(birth, aes(site, birth)) +
  geom_point(size = 2) +
  ggrepel::geom_label_repel(aes(label = specimen), force = 10, nudge_x = 0.3, size = 1.8) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous(
    limits = c(0, 1), expand = c(0, 0),
    breaks = seq(0, 1, by = .1), name = "Birth Seasonality",
    minor_breaks = NULL
  ) +
  scale_x_discrete(limits = rev(levels(birth$site)), name = "")

ggsave(
  "plots/birth_seasonality_plot_caprines.png",
  birth_plot,
  dpi = 300,
  width = 40, height = 30, units = c("cm"), scale = .35
)
