library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(docxtractr)
library(tibble)
library(readr)
library(readxl)
library(janitor)

# import data from Leibniz Labor Kiel Univeristy; data sheet prepared by Dr. Nils Andersen
leibniz_data <- readxl::read_excel("~/Dropbox (MPI SHH)/Margins or Nodes/Chap tooth SIA/01-data/zTHEC2007_IV_data.xlsx",
                                   range = "A10:X231") %>% 
  clean_names() %>% 
  rename(
    run_counter = 1,
    period = 8,
    run_ID = 15,
    subsample_analyzed = 16,
    d13C = 18,
    errorc = 19,
    d18O = 20,
    erroro = 21) %>% 
  select(-10,-11) %>% 
  remove_empty("rows")

# split up isotope data from animal teeth and isotope data for lab standards 
# line to cut df
n = 163

study_data <- leibniz_data[row.names(leibniz_data) %in% 1:n, ] 
 
standards_stats <- leibniz_data[row.names(leibniz_data) %in% (n+2):nrow(leibniz_data), ] %>% 
  remove_empty("cols") %>% 
  filter(x44_m_v_pr_1_cycle>1000) %>% 
  group_by(run_ID) %>% 
  summarize(mean_d13C = mean(d13C),
            sd_d13C = sd(d13C),
            mean_d18O = mean(d18O),
            sd_d18O = sd(d18O),
            n = n()
            )

dir.create("tables")
write_csv(standards_stats, "tables/isotope_lab_run_standards_THEC2007_Kiel_Leibniz_labor.csv")


