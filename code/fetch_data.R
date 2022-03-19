# This script will import isotope data from an Excel spreadsheet prepared by
# the Leibniz labor at Kiel University
# This script also downloads comparative isotope data from 
# Ventresca-Miller et al. 2020, doi:10.1080/20548923.2020.1759316
# Hermes et al. 2019, doi:10.1098/rspb.2019.1273
# and generate two CSV files to hold 1) data and 2) metadata (linked by specimen ID),
# which is then combined with the new data

library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(docxtractr)
library(tibble)
library(readr)
library(readxl)
library(janitor)

#### import new data from Leibniz Labor Kiel Univeristy ####
# data sheet prepared by Dr. Nils Andersen

leibniz_data <- readxl::read_excel("data/raw/zTHEC2007_IV_data.xlsx",
                                   range = "A10:X231") %>% 
  janitor::clean_names() %>% 
  dplyr::rename(
    run_counter = 1,
    specimen = tooth_id,
    measure = meas,
    period = chronology,
    run_ID = 15,
    subsample_analyzed = 16,
    d13C = 18,
    errorc = 19,
    d18O = 20,
    erroro = 21) %>% 
  dplyr::select(-10,-11) %>% 
  janitor::remove_empty("rows") %>% 
  dplyr::mutate(specimen = gsub('-', '', specimen)) %>% 
  dplyr::mutate(site = gsub('-I', '', site))

# split up isotope data from animal teeth and isotope data for lab standards 
#
n = 163   # line to cut df
study_data_ <- leibniz_data[row.names(leibniz_data) %in% 1:n, ] %>% 
  dplyr::mutate(period = stringr::str_to_title(period))

#### Produce summary stats for isotope lab runs of standards ####
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
readr::write_csv(standards_stats %>% 
                   dplyr::mutate(
                     dplyr::across(
                       tidyselect:::where(is.numeric),
                       round, digits = 4
                     )
                   ), "tables/Table_S2_isotope_lab_run_standards_THEC2007_Kiel_Leibniz_labor.csv")


#### Import new metadata for samples ####

new_metadata_ <- readxl::read_excel("data/raw/Chap_analyzed_teeth_2020_07.xlsx") %>% 
  clean_names() %>% 
  rename(specimen = 1) %>% 
  mutate(specimen = gsub('Chap-', "CHP", specimen)) %>% 
  mutate(specimen = gsub('Ibex-', "IBX", specimen))
readr::write_csv(new_metadata_, "tables/Table_S1_metadata_newly_analyzed_caprine_teeth.csv")

#### Download comparative data from Ventresca-Miller et al. 2020 ####
# This data set is archived as two tables embedded in Microsoft Word documents

# name of dir and file
#
dir.create("data/comparative/vm", recursive = T)
vm_supp_files <- "data/comparative/vm/vm_supp_data.zip"


# download file (Web server does not like wget or curl request, but can be overcome with unique User-Agent)
#
download.file("https://www.tandfonline.com/doi/suppl/10.1080/20548923.2020.1759316/suppl_file/ysta_a_1759316_sm2281.zip",
              vm_supp_files, 
              method = "curl", 
              extra = paste0("-L -H ", '"User-Agent: Mozilla ', date(), '"'))

# directory to hold files in zip archive
#
vm_dir <- "data/comparative/vm"
dir.create(vm_dir, recursive = T)

# unzip
#
unzip(vm_supp_files, 
      exdir = vm_dir)

# function to read in data from tables in MS Word docs, rm blank rows, and insert site name fetched from file names
#
vm_tbl_manip <- function(x) { 
  read_docx(x) %>% 
    docx_extract_tbl %>%
    dplyr::na_if("") %>% 
    na.omit() %>% 
    dplyr::mutate(site = str_extract(x,'(?<=_)[^_]+(?=\\.)'))   # Regex to cut out site name from file names containing data
}

# acquire data from MS Word tables, list collapse, standardize column names and order
#
vm_comp_data_ <- list.files(vm_dir, full.names = T, pattern = "\\.docx$") %>%
  purrr::map(vm_tbl_manip) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate_at(vars(matches("Sample.Number.by.Increment")), 
                   funs(as.numeric(str_trunc(.,width = 2, side = "left", ellipsis = '')))) %>% 
  dplyr::relocate(increment = `Sample.Number.by.Increment`,
                specimen = `Sample.number`,
                d13C = `δ13C..permil.VPDB.`,
                d18O = `δ18O..permil.VPDB.`,
                measure = `Distance.from.ERJ..mm.`) %>%
  group_by(specimen) %>% 
#  filter(n() > 4) %>%             # specimens 5582 & 5583 have 4 measurements, not enough for curve fitting
  readr::type_convert()

#### Download comparative data from Hermes et al. 2019 ####
#

# name of file
#
h_supp_data <- "h_supp_data.csv"

# directory to hold data
#
h_dir <- "data/comparative/h"
dir.create(h_dir, recursive = T)

# download file (Webserver does not like wget or curl request, but can be overcome with unique User-Agent)
# download goes into h_dir
#
download.file("https://royalsocietypublishing.org/action/downloadSupplement?doi=10.1098%2Frspb.2019.1273&file=rspb20191273supp4.csv",
              paste0(h_dir,"/",h_supp_data), 
              method = "curl", 
              extra = paste0("-L -H ", '"User-Agent: Mozilla ', date(), '"'))

# read data and filter for second mandibular molars of caprines
#
h_comp_data_ <- read_csv(paste0(h_dir,"/",h_supp_data)) %>% 
  filter(element == "M/2" & taxon == "caprine")


#### Prepare metadata ####
#

# columns to keep for metadata
#
metadata_cols <- c("specimen",
                   "site",
                   "element",
                   "symmetry",
                   "individual",
                   "taxon",
                   "period")

new_metadata <- study_data_ %>% left_join(
  dplyr::select(new_metadata_, c("specimen",
                          "element",
                          "symmetry",
                          "individual")), by = "specimen") %>% 
  dplyr::select(all_of(metadata_cols)) %>% 
  unique()

# subset h metadata
#
h_comp_metadata <- h_comp_data_ %>%
  dplyr::group_by(specimen) %>% 
  dplyr::mutate(individual = cur_group_id()) %>% 
  dplyr::select(all_of(metadata_cols)) %>% 
  unique()

# prepare vm metadata from source publication and subset 
#
vm_comp_metadata <- vm_comp_data_ %>% 
  dplyr::mutate(element = "M/2", 
                symmetry = "left",
                taxon = "sheep",
                period = case_when(
                  site == "Kent" ~ "Late Bronze Age",
                  site == "Turgen" ~ "Late Bronze Age/Iron Age")
            ) %>% 
  dplyr::group_by(specimen) %>% 
  dplyr::mutate(individual = cur_group_id()) %>% 
  dplyr::select(all_of(metadata_cols)) %>% 
  unique() %>% 
  dplyr::mutate(specimen = as.character(specimen))

# merge metadata and create comparative column
#
comp_metadata <- base::rbind(vm_comp_metadata, h_comp_metadata) %>% 
  dplyr::mutate(comparative = "TRUE") %>% 
  rbind(unique(mutate(new_metadata, comparative = "FALSE")))


#### Prepare comparative isotope data file ####

# columns to keep for data
#
data_cols <- c("increment",
               "specimen",
               "d13C",
               "d18O",
               "measure")

study_data <- study_data_ %>% 
  dplyr::select(all_of(data_cols)) %>% 
  dplyr::filter(!row_number() %in% c(152, 163))  
  #                                  ^    ^
  # Remove final measurement for IBX03-04 due to sharp change in values near root-enamel junction

# combine comparative data
#
all_data <- data.table::rbindlist(list(study_data, h_comp_data_, vm_comp_data_), fill=T) %>% 
  dplyr::select(all_of(data_cols)) %>% 
  tibble::as_tibble()

#### Output ####
#

# directory to hold final files
#
out_data_dir <- "data/input/"

# outputs individual CSVs for each tooth sequence
# group_by(comp_data, specimen) %>%
#   do(write_csv(., paste0(out_data_dir, "isodata/", unique(.$specimen), "_C_O_meas.csv")))

# makes list of tables
# comp_data %>% group_split(specimen)

# outputs one CSV for all data
data.table::fwrite(all_data,
       paste0(out_data_dir,"all_data.csv"))

data.table::fwrite(comp_metadata,
          paste0(out_data_dir,"specimen.csv"))
