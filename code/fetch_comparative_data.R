# This script will download comparative isotope data from 
# Ventresca-Miller et al. 2020, doi:10.1080/20548923.2020.1759316
# Hermes et al. 2019, doi:10.1098/rspb.2019.1273
# and generate two CSV files to hold 1) data and 2) metadata (linked by specimen ID)

library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(docxtractr)
library(tibble)
library(readr)

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


#### Generate metadata ####
#

# columns to keep for metadata
#
metadata_cols <- c("specimen",
                   "site",
                   "element",
                   "symmetry",
                   "taxon",
                   "period")

# subset h metadata
#
h_comp_metadata <- h_comp_data_ %>% 
  select(all_of(metadata_cols)) %>% 
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
  select(all_of(metadata_cols)) %>% 
  unique()

# merge metadata and create comparative column
#
comp_metadata <- rbind(vm_comp_metadata, h_comp_metadata) %>% 
  mutate(comparative = "TRUE")


#### Prepare comparative data file ####

# columns to keep for data
#
data_cols <- c("increment",
               "specimen",
               "d13C",
               "d18O",
               "measure")

# combine data
#
comp_data <- rbindlist(list(h_comp_data_, vm_comp_data_), fill=T) %>% 
  select(all_of(data_cols)) %>% 
  as_tibble()


#### Output ####
#

# directory to hold final files
#
data_dir <- "data/comparative"

write_csv(comp_data,
          paste0(data_dir,"/all_comparative_data.csv"))

write_csv(comp_metadata,
          paste0(data_dir,"/all_comparative_metadata.csv"))
