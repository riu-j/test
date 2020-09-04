# import library
library("tidyverse")
library("genetics")

# import source code


# set working directory
{
  drive <- getwd() %>% substring(1, 1)
  if(drive=="D"){
    #Desktop PC
    setwd("D:/Users/Ryota/OneDrive - 千葉大学/Project/src_r")
  }
  if(drive=="C"){
    #Note PC
    setwd("C:/Users/Ryota Jin/OneDrive - 千葉大学/Project/src_r")
  }
}

# define global variable
proj <- "s139"
process_time <- Sys.time() %>% substring(3, 19) %>% gsub("-", "", .) %>% gsub(":", "", .) %>% gsub(" ", "_", .)

# process
source("01_datamaker.R", encoding="UTF-8")
source("02_id_selection.R", encoding="UTF-8")
source("03_HWE_test.R", encoding="UTF-8")
source("04_converter_mean_2.0.R", encoding="UTF-8")
source("05_nmctl.R", encoding="UTF-8")
source("06_scatter_matrix.R", encoding="UTF-8")
source("07_spaghetti.R", encoding="UTF-8")

