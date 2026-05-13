## Land composition ###

# libraries ---------------------------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(reshape2)
library(ggplot2)
library(stringr)
library(conflicted)
library(writexl)
library(openxlsx)
library(stats)
library(zoo)
library(cluster)
library(scales)

here()
conflicted::conflict_prefer("rename", "dplyr")
conflicted::conflict_prefer("mutate", "dplyr")
conflicted::conflict_prefer("summarise", "dplyr")
conflicts_prefer(dplyr::filter)
here()

grid50_equal_area <- read_csv("~/GitHub/DownscalingFABLE/Data/grid50_equal_area.csv")
EF_Pools_transition_Ecoregion <- read_excel("~/GitHub/DownscalingFABLE/Data/EF_Pools_transition_Ecoregion.xlsx")

df_region <-grid50_equal_area %>%
  select(id_c, ECO_NAME)

df <-df_region %>%
  left_join(EF_Pools_transition_Ecoregion) %>%
  distinct() %>%
  filter(ECO_NAME !="NA")

readr::write_csv(df, here::here( "EF_Pools_transition_Ecoregion.csv"))

