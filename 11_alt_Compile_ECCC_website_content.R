
## compile trend and index files for
### 1 - website
YYYY <- 2024

webmaps <- FALSE # set to true if needing to create all map images for ECCC website

library(bbsBayes2)
library(tidyverse)

source("functions/mapping.R")
source("functions/loess_func.R")
# custom functions to calculate reliability categories and determine website inclusion


sp_list <- readRDS("sp_list_w_generations.rds") %>%
  filter(model == TRUE)

trends <- readRDS(paste0("Website/All_BBS_Trends_",YYYY,".rds"))
indices <- readRDS(paste0("Website/All_BBS_Full_Indices_",YYYY,".rds"))
indices_smooth <- readRDS(paste0("Website/All_BBS_Smoothed_Indices_",YYYY,".rds"))


# Website trends ----------------------------------------------------------


web <- trends %>%
  filter(trend_time != "Three-generation",
         !(region == "BCR7" & region_type == "prov_state")) %>%
  mutate(region = ifelse(region_type == "bcr",
                         paste0("BCR_",region),
                         region),
         region = ifelse(region == "continent", "Survey-wide",region),
         region_type = as.character(region_type),
         region_type = ifelse(region_type == "continent", "survey-wide",region_type),
         prob_decrease_0_25_percent = prob_decrease_0_percent-prob_decrease_25_percent,
         prob_decrease_25_50_percent = prob_decrease_0_percent - (prob_decrease_0_25_percent + prob_decrease_50_percent),
         prob_increase_0_33_percent = prob_increase_0_percent-prob_increase_33_percent,
         prob_increase_33_100_percent = prob_increase_0_percent - (prob_increase_0_33_percent + prob_increase_100_percent),
         mapfile = paste(bbs_num,region,trend_time,"map.png",sep = "_"),
         strata_included = paste(strata_included,strata_excluded,sep = " ; "),
         strata_excluded = "")





# web_species <- read.csv("required_data/BBS_AvianCore.csv")

avian_core <- readxl::read_xlsx("data/Avian_Core_20251124.xlsx") %>%
  filter(Full_Species__Espèce_complète == "Yes - Oui" |
           Species_ID_Espèce == "WAVI") %>%
  select(Species_ID_Espèce, BBS_Number__Numéro_BBS, Sort_Order__Ordre_de_tri) %>%
  rename_with(.,.fn = ~paste0(.x,"_core")) %>%
  distinct() %>%
  mutate(aou = as.integer(BBS_Number__Numéro_BBS_core),
         aou = ifelse(aou == 5280, 5275,aou),
         aou = ifelse(aou == 4641, 4642,aou),
         aou = ifelse(Species_ID_Espèce_core == "AHGU", 510,aou),
         aou = ifelse(Species_ID_Espèce_core == "WAVI", 6270,aou))



names_match <- web %>%
  select(species,espece,bbs_num) %>%
  distinct()%>%
  left_join(avian_core, by = c("bbs_num" = "aou"))

sp_not_on_website <- names_match %>%
  filter(is.na(Species_ID_Espèce_core))

paste(sp_not_on_website$species,collapse = "; ")
### check to make sure species have not been excluded due to
### out of sync changes to core and BBS species names

names_match <- web %>%
  select(species,espece,bbs_num) %>%
  distinct()%>%
  left_join(avian_core, by = c("bbs_num" = "aou")) %>%
  filter(!is.na(Species_ID_Espèce_core))

web <- web %>%
  filter(bbs_num %in% names_match$bbs_num)


# Table of regions --------------------------------------------------------

regions_w_trends <- web %>%
  select(region,region_type) %>%
  distinct() %>%
  arrange(region_type,region)

bcrs <- sf::read_sf("data/bcr_2025_lakes12.gpkg") %>%
  mutate(region_type = "bcr",
         region = paste0("BCR_",bcr_label),
         region_name_en = bcr_name_en,
         region_name_fr = bcr_name_fr) %>%
  sf::st_drop_geometry() %>%
  select(region_type,region,region_name_en,region_name_fr) %>%
  distinct()



strats <- sf::read_sf("data/bcr_2025_lakes12_statprov3.gpkg") %>%
  filter(statprov_name != "Newfoundland") %>%
  mutate(country_name = ifelse(country_name == "United States",
                               "United States of America",
                               country_name),
         region_type = "stratum",
         region = paste0(country_code,"-",statprov_code,"-",bcr_label),
         region_name_en = paste0(country_name,"-",statprov_name,"-",bcr_name_en),
         #region_name_fr = paste0(country_name,"-",statprov_name,"-",bcr_name_fr),
         region_name_fr = "") %>%
  sf::st_drop_geometry() %>%
  select(region_type,region,region_name_en,region_name_fr) %>%
  distinct()


bcr_country <- sf::read_sf("data/bcr_2025_lakes12_statprov3.gpkg") %>%
  mutate(country_name = ifelse(country_name == "United States",
                               "United States of America",
                               country_name),
         region_type = "bcr_by_country",
         region = paste0(country_name,"-","BCR_",bcr_label),
         region_name_en = paste0(country_name,"-",bcr_name_en),
         #region_name_fr = paste0(country_name,"-",bcr_name_fr),
         region_name_fr = "") %>%
  sf::st_drop_geometry() %>%
  select(region_type,region,region_name_en,region_name_fr) %>%
  distinct()

countries <- sf::read_sf("data/bcr_2025_lakes12_statprov3.gpkg") %>%
  mutate(country_name = ifelse(country_name == "United States",
                               "United States of America",
                               country_name),
         region_type = "country",
         region = paste0(country_name),
         region_name_en = paste0(country_name),
         #region_name_fr = paste0(country_name),
         region_name_fr = "") %>%
  filter(country_code %in% c("US","CA")) %>%
  sf::st_drop_geometry() %>%
  select(region_type,region,region_name_en,region_name_fr) %>%
  distinct()


prov_state <- sf::read_sf("data/bcr_2025_lakes12_statprov3.gpkg") %>%
  mutate(country_name = ifelse(country_name == "United States",
                               "United States of America",
                               country_name),
         region_type = "prov_state",
         region = paste0(statprov_code),
         region_name_en = paste0(statprov_name),
         #region_name_fr = paste0(statprov_name),
         region_name_fr = "") %>%
  filter(country_code %in% c("US","CA"),
         region_name_en != "Newfoundland") %>%
  sf::st_drop_geometry() %>%
  select(region_type,region,region_name_en,region_name_fr) %>%
  distinct()

survey_w <- data.frame(region_type = "survey-wide",
                       region = "Survey-wide",
                       region_name_en = "Survey-wide",
                       region_name_fr = "")

regions <- bind_rows(survey_w,
                     countries,
                     prov_state,
                     bcrs,
                     bcr_country,
                     strats) %>%
  distinct() %>%
  inner_join(regions_w_trends,
             by = c("region_type","region"))

write_excel_csv(regions,"website/list_of_regions.csv")


# miss_bbs_num <- names_match %>%
#   select(bbs_num,species) %>%
#   left_join(.,
#             web_species,
#             by = c("bbs_num" = "bbsNumber"),
#             multiple = "all") %>%
#   arrange(bbs_num) %>%
#   filter(is.na(commonNameE))
#
# if(nrow(miss_bbs_num) > 0){
#   warning("At least one bbs number is missing from Avian Core")
#
#   print(paste("Avian core is missing",
#               paste(miss_bbs_num$bbs_num,
#                     collapse = ", ")))
#   web <- web %>%
#     filter(bbs_num %in% web_species$bbsNumber)
# }
#
#
# miss_english_names <- names_match %>%
#   select(bbs_num,species,espece) %>%
#   left_join(.,
#             web_species,
#             by = c("species" = "commonNameE"),
#             multiple = "all") %>%
#   arrange(bbs_num) %>%
#   filter(is.na(bbs_num))
#
#
# if(nrow(miss_english_names) > 0){
#   warning("At least one species name is missing from Avian Core")
#   web <- web %>%
#     filter(species %in% web_species$commonNameE)
#
# }

# generate maps for CWS website -------------------------------------------

if(webmaps){


#loading base maps for the website plots
canmap <- bbsBayes2::load_map("bbs")


## looping through species to create the maps in folder webmaps
sp_loop <- unique(web$bbs_num)
for(sp in sp_loop){
  dft <- web %>%
    filter(bbs_num == sp,
           trend_time == "Long-term")


  generate_web_maps(dft,
                    canmap = canmap,
                    basemap = basemap)
  dft <- web %>%
    filter(bbs_num == sp,
           trend_time == "Short-term")
  generate_web_maps(dft,
                    canmap = canmap,
                    basemap = basemap)

  print(round(which(sp_loop == sp)/length(sp_loop),2))
}
## maps have Canadian regions only and show the regions included in each trend

map_dup_test <- any(duplicated(web$mapfile))

test_map <- any(!file.exists(paste0("website/webmaps/",web$mapfile)))

if(test_map){
  map_test <- paste0("website/webmaps/",web$mapfile)

  w_miss <- map_test[which(!file.exists(map_test))]

  stop("At least one map is missing")

}

maps <- list.files("website/webmaps/",
                   pattern = ".png")
test_map_extra <- any((web$mapfile %in% maps) == FALSE)

if(test_map_extra){
  warning("There are extra maps in the webmaps folder")
}

}#end webmaps

clout = c("bbs_num",
          "species",
          "espece",
          "region",
          "trend_time",
          "start_year",
          "end_year",
          "trend",
          "trend_q_0.025",
          "trend_q_0.975",
          "reliability",
          "width_of_95_percent_credible_interval",
          "reliab.cov",
          "backcast_flag",
          "prob_decrease_0_percent",
          "prob_increase_0_percent",
          "prob_decrease_50_percent",
          "prob_decrease_25_50_percent",
          "prob_decrease_0_25_percent",
          "prob_increase_0_33_percent",
          "prob_increase_33_100_percent",
          "prob_increase_100_percent",
          "percent_change",
          "percent_change_q_0.025",
          "percent_change_q_0.975",
          "n_routes",
          "strata_included",
          "strata_excluded",
          "mapfile")

clnms = c("sp","species","espece","geo.area","trendtype",
          "startyear","endyear","trend",
          "llimit","ulimit","reliab.over",
          "reliab.prec","reliab.cov","reliab.pool",
          "p.decrease","p.increase","p.d50","pd50.25",
          "pd25.0","pi0.33","pi33.100","pi100",
          "percent.change","percent.change.llimit",
          "percent.change.ulimit","nroutesduringtrend",
          "strata.inc","st.excl.long","mapfile")

if(any(!clout %in% names(web))){
  print(clout[which(!clout %in% names(web))])
}

web = web[,clout]
names(web) = clnms


readr::write_excel_csv(web, paste0("website/",YYYY," BBS trends for website w US.csv"))











# Indices for website -----------------------------------------------------



webi_short <- indices %>%
  mutate(region_type = as.character(region_type),
    region = ifelse(region_type == "bcr",
                         paste0("BCR_",region),
                         region)) %>%
  filter(year > (YYYY-11),
         trend_time == "Short-term")

webi <- indices %>%
  mutate(region = ifelse(region_type == "bcr",
                         paste0("BCR_",region),
                         region)) %>%
  filter(year > 1969,
         trend_time == "Long-term") %>%
  bind_rows(.,webi_short)

webi <- webi %>%
  filter(bbs_num %in% names_match$bbs_num) %>%
  mutate(region = ifelse(region == "continent", "Survey-wide",region),
         region_type = ifelse(region_type == "continent", "survey-wide",region_type))

clouti =  c("bbs_num",
            "species",
            "espece",
            "region",
            "trend_time",
            "year",
            "index",
            "index_q_0.05",
            "index_q_0.95")
clnmsi = c("sp","species","espece","geo.area","trendtype",
           "year","an.index",
           "llimit","ulimit")


if(any(!clouti %in% names(webi))){
  print(clouti[which(!clouti %in% names(webi))])
}


webi = webi[,clouti]
names(webi) <- clnmsi



tshort = web %>%
  filter(trendtype == "Short-term")
tlong= web %>%
  filter(trendtype == "Long-term")
index_short <- webi %>%
  filter(year == YYYY-10,
         trendtype == "Short-term")
index_long <- webi %>%
  filter(year == YYYY-10,
         trendtype == "Long-term")

if(nrow(index_long) != nrow(tlong) |
   nrow(index_short) != nrow(tshort)){
warning("The number of indices and trends don't match \n explore index_trend_test")

index_trend_test <- webi %>%
  filter(year == YYYY-10) %>%
  select(species,trendtype,geo.area,an.index) %>%
  full_join(.,web,
            by = c("species","trendtype","geo.area"))
}

readr::write_excel_csv(webi, paste0("website/",YYYY," BBS indices for website w US.csv"))













