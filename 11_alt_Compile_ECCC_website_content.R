
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
  mutate(
    # region = ifelse(region_type == "bcr",
    #                      paste0("BCR_",region),
    #                      region),
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
   filter(!Species_ID_Espèce %in% c("RBME_EAS","RBME_WES",
                                   "CORE","HORE")) %>%
  select(Species_ID_Espèce, BBS_Number__Numéro_BBS, Sort_Order__Ordre_de_tri) %>%
  rename_with(.,.fn = ~paste0(.x,"_core")) %>%
  distinct() %>%
  mutate(aou = as.integer(BBS_Number__Numéro_BBS_core),
         aou = ifelse(aou == 5280, 5275,aou),
         aou = ifelse(aou == 4641, 4642,aou),
         aou = ifelse(Species_ID_Espèce_core == "AHGU", 510,aou),
         aou = ifelse(Species_ID_Espèce_core == "WAVI", 6270,aou))


avian_core_table <- readxl::read_xlsx("data/Avian_Core_20251124.xlsx") %>%
  filter(!Species_ID_Espèce %in% c("RBME_EAS","RBME_WES",
                                   "CORE","HORE")) %>%
  select(Species_ID_Espèce, Full_Species__Espèce_complète, BBS_Number__Numéro_BBS,
         English_Name__Nom_Anglais,French_Name__Nom_Français,
         Scientific_Name__Nom_Scientifique, Sort_Order__Ordre_de_tri,
         Family_CmName_EN__Nom_Famille_Anglais,Family_CmName_FR__Nom_Famille_Français) %>%
  #rename_with(.,.fn = ~paste0(.x,"_core")) %>%
  distinct() %>%
  mutate(aou = as.integer(BBS_Number__Numéro_BBS),
         aou = ifelse(aou == 5280, 5275,aou),
         aou = ifelse(aou == 4641, 4642,aou),
         aou = ifelse(Species_ID_Espèce == "AHGU", 510,aou),
         aou = ifelse(Species_ID_Espèce == "WAVI", 6270,aou))



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


web_sp <- web %>%
  select(species,espece,bbs_num) %>%
  distinct() %>%
  rename(sp = species,
         es = espece)


BBS_AVIANCORE <- avian_core_table %>%
  select(Sort_Order__Ordre_de_tri,
         Full_Species__Espèce_complète,
         Species_ID_Espèce,
         English_Name__Nom_Anglais,
         French_Name__Nom_Français,
         Scientific_Name__Nom_Scientifique,
         Family_CmName_EN__Nom_Famille_Anglais,
         Family_CmName_FR__Nom_Famille_Français,
         aou) %>%
  filter(aou %in% names_match$bbs_num) %>%
  rename(Sort_Order = Sort_Order__Ordre_de_tri,
         Full_Species = Full_Species__Espèce_complète,
         Species_ID = Species_ID_Espèce,
         English_Name = English_Name__Nom_Anglais,
         French_Name = French_Name__Nom_Français,
         Scientific_Name = Scientific_Name__Nom_Scientifique,
         Family_CmName_EN = Family_CmName_EN__Nom_Famille_Anglais,
         Family_CmName_FR = Family_CmName_FR__Nom_Famille_Français,
         BBS_Number = aou) %>%
  mutate(Full_Species = ifelse(Full_Species == "Yes - Oui",
                               "Yes","No")) %>%
  inner_join(web_sp, by = c("BBS_Number" = "bbs_num")) %>%
  mutate(English_Name = ifelse(grepl(pattern = "identif",x = English_Name),sp,English_Name),
         French_Name = ifelse(grepl(pattern = "identif",x = French_Name),es,French_Name),
         English_Name = ifelse(!grepl(pattern = "(",x = English_Name, fixed = TRUE) &
                                 Full_Species == "No",sp,English_Name),
         English_Name = ifelse(grepl(pattern = "Sapsuckers (",x = English_Name, fixed = TRUE),
                               "Sapsuckers (Yellow-bellied/Red-naped/Red-breasted)",English_Name),
         French_Name = ifelse(!grepl(pattern = "(",x = French_Name, fixed = TRUE) &
                                Full_Species == "No",es,French_Name),
         test_en = ifelse(English_Name == sp,TRUE,FALSE),
         test_fr = ifelse(French_Name == es,TRUE,FALSE))

en_fail <- BBS_AVIANCORE %>%
  filter(!test_en)

fr_fail <- BBS_AVIANCORE %>%
  filter(!test_fr)

BBS_AVIANCORE <- BBS_AVIANCORE %>%
  select(-c(sp,es,test_en,test_fr,Full_Species))

write_excel_csv(BBS_AVIANCORE,paste0("website/bbs_aviancore_",YYYY+1,".csv"))

bbs_sp_repl <- BBS_AVIANCORE %>%
  select(BBS_Number,English_Name,French_Name) %>%
  distinct()

web <- web %>%
  inner_join(bbs_sp_repl,
             by = c("bbs_num" = "BBS_Number")) %>%
  mutate(species = English_Name,
         espece = French_Name) %>%
  select(-c(English_Name,French_Name))


families <- BBS_AVIANCORE %>%
  select(Sort_Order,Family_CmName_EN,Family_CmName_FR) %>%
  group_by(Family_CmName_EN,Family_CmName_FR) %>%
  summarise(sort_min = min(Sort_Order)) %>%
  arrange(sort_min) %>%
  ungroup() %>%
  mutate(familyID = row_number()) %>%
  rename(familyNameE = Family_CmName_EN,
         familyNameF = Family_CmName_FR) %>%
  select(familyID,familyNameE,familyNameF) %>%
  mutate(sortE = familyID,
         sortF = familyID)



write_excel_csv(families,paste0("website/bbs_families_",YYYY+1,".csv"))




table(web$species, useNA = "always")

# Table of regions --------------------------------------------------------

regions_w_trends <- web %>%
  select(region,region_type) %>%
  distinct() %>%
  arrange(region_type,region)

states <- rnaturalearth::ne_states(country = c("Canada","United States of America")) %>%
  sf::st_drop_geometry() %>% select(postal,name_fr) %>%
  distinct() %>%
  rename(statprov_name_fr = name_fr)

# countries <- rnaturalearth::ne_countries() %>%
#   sf::st_drop_geometry() %>%
#   filter(admin %in% c("Canada","United States of America")) %>%
#   select(postal, name_fr)
#
# country_states <- countries %>% bind_rows(states)

bcrs <- sf::read_sf("data/bcr_2025_lakes12.gpkg") %>%
  mutate(region_type = "bcr",
         region = paste0("BCR_",bcr_label),
         geo.area = region,
         region_name_en = bcr_name_en,
         region_name_fr = bcr_name_fr) %>%
  sf::st_drop_geometry() %>%
  select(region_type,geo.area,region,region_name_en,region_name_fr) %>%
  distinct()



strats <- sf::read_sf("data/bcr_2025_lakes12_statprov3.gpkg") %>%
  filter(statprov_name != "Newfoundland") %>%
  left_join(states, by = c("statprov_code" = "postal")) %>%
  rowwise() %>%
  mutate(country_name = ifelse(country_name == "United States",
                               "United States of America",
                               country_name),
         country_name_fr = ifelse(country_name == "United States of America",
                               "États-Unis",
                               "Canada"),
         region_type = "stratum",
         region = paste0(country_code,"-",statprov_code,"-",bcr_label),
         geo.area = region,
         region_name_en = paste0(country_name,"-",statprov_name,"-",bcr_name_en),
         region_name_fr = paste0(country_name_fr,"-",statprov_name_fr,"-",bcr_name_fr)) %>%
  sf::st_drop_geometry() %>%
  select(region_type,geo.area,region,region_name_en,region_name_fr) %>%
  distinct()


bcr_country <- sf::read_sf("data/bcr_2025_lakes12_statprov3.gpkg") %>%
  left_join(states, by = c("statprov_code" = "postal")) %>%
  rowwise() %>%
  mutate(country_name = ifelse(country_name == "United States",
                               "United States of America",
                               country_name),
         country_name_fr = ifelse(country_name == "United States of America",
                                  "États-Unis",
                                  "Canada"),
         region_type = "bcr_by_country",
         region = paste0(country_name,"-","BCR_",bcr_label),
         region_name_en = paste0(country_name,"-",bcr_name_en),
         region_name_fr = paste0(country_name_fr,"-",bcr_name_fr),
         geo.area = ifelse(country_name == "United States of America",
                           paste0("BCR_",bcr_label,"U"),
                           paste0("BCR_",bcr_label,"C"))) %>%
  sf::st_drop_geometry() %>%
  select(region_type,geo.area,region,region_name_en,region_name_fr) %>%
  distinct()

countries <- sf::read_sf("data/bcr_2025_lakes12_statprov3.gpkg")  %>%
  rowwise() %>%
  mutate(country_name = ifelse(country_name == "United States",
                               "United States of America",
                               country_name),
         country_name_fr = ifelse(country_name == "United States of America",
                                  "États-Unis",
                                  "Canada"),
         region_type = "country",
         region = paste0(country_name),
         geo.area = ifelse(region == "United States of America",
                           "UU","CC"),
         region_name_en = paste0(country_name),
         region_name_fr = paste0(country_name_fr)) %>%
  filter(country_code %in% c("US","CA")) %>%
  sf::st_drop_geometry() %>%
  select(region_type,geo.area,region,region_name_en,region_name_fr) %>%
  distinct()



prov_state <- sf::read_sf("data/bcr_2025_lakes12_statprov3.gpkg")  %>%
  left_join(states, by = c("statprov_code" = "postal")) %>%
  rowwise() %>%
  mutate(country_name = ifelse(country_name == "United States",
                               "United States of America",
                               country_name),
         country_name_fr = ifelse(country_name == "United States of America",
                                  "États-Unis",
                                  "Canada"),
         region_type = "prov_state",
         region = paste0(statprov_code),
         geo.area = region,
         region_name_en = paste0(statprov_name),
         region_name_fr = paste0(statprov_name_fr)) %>%
  filter(country_code %in% c("US","CA"),
         region_name_en != "Newfoundland") %>%
  sf::st_drop_geometry() %>%
  select(region_type,geo.area,region,region_name_en,region_name_fr) %>%
  distinct()

survey_w <- data.frame(region_type = "survey-wide",
                       region = "Survey-wide",
                       geo.area = "SW",
                       region_name_en = "Survey-wide",
                       region_name_fr = "Zone complète du relevé")

regional_prefixes_sorts <- readxl::read_excel("data/State_prefixes.xlsx")

regions <- bind_rows(survey_w,
                     countries,
                     prov_state) %>%
  left_join(regions_w_trends,
            by = c("region","region_type")) %>%
  distinct()

regions_link <- regional_prefixes_sorts  %>%
  mutate(geo.area = State_code) %>%
  left_join(regions,
            by = c("geo.area")) %>%
  mutate(sortE = sort,
         sortF = sort) %>%
  rename(provID = sort,
         shortCode = State_code,
         descriptionE = region_name_en,
         descriptionF = region_name_fr,
         prefixF = prefix) %>%
  select(provID,
           shortCode,
           descriptionE,
           descriptionF,
           sortE,
           sortF,
           prefixF,
         geo.area,
         region_type,
         region) |>
  arrange(provID)

regions <- regions_link %>%
  select(-c(geo.area,region_type,region))
write_excel_csv(regions,paste0("website/bbs_province_",YYYY+1,".csv"))


bcr_sorts <- readxl::read_excel("data/BCR_prefixes.xlsx")

sub_0 <- function(x){
  dig <- paste0(0,str_extract(x,
                     "[[:digit:]]"))
  new <- str_replace(x,
              pattern = "[[:digit:]]",
              replacement = dig)
}
bcr_regions <- bind_rows(bcrs,
                     bcr_country) %>%
  distinct() %>%
  inner_join(regions_w_trends,
             by = c("region_type","region")) %>%
  mutate(shortCode = ifelse(str_detect(geo.area,
                                       pattern = "_[[:digit:]]{1}$|_[[:digit:]]{1}[[:alpha:]]{1}"),
                            sub_0(geo.area),
                            geo.area),
         shortCode = str_replace(shortCode,
                                 "_",
                                 ""))

bcr_regions_link <- bcr_sorts %>%
  left_join(bcr_regions,
            by = c("shortCode")) %>%
  select(bcrID,
         shortCode,
         bcrE,
         bcrF,
         bcrNameE,
         bcrNameF,
         sortE,
         sortF,
         prefixF,
         geo.area,
         region_type,
         region) |>
  arrange(bcrNameE) |>
  ungroup() |>
  mutate(sortE = row_number(bcrNameE))  |>
  arrange(bcrNameF) |>
  ungroup() |>
  mutate(sortF = row_number(bcrNameF)) %>%
  arrange(bcrID)

bcr_regions <- bcr_regions_link %>%
  select(-c(geo.area,region_type,region))

write_excel_csv(bcr_regions,paste0("website/bbs_bcr_",YYYY+1,".csv"))






strata_sorts <- read_csv("data/bbs_asr_lookups_20260514.csv") |>
  select(-c(9:10)) |>
  mutate(shortCode = asrID)


strat_regions_link <- strats %>%
  distinct() %>%
  inner_join(regions_w_trends,
             by = c("region_type","region")) |>
  select(-c(region_name_en,region_name_fr)) |>
  mutate(asrID = geo.area) |>
  full_join(strata_sorts,
            by = c("asrID")) |>
  arrange(asrNameE) |>
  ungroup() |>
  mutate(sortE = row_number(asrNameE))  |>
  arrange(asrNameF) |>
  ungroup() |>
  mutate(sortF = row_number(asrNameF))


strat_regions <- strat_regions_link |>
  select(-c(geo.area,region_type,region))


write_excel_csv(strat_regions,paste0("website/bbs_asr_",YYYY+1,".csv"))


regions <- bind_rows(regions_link,
                     bcr_regions_link,
                     strat_regions_link)

saveRDS(regions,"required_data/regions_link_translate.rds")

regs_repl <- regions %>%
  select(region,shortCode,region_type)

web <- web %>%
  left_join(regs_repl, by = c("region","region_type"))

regs <- unique(regions$region)
regs_web <- unique(web$region)

if( any(!regs_web %in% regs)){
  paste("missing",regs_web[which(!regs_web %in% regs)])
  stop("missing regions from the regional tables for website")
}


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

if(!dir.exists(paste0("website/WebMaps/"))){
  dir.create(paste0("website/WebMaps/"))
}
#loading base maps for the website plots
canmap <- bbsBayes2::load_map("bbs")


## looping through species to create the maps in folder webmaps
sp_loop <- unique(web$bbs_num)

n_files <- 0
j = 1
dir.create(paste0("website/WebMaps/",j))

for(sp in sp_loop){

  dfta <- web %>%
    filter(bbs_num == sp)

  if(n_files + nrow(dfta) > 4999){
    j <- j+1
    n_files <- 0
    dir.create(paste0("website/WebMaps/",j))
    dir_map_tmp <- paste0("website/WebMaps/",j,"/")
  }else{
    dir_map_tmp <- paste0("website/WebMaps/",j,"/")
    n_files <- n_files+nrow(dfta)
  }

  dft <- web %>%
    filter(bbs_num == sp,
           trend_time == "Long-term")

if(!all(file.exists(paste0(paste0(dir_map_tmp,dft$mapfile))))){
  generate_web_maps(dft,
                    canmap = canmap,
                    basemap = NULL,
                    n_cores = 6)
}

  dft <- web %>%
    filter(bbs_num == sp,
           trend_time == "Short-term")

if(!all(file.exists(paste0(paste0(dir_map_tmp,dft$mapfile))))){
  generate_web_maps(dft,
                    canmap = canmap,
                    basemap = NULL,
                    n_cores = 6)
}
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


# temporary reorder maps into smaller folders
#
# library(fs)
# sp_loop <- unique(web$bbs_num)
#
# n_files <- 0
# j = 1
# dir.create(paste0("website/WebMaps/",j))
#
# for(sp in sp_loop){
#
#   dfta <- web %>%
#     filter(bbs_num == sp)
#
#   if(n_files + nrow(dfta) > 9999){
#     j <- j+1
#     n_files <- nrow(dfta)
#     dir.create(paste0("website/WebMaps/",j))
#     dir_map_tmp <- paste0("website/WebMaps/",j,"/")
#   }else{
#     dir_map_tmp <- paste0("website/WebMaps/",j,"/")
#     n_files <- n_files+nrow(dfta)
#   }
#
#   files_cut <- paste0("website/WebMaps/",dfta$mapfile)
#   files_paste <- paste0(dir_map_tmp,dfta$mapfile)
#
#
# file_move(files_cut,files_paste)
#
# print(paste(n_files,round(which(sp_loop == sp)/length(sp_loop),2)))
# }



# Export final csv files --------------------------------------------------


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
          "mapfile",
          "shortCode")

clnms = c("sp","species","espece","geo.area","trendtype",
          "startyear","endyear","trend",
          "llimit","ulimit","reliab.over",
          "reliab.prec","reliab.cov","reliab.pool",
          "p.decrease","p.increase","p.d50","pd50.25",
          "pd25.0","pi0.33","pi33.100","pi100",
          "percent.change","percent.change.llimit",
          "percent.change.ulimit","nroutesduringtrend",
          "strata.inc","st.excl.long","mapfile","shortCode")

if(any(!clout %in% names(web))){
  print(clout[which(!clout %in% names(web))])
}

web = web[,clout]
names(web) = clnms


readr::write_excel_csv(web, paste0("website/",YYYY," BBS trends for website w US.csv"))











# Indices for website -----------------------------------------------------


webi_short <- indices %>%
  inner_join(bbs_sp_repl,
             by = c("bbs_num" = "BBS_Number")) %>%
  mutate(species = English_Name,
         espece = French_Name) %>%
  select(-c(English_Name,French_Name))%>%
  mutate(region_type = as.character(region_type),
    region = ifelse(region_type == "bcr",
                         paste0("BCR_",region),
                         region),
    region = ifelse(region == "continent", "Survey-wide",region),
    region_type = ifelse(region_type == "continent", "survey-wide",region_type)) %>%
  left_join(regs_repl, by = c("region","region_type")) %>%
  filter(year > (YYYY-11),
         trend_time == "Short-term")

webi <- indices %>%
  inner_join(bbs_sp_repl,
             by = c("bbs_num" = "BBS_Number")) %>%
  mutate(species = English_Name,
         espece = French_Name) %>%
  select(-c(English_Name,French_Name)) %>%
  mutate(region_type = as.character(region_type),
         region = ifelse(region_type == "bcr",
                         paste0("BCR_",region),
                         region),
         region = ifelse(region == "continent", "Survey-wide",region),
         region_type = ifelse(region_type == "continent", "survey-wide",region_type)) %>%
  left_join(regs_repl, by = c("region","region_type")) %>%
  filter(year > 1969,
         trend_time == "Long-term") %>%
  bind_rows(.,webi_short)


clouti =  c("bbs_num",
            "species",
            "espece",
            "region",
            "trend_time",
            "year",
            "index",
            "index_q_0.05",
            "index_q_0.95",
            "shortCode")
clnmsi = c("sp","species","espece","geo.area","trendtype",
           "year","an.index",
           "llimit","ulimit","shortCode")


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

table(webi$geo.area,useNA = "always")











