
### Generates the csv tables of trends and annual indices for sharing through the
### Google Drive and for upload into the NatureCounts database to support
###  - State of Canada's Birds, etc.
###
YYYY <- 2024


##################



webmaps <- FALSE # set to true if needing to create all map images for ECCC website

library(bbsBayes2)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(readxl)
#setwd("C:/GitHub/CWS_2023_BBS_Analyses")

external_dir <- "d:/BBS_Trends_CWS"
output_dir <- "d:/BBS_Trends_CWS/output"

source("functions/mapping.R")
source("functions/loess_func.R")
# custom functions to calculate reliability categories and determine website inclusion


sp_list <- readRDS("sp_list_w_generations.rds") %>%
  filter(model == TRUE)

#
# avian_core <- read_csv("data/ECCC Avian Core 20230601.csv") %>%
#   rename_with(.,.fn = ~paste0(.x,"_core")) %>%
#   mutate(aou = as.integer(BBS_Number_core))
#
# rep_aou <- avian_core %>% group_by(aou) %>% summarise(n = n()) %>% filter(n > 1, !is.na(aou))
# rep_core <- avian_core %>%
#   filter(aou %in% rep_aou$aou)
#
# avian_core <- avian_core %>%
#   filter(!(aou %in% rep_aou$aou & Full_Species_core == "No"))
#
# sp_list <- sp_list %>%
#   inner_join(.,avian_core,by = "aou")
#
# naturecounts_codes <- naturecounts::meta_species_codes() %>%
#   filter(authority == "BBS2") %>%
#   select(species_id2,species_code) %>%
#   mutate(aou = as.integer(species_code),
#          naturecounts_species_id = species_id2) %>%
#   distinct() %>%
#   select(-c(species_id2,species_code))
#
# sp_list <- sp_list %>%
#   left_join(.,naturecounts_codes,
#             by = "aou")
#
# tmp2 <- sp_list %>% filter(is.na(naturecounts_species_id)) %>%
#   ungroup() %>%
#   select(aou,english,BBS_Number_core,naturecounts_species_id)
# if(nrow(tmp2) > 0){
#   stop("Species don't match with nature counts")
# }



re_collect <- FALSE
# Compile all trends and indices ------------------------------------------------------

if(re_collect){
trends <- NULL
indices <- NULL
indices_smooth <- NULL

for(i in 1:nrow(sp_list)){


  aou <- as.integer(sp_list[i,"aou"])

  if(file.exists(paste0(external_dir,"/Indices/list_",aou,"_indices.rds"))){
  inds_1 <- readRDS(paste0(external_dir,"/Indices/list_",aou,"_indices.rds"))

  trends_1 <- readRDS(paste0(external_dir,"/Trends/",aou,"_trends.rds"))

  trends <- bind_rows(trends,trends_1)

  inds_1s <- inds_1 %>%
    filter(indices_type == "smooth")
  inds_1f <- inds_1 %>%
    filter(indices_type == "full")

  indices_smooth <- bind_rows(indices_smooth,inds_1s)
  indices <- bind_rows(indices,inds_1f)
}


  print(round(i/nrow(sp_list),2))



}

saveRDS(trends,"output/trends_collected.rds")
saveRDS(indices,"output/indices_collected.rds")
saveRDS(indices_smooth,"output/indices_smooth_collected.rds")

}else{


  trends <- readRDS("output/trends_collected.rds")%>%
    mutate(region = ifelse(region == "continent",
                           "Survey-wide",region),
           region_type = ifelse(region_type == "continent",
                                "survey-wide",region_type))
  indices <- readRDS("output/indices_collected.rds")%>%
    mutate(region = ifelse(region == "continent",
                           "Survey-wide",region),
           region_type = ifelse(region_type == "continent",
                                "survey-wide",region_type))
  indices_smooth <- readRDS("output/indices_smooth_collected.rds")%>%
    mutate(region = ifelse(region == "continent",
                           "Survey-wide",region),
           region_type = ifelse(region_type == "continent",
                                "survey-wide",region_type))

}


# Checking for long-term hind-casting ------------------------------------
#
# back_cast <- trends %>%
#   filter(#backcast_flag < 0.67,
#          trend_time == "Long-term",
#          region_type == "stratum")
#
# strat_back <- back_cast %>%
#   group_by(region,backcast_reliab) %>%
#   summarise(n_species = n(),
#             .groups = "drop") %>%
#   pivot_wider(id_cols = c(region),
#               names_from = backcast_reliab,
#               values_from = n_species) %>%
#   rowwise() %>%
#   mutate(Low = ifelse(is.na(Low),0,Low),
#          Medium = ifelse(is.na(Medium),0,Medium),
#          High = ifelse(is.na(High),0,High),
#          n_sp = sum(Low,Medium,High),
#          p_low = Low/n_sp)
#
# strat_to_bolster <- strat_back %>%
#   filter(p_low > 0.66) %>%
#   select(region) %>%
#   unlist()
#
# strat_to_bolster
# # region1    region2    region3
# # "CA-AB-6N" "CA-NT-6N" "CA-NT-7W"
# #
#
# sp_to_run_w_voronoi <- trends %>%
#   filter(region %in% strat_to_bolster) %>%
#   select(species,bbs_num) %>%
#   distinct()
#
# saveRDS(sp_to_run_w_voronoi,"species_to_consider_voronoi.rds")

# Compare to last year's trends -------------------------------------------

avian_core <- readxl::read_xlsx("data/Avian_Core_20251124.xlsx") %>%
  filter(Full_Species__Espèce_complète == "Yes - Oui") %>%
  select(Species_ID_Espèce, BBS_Number__Numéro_BBS, Sort_Order__Ordre_de_tri) %>%
  rename_with(.,.fn = ~paste0(.x,"_core")) %>%
  distinct() %>%
  mutate(aou = as.integer(BBS_Number__Numéro_BBS_core))

core_link <- sp_list %>%
  ungroup() %>%
  select(naturecounts_sort_order,aou,naturecounts_species_id) %>%
  left_join(avian_core, by = "aou")



lastyear = read_csv("data/All_BBS_trends_2023.csv")

ly_trends <- lastyear[,c("species","bbs_num","region","trend_time",
                         "n_strata_included","n_routes",
                         "trend",
                         "trend_q_0.05","trend_q_0.95",
                         "width_of_95_percent_credible_interval",
                         "reliab.cov","coverage")] %>%
  rename(trend_2023 = trend,
         trend_q_0.05_2023 = trend_q_0.05,
         trend_q_0.95_2023 = trend_q_0.95,
         number_of_strata_2023 = n_strata_included,
         number_of_routes_2023 = n_routes,
         CI_2023 = width_of_95_percent_credible_interval,
         reliab_cov_2023 = reliab.cov,
         coverage_2023 = coverage) %>%
  mutate(region = ifelse(region == "continent",
                         "Survey-wide",region)) |>
  filter(region %in% c("Survey-wide","Canada","United States of America")) %>%
  select(-c(species))



trends_comp <- trends %>%
  inner_join(.,ly_trends,
             by = c("bbs_num",
                    "region",
                    "trend_time")) %>%
  left_join(.,core_link,by = c("bbs_num" = "aou")) %>%
  mutate(diff_trend = trend - trend_2023,
         diff_coverage = reliab.cov - reliab_cov_2023) %>%
  rename(CI = width_of_95_percent_credible_interval)


comp_xy <- ggplot(data = trends_comp,
                  aes(x = trend_2023,
                      y = trend,
                      alpha = 1/CI))+
  geom_point()+
  geom_abline(intercept = 0,slope = 1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  facet_grid(cols = vars(trend_time),
             rows = vars(region),
             scales = "free")+
  scale_colour_viridis_d(begin = 0.1, end = 0.9)

pdf("trends_comparison.pdf",
    width = 11, height = 8.5)
print(comp_xy)

trends_comp_can <- trends_comp %>%
  filter(region == "Canada",
         trend_time == "Long-term") %>%
  mutate(non_overlap = ifelse((diff_trend > (trend-trend_q_0.05)) |
                            (diff_trend < (trend-trend_q_0.95)),
                          TRUE,FALSE))
sp_labs <- trends_comp_can %>%
  filter(non_overlap)

comp_xy <- ggplot(data = trends_comp_can,
                  aes(x = trend_2023,
                      y = trend,
                      alpha = 1/CI,
                      colour = non_overlap))+
  geom_point()+
  geom_errorbar(aes(ymin = trend_q_0.05,
                    ymax = trend_q_0.95,
                    colour = non_overlap),
                alpha = 0.2,
                width = 0)+
  geom_errorbar(aes(xmin = trend_q_0.05_2023,
                    xmax = trend_q_0.95_2023,
                    colour = non_overlap),
                alpha = 0.2,
                width = 0)+
  geom_abline(intercept = 0,slope = 1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  ggrepel::geom_text_repel(data = sp_labs,aes(label = Species_ID_Espèce_core),
                           min.segment.length = 0)+
  coord_cartesian(xlim = c(-8,8),
                  ylim = c(-8,8))+
  theme_bw()+
  scale_colour_viridis_d(end = 0.8, direction = -1)
print(comp_xy)
dev.off()




comp_xy_cov <- ggplot(data = trends_comp,
                  aes(x = reliab_cov_2023,
                      y = reliab.cov,
                      colour = coverage))+
  geom_point(alpha = 0.2)+
  geom_abline(intercept = 0,slope = 1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  facet_grid(cols = vars(trend_time),
             rows = vars(region),
             scales = "free")+
  theme_bw()+
  scale_colour_viridis_d(begin = 0.1, end = 0.9)

comp_xy_cov

pdf("coverage_comparison.pdf",
    width = 11, height = 8.5)
print(comp_xy_cov)
dev.off()


trends_comp_sel <- trends_comp %>%
  filter(diff_trend > 1 | diff_trend < -1,
         region_type == "survey-wide")


comp_xy_sel <- ggplot(data = trends_comp_sel,
                  aes(x = trend_2023,
                      y = trend,
                      colour = factor(bbs_num)))+
  geom_point()+
  geom_errorbar(aes(ymin = trend_q_0.05, ymax = trend_q_0.95),
                alpha = 0.4)+
  geom_errorbar(aes(xmin = trend_q_0.05_2023, xmax = trend_q_0.95_2023),
                 alpha = 0.4,
                orientation = "y")+
  geom_abline(intercept = 0,slope = 1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = "none")+
  geom_text_repel(aes(label = factor(bbs_num)))+
  facet_grid(rows = vars(trend_time),
             cols = vars(region),
             scales = "free")

comp_xy_sel



# Compile and organize trends and indices for SOCB ------------------------





states <- rnaturalearth::ne_states(country = c("Canada","United States of America")) %>%
  sf::st_drop_geometry() %>% select(postal,name_fr) %>%
  distinct() %>%
  rename(statprov_name_fr = name_fr)


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


bcrs <- sf::read_sf("data/bcr_2025_lakes12.gpkg") %>%
  mutate(region_type = "bcr",
         region = paste0("BCR_",bcr_label),
         geo.area = region,
         region_name_en = bcr_name_en,
         region_name_fr = bcr_name_fr) %>%
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
                       region_name_fr = "Zone complète de l'enquête")

regions <- bind_rows(survey_w,
                     countries,
                     prov_state,
                     bcrs,
                     bcr_country,
                     strats) |>
  select(-geo.area) |>
  rename(region_en = region_name_en,
         region_fr = region_name_fr)




## translate columns with region, region_type, trend_time, precision, coverage backcast_reliability, reliability,
## join tables for CWS website should provide translations
## add these new french columns to the header sheets
## DROP for_web
# reconcile with template and Catherine's email


trends <- trends %>%
  mutate(prob_LD = prob_decrease_50_percent,
         prob_MD = prob_decrease_25_percent - prob_decrease_50_percent,
         prob_LC = (prob_decrease_0_percent-prob_decrease_25_percent)+(prob_increase_0_percent-prob_increase_33_percent) ,
         prob_MI = prob_increase_33_percent - prob_increase_100_percent,
         prob_LI = prob_increase_100_percent,
         region_type = factor(region_type,
                              levels = c("survey-wide","country","prov_state","bcr","bcr_by_country","stratum"),
                              ordered = TRUE),
         region = ifelse(region_type == "bcr",
                         paste0("BCR_",region),
                         region))

test_probs <- trends %>%
  mutate(prob_test = prob_LD+prob_MD+prob_LC+prob_MI+prob_LI)

if(any(round(test_probs$prob_test,2) != 1)){warning("probabilites of change categories don't sum properly")}




cat_translate <- function(x){
  y <- as.character(x)
  r1 <- which(y == "Medium")
  y[r1] <- rep("Moyenne",length = length(r1))
  r1 <- which(y == "High")
  y[r1] <- rep("Élevée",length = length(r1))
  r1 <- which(y == "Low")
  y[r1] <- rep("Faible",length = length(r1))

  r1 <- which(y == "Long-term")
  y[r1] <- rep("a long term",length = length(r1))
  r1 <- which(y == "Short-term")
  y[r1] <- rep("a courte term",length = length(r1))
  r1 <- which(y == "Three-generation")
  y[r1] <- rep("trois générations",length = length(r1))

  r1 <- which(y == "bcr")
  y[r1] <- rep("rco",length = length(r1))
  r1 <- which(y == "bcr_by_country")
  y[r1] <- rep("rco_par_pays",length = length(r1))
  r1 <- which(y == "survey-wide")
  y[r1] <- rep("zone complète de l'enquête",length = length(r1))
  r1 <- which(y == "country")
  y[r1] <- rep("pays",length = length(r1))
  r1 <- which(y == "prov_state")
  y[r1] <- rep("prov_etat",length = length(r1))
  r1 <- which(y == "stratum")
  y[r1] <- rep("strate",length = length(r1))

  r1 <- which(y == "smooth")
  y[r1] <- rep("lisse",length = length(r1))

  r1 <- which(y == "full")
  y[r1] <- rep("complet",length = length(r1))

  return(y)
}

trends_bil <- trends %>%
  left_join(.,core_link,by = c("bbs_num" = "aou")) %>%
  left_join(regions, by = c("region","region_type")) |>
  mutate(type_region = cat_translate(region_type),
         period = cat_translate(trend_time),
         couverture = cat_translate(coverage),
         precision_fr = cat_translate(precision),
         poids_des_donnees_locales = cat_translate(backcast_reliab),
         fiabilite = cat_translate(reliability),
         strat_incl = paste0(strata_included,"; ",strata_excluded)) |>
  relocate(region_en,region_fr,region_type,type_region,species,espece,trend_time,period,start_year,end_year,
           starts_with("trend"),
           starts_with("percent"),
           width_of_95_percent_credible_interval,
           starts_with("prob_"),
           rel_abundance, n_routes, mean_n_routes, n_strata_included, backcast_flag,
           precision,precision_fr,coverage,couverture,backcast_reliab,poids_des_donnees_locales,
           reliability, fiabilite) %>%
  arrange(naturecounts_sort_order,region_type,region,start_year) |>
  rename(from_de = start_year,
         to_a = end_year,
         prob_LD_DI = prob_LD,
         prob_MD_DM = prob_MD,
         prob_LC_PdC = prob_LC,
         prob_MI_AM = prob_MI,
         prob_LI_AI = prob_LI,
         rel_ab = rel_abundance,
         rel_ab_obs = obs_rel_abundance,
         n_site = n_routes,
         mean_n_site_moyenne = mean_n_routes,
         n_strat_incl = n_strata_included,
         sequence_naturecounts = naturecounts_sort_order,
         width_CI_largeur_IC = width_of_95_percent_credible_interval,
         local_data_donnees_locales = backcast_flag) |>
  rename_with(.fn = ~gsub("trend","trend_tendence",.x),
              .cols = starts_with("trend")) |>
  rename_with(.fn = ~gsub("percent_change","percent_ch_pourcentage",.x),
              .cols = starts_with("percent_change")) |>
  rename_with(.fn = ~gsub("prob_decrease","prob_decrease_diminution",.x),
              .cols = starts_with("prob_decrease"))|>
  rename_with(.fn = ~gsub("prob_increase","prob_increase_augmentation",.x),
              .cols = starts_with("prob_increase")) |>
  rename(trend_time = trend_tendence_time) |>
  select(-c(strata_included,strata_excluded,region,for_web,
            bbs_num))




# Indices reorder ---------------------------------------------------------


indices_smooth <- indices_smooth %>%
  left_join(.,core_link,by = c("bbs_num" = "aou"))|>
  left_join(regions, by = c("region","region_type")) %>%
  mutate(region_type = factor(region_type,
                              levels = c("survey-wide","country","prov_state","bcr","bcr_by_country","stratum"),
                              ordered = TRUE),
         trend_time = factor(trend_time,
                             levels = c("Long-term","Three-generation","Short-term"),
                             ordered = TRUE)) %>%
  relocate(region,region_type,year,species,espece,trend_time,indices_type,
           starts_with("index"),
           starts_with("n_"),
           obs_mean, backcast_flag)%>%
    arrange(naturecounts_sort_order,region_type,region,trend_time,year)

indices_smooth_bil <- indices_smooth  |>
  mutate(type_region = cat_translate(region_type),
         period = cat_translate(trend_time),
         type_ind = cat_translate(indices_type),
         strat_incl = paste0(strata_included,"; ",strata_excluded)) |>
  relocate(region_en,region_fr,region_type,type_region,species,espece,trend_time,period,year,
           starts_with("index"),
           obs_mean, n_routes, n_routes_total, n_non_zero, backcast_flag) %>%
  arrange(naturecounts_sort_order,region_type,region,year) |>
  rename(year_an = year,
         rel_ab_obs = obs_mean,
         ind_type = indices_type,
         n_site = n_routes,
         n_site_total = n_routes_total,
         n_non_zero_pas_zero = n_non_zero,
         local_data_donnees_locales = backcast_flag) |>
  rename_with(.fn = ~gsub("index","ind",.x),
              .cols = starts_with("index"))|>
  select(-c(strata_included,strata_excluded,region,for_web,
            bbs_num))






indices <- indices %>%
  left_join(.,core_link,by = c("bbs_num" = "aou")) %>%
  left_join(regions, by = c("region","region_type")) %>%
  mutate(region_type = factor(region_type,
                              levels = c("survey-wide","country","prov_state","bcr","bcr_by_country","stratum"),
                              ordered = TRUE),
         trend_time = factor(trend_time,
                             levels = c("Long-term","Three-generation","Short-term"),
                             ordered = TRUE)) %>%
  relocate(region,region_type,year,species,espece,trend_time,indices_type,
           starts_with("index"),
           starts_with("n_"),
           obs_mean, backcast_flag)%>%
  arrange(naturecounts_sort_order,region_type,region,trend_time,year)

indices_bil <- indices|>
  mutate(type_region = cat_translate(region_type),
         period = cat_translate(trend_time),
         type_ind = cat_translate(indices_type),
         strat_incl = paste0(strata_included,"; ",strata_excluded)) |>
  relocate(region_en,region_fr,region_type,type_region,species,espece,trend_time,period,year,
           starts_with("index"),
           obs_mean, n_routes, n_routes_total, n_non_zero, backcast_flag) %>%
  arrange(naturecounts_sort_order,region_type,region,year) |>
  rename(year_an = year,
         rel_ab_obs = obs_mean,
         ind_type = indices_type,
         n_site = n_routes,
         n_site_total = n_routes_total,
         n_non_zero_pas_zero = n_non_zero,
         local_data_donnees_locales = backcast_flag) |>
  rename_with(.fn = ~gsub("index","ind",.x),
              .cols = starts_with("index"))|>
  select(-c(strata_included,strata_excluded,region,for_web,
            bbs_num))

# csv files with trends and indices for Google Drive ----------------------

saveRDS(indices,paste0(external_dir,"/Website/All_BBS_Full_Indices_",YYYY,".rds"))
saveRDS(indices_smooth,paste0(external_dir,"/Website/All_BBS_Smoothed_Indices_",YYYY,".rds"))
saveRDS(trends,paste0(external_dir,"/Website/All_BBS_Trends_",YYYY,".rds"))



saveRDS(indices_bil,paste0(external_dir,"/Website/All_BBS_Full_Indices_bil_",YYYY,".rds"))
saveRDS(indices_smooth_bil,paste0(external_dir,"/Website/All_BBS_Smoothed_Indices_bil_",YYYY,".rds"))
saveRDS(trends_bil,paste0(external_dir,"/Website/All_BBS_Trends_bil_",YYYY,".rds"))

# trends <- readRDS(paste0(external_dir,"/Website/All_BBS_Trends_",YYYY,".rds"))
# indices_smooth <- readRDS(paste0(external_dir,"/Website/All_BBS_Smoothed_Indices_",YYYY,".rds"))
# indices <- readRDS(paste0(external_dir,"/Website/All_BBS_Full_Indices_",YYYY,".rds"))
#
# tst <- indices_smooth %>% filter(region == "survey-wide", trend_time == "Long-term") %>% group_by(species) %>% summarise(n_test = n())

readr::write_excel_csv(indices_bil,paste0(external_dir,"/Website/All_Toutes_BBS_Full_Indices_Complets_",YYYY,".csv"))
readr::write_excel_csv(indices_smooth_bil,paste0(external_dir,"/Website/All_Toutes_BBS_Smoothed_Indices_Lisses_",YYYY,".csv"))
readr::write_excel_csv(trends_bil,paste0(external_dir,"/Website/All_Toutes_Tendance_BBS_Trends_",YYYY,".csv"))




inds_select <- indices_bil %>%
  filter(region_type %in% c("survey-wide","country"))
readr::write_excel_csv(inds_select,paste0(external_dir,"/Website/BBS_Full_Indices_Complets_survey-wide_zone_complete_country_pays_",YYYY,".csv"))

inds_select <- indices_bil %>%
  filter(region_type %in% c("prov_state"))
readr::write_excel_csv(inds_select,paste0(external_dir,"/Website/BBS_Full_Indices_Complets_prov_state_etat_",YYYY,".csv"))


inds_select <- indices_bil %>%
  filter(region_type %in% c("bcr"))
readr::write_excel_csv(inds_select,paste0(external_dir,"/Website/BBS_Full_Indices_Complets_bcr_rco_",YYYY,".csv"))

inds_select <- indices_bil %>%
  filter(region_type %in% c("bcr_by_country"))
readr::write_excel_csv(inds_select,paste0(external_dir,"/Website/BBS_Full_Indices_Complets_bcr_by_country_rco_par_pays_",YYYY,".csv"))

inds_select <- indices_bil %>%
  filter(region_type %in% c("stratum"))
readr::write_excel_csv(inds_select,paste0(external_dir,"/Website/BBS_Full_Indices_Complets_strata_strates_",YYYY,".csv"))


inds_select <- indices_smooth_bil %>%
  filter(region_type %in% c("survey-wide","country"))
readr::write_excel_csv(inds_select,paste0(external_dir,"/Website/BBS_Smoothed_Indices_Lisses_survey-wide_zone_complete_country_pays_",YYYY,".csv"))

inds_select <- indices_smooth_bil %>%
  filter(region_type %in% c("prov_state"))
readr::write_excel_csv(inds_select,paste0(external_dir,"/Website/BBS_Smoothed_Indices_Lisses_prov_state_etat_",YYYY,".csv"))


inds_select <- indices_smooth_bil %>%
  filter(region_type %in% c("bcr"))
readr::write_excel_csv(inds_select,paste0(external_dir,"/Website/BBS_Smoothed_Indices_Lisses_bcr_rco_",YYYY,".csv"))

inds_select <- indices_smooth_bil %>%
  filter(region_type %in% c("bcr_by_country"))
readr::write_excel_csv(inds_select,paste0(external_dir,"/Website/BBS_Smoothed_Indices_Lisses_bcr_by_country_rco_par_pays_",YYYY,".csv"))


inds_select <- indices_smooth_bil %>%
  filter(region_type %in% c("stratum"))
readr::write_excel_csv(inds_select,paste0(external_dir,"/Website/BBS_Smoothed_Indices_Lisses_strata_strates_",YYYY,".csv"))



#
# trends_select <- trends_bil %>%
#   filter(region_type %in% c("survey-wide","country","prov_state"))
# readr::write_excel_csv(trends_select,paste0(external_dir,"/Website/BBS_Trends_survey-wide_country_prov_state_",YYYY,".csv"))
# trends_select <- trends_bil %>%
#   filter(region_type %in% c("bcr","bcr_by_country"))
# readr::write_excel_csv(trends_select,paste0(external_dir,"/Website/BBS_Trends_bcr_bcr_by_country_",YYYY,".csv"))
# trends_select <- trends_bil %>%
#   filter(region_type %in% c("stratum"))
# readr::write_excel_csv(trends_select,paste0(external_dir,"/Website/BBS_Trends_strata_",YYYY,".csv"))






sp_no_coverage <- trends %>%
  mutate(w_cov = ifelse(is.na(reliab.cov),FALSE,TRUE)) %>%
  group_by(species,trend_time,w_cov) %>%
  summarise(n_regions = n(),
            .groups = "drop") %>%
  pivot_wider(.,values_from = n_regions,names_from = w_cov,names_prefix = "cov") %>%
  mutate(covTRUE = ifelse(is.na(covTRUE),0,covTRUE),
         covFALSE = ifelse(is.na(covFALSE),0,covFALSE),
         p_missing_coverage = covFALSE/covTRUE)


saveRDS(sp_no_coverage,"species_coverage_summary.rds")













# SOCB upload files -------------------------------------------------------


socb_areas <- read_csv("data/SOCB_regions.csv") %>%
  filter(results_code == "BBS") %>%
  select(area_code,area_name)

#
# trends_out <- trends %>%
#   filter((for_web == TRUE | region %in% c("survey-wide","United States of America")))

trends_out2 <- trends  %>%
  mutate(years = paste(start_year,end_year,sep = "-"),
         results_code = "BBS",
         season = "breeding",
         version = YYYY,
         area_code = ifelse(region == "survey-wide","Survey-wide",region),
         area_code = gsub(area_code,pattern = "United States of America",
                          replacement = "USA"),
         area_code = ifelse(region_type == "bcr",paste0("BCR_",region),area_code),
         model_type = "GAMYE",
         index_type = "mean_predicted_count",
         sample_size_units = "number of routes",
         trend_time = ifelse(trend_time == "Three-generation","3Gen-Recent",trend_time)) %>%
  left_join(.,socb_areas, by = "area_code") %>%
  #select(-area_code) %>%
  rename(species_name = species,
         species_code = Species_ID_core,
         species_id = naturecounts_species_id,
         period = trend_time,
         year_start = start_year,
         year_end = end_year,
         trnd = trend,
         lower_ci = trend_q_0.025,
         upper_ci = trend_q_0.975,
         percent_change = percent_change,
         percent_change_low = percent_change_q_0.025,
         percent_change_high = percent_change_q_0.975,
         prob_decrease_0 = prob_decrease_0_percent,
         prob_decrease_25 = prob_decrease_25_percent,
         prob_decrease_30 = prob_decrease_30_percent,
         prob_decrease_50 = prob_decrease_50_percent,
         prob_increase_0 = prob_increase_0_percent,
         prob_increase_33 = prob_increase_33_percent,
         prob_increase_100 = prob_increase_100_percent,
         precision_num = width_of_95_percent_credible_interval,
         precision_cat = precision,
         coverage_num = reliab.cov,
         coverage_cat = coverage,
         sample_size_alt = mean_n_routes,
         sample_size = n_routes,
         prob_LD = prob_LD,
         prob_MD = prob_MD,
         prob_LC = prob_LC,
         prob_MI = prob_MI,
         prob_LI = prob_LI)




trends_socb <- trends_out2 %>%
  select(results_code,
           version,
           area_code,
           area_name,
           season,
           period,
           species_name,
           species_code,
           species_id,
           years,
           year_start,
           year_end,
           trnd,
           lower_ci,
           upper_ci,
           index_type,
           model_type,
           percent_change,
           percent_change_low,
           percent_change_high,
           prob_decrease_0,
           prob_decrease_25,
           prob_decrease_30,
           prob_decrease_50,
           prob_increase_0,
           prob_increase_33,
           prob_increase_100,
           precision_num,
           precision_cat,
           coverage_num,
           coverage_cat,
           sample_size,
           sample_size_units,
           prob_LD,
           prob_MD,
           prob_LC,
           prob_MI,
           prob_LI)

readr::write_excel_csv(trends_socb,
                       paste0(external_dir,"/website/BBS_",YYYY,"_trends_for_socb.csv"))

# tmp <- trends_socb %>%
#   filter(species_name == "Killdeer")
# readr::write_excel_csv(tmp,
#                        paste0("website/BBS_",YYYY,"_trends_for_Killdeer.csv"))
#


# SOCB extra trends -------------------------------------------------------
#
# trends_out <- trends %>%
#   filter((for_web == FALSE & !(region %in% c("survey-wide","United States of America"))))
#
# trends_out2 <- trends_out  %>%
#   mutate(years = paste(start_year,end_year,sep = "-"),
#          results_code = "BBS",
#          season = "breeding",
#          version = YYYY,
#          area_code = ifelse(region == "survey-wide","survey-wideal",region),
#          area_code = gsub(area_code,pattern = "United States of America",
#                           replacement = "USA"),
#          area_code = ifelse(region_type == "bcr",paste0("BCR_",region),area_code),
#          model_type = "GAMYE",
#          index_type = "mean_predicted_count",
#          sample_size_units = "number of routes",
#          trend_time = ifelse(trend_time == "Three-generation","3Gen-Recent",trend_time)) %>%
#   left_join(.,socb_areas, by = "area_code") %>%
#   #select(-area_code) %>%
#   rename(species_name = species,
#          species_code = Species_ID_core,
#          species_id = naturecounts_species_id,
#          period = trend_time,
#          year_start = start_year,
#          year_end = end_year,
#          trnd = trend,
#          lower_ci = trend_q_0.025,
#          upper_ci = trend_q_0.975,
#          percent_change = percent_change,
#          percent_change_low = percent_change_q_0.025,
#          percent_change_high = percent_change_q_0.975,
#          prob_decrease_0 = prob_decrease_0_percent,
#          prob_decrease_25 = prob_decrease_25_percent,
#          prob_decrease_30 = prob_decrease_30_percent,
#          prob_decrease_50 = prob_decrease_50_percent,
#          prob_increase_0 = prob_increase_0_percent,
#          prob_increase_33 = prob_increase_33_percent,
#          prob_increase_100 = prob_increase_100_percent,
#          precision_num = width_of_95_percent_credible_interval,
#          precision_cat = precision,
#          coverage_num = reliab.cov,
#          coverage_cat = coverage,
#          sample_size_alt = mean_n_routes,
#          sample_size = n_routes,
#          prob_LD = prob_LD,
#          prob_MD = prob_MD,
#          prob_LC = prob_LC,
#          prob_MI = prob_MI,
#          prob_LI = prob_LI)
#
#
#
#
# trends_socb <- trends_out2 %>%
#   select(results_code,
#          version,
#          area_code,
#          area_name,
#          season,
#          period,
#          species_name,
#          species_code,
#          species_id,
#          years,
#          year_start,
#          year_end,
#          trnd,
#          lower_ci,
#          upper_ci,
#          index_type,
#          model_type,
#          percent_change,
#          percent_change_low,
#          percent_change_high,
#          prob_decrease_0,
#          prob_decrease_25,
#          prob_decrease_30,
#          prob_decrease_50,
#          prob_increase_0,
#          prob_increase_33,
#          prob_increase_100,
#          precision_num,
#          precision_cat,
#          coverage_num,
#          coverage_cat,
#          sample_size,
#          sample_size_units,
#          prob_LD,
#          prob_MD,
#          prob_LC,
#          prob_MI,
#          prob_LI)
#
# readr::write_excel_csv(trends_socb,
#                        paste0(external_dir,"/website/BBS_",YYYY,"_extra_trends_for_socb.csv"))
#
#



# SOCB indices ------------------------------------------------------------


smooth_join <- indices_smooth %>%
  select(species,region,region_type,trend_time,
         year,index) %>%
  rename(smooth_index = index)

indices_socb <- indices %>%
  #filter((for_web == TRUE | region %in% c("survey-wide","United States of America"))) %>%
  inner_join(.,smooth_join,
             by = c("species",
                    "region",
                    "region_type",
                    "trend_time",
                    "year")) %>%
  group_by(species,region,region_type,trend_time) %>%
  mutate(LOESS_index = loess_func(index,year),
         results_code = "BBS",
         season = "breeding",
         version = YYYY,
         area_code = ifelse(region == "survey-wide","Survey-wide",region),
         area_code = gsub(area_code,pattern = "United States of America",
                          replacement = "USA"),
         area_code = ifelse(region_type == "bcr",paste0("BCR_",region),area_code),
         trend_time = as.character(trend_time),
         trend_time = ifelse(trend_time == "Three-generation","3Gen-Recent",trend_time)) %>%
  ungroup() %>%
  left_join(.,socb_areas, by = "area_code") %>%
  #select(-area_code) %>%
  rename(species_name = species,
         species_code = Species_ID_core,
         species_id = naturecounts_species_id,
         period = trend_time,
         upper_ci = index_q_0.95,
         lower_ci = index_q_0.05) %>%
  select(-c(index_q_0.025,
            index_q_0.975,
            region_type,
            region)) %>%
  relocate(results_code,
           version,
           area_code,
           season,
           period,
           species_name,
           species_code,
           species_id,
           year,
           index,
           upper_ci,
           lower_ci,
           LOESS_index,
           smooth_index) %>%
  select(results_code,
         version,
         area_code,
         species_id,
         area_name,
         season,
         period,
         species_name,
         species_code,
         year,
         index,
         upper_ci,
         lower_ci,
         LOESS_index,
         smooth_index)

readr::write_excel_csv(indices_socb,
                       file = paste0(external_dir,"/website/BBS_",YYYY,"_annual_indices_for_socb.csv"))


# tmp <- indices_socb %>%
#   filter(species_name == "Killdeer")
# readr::write_excel_csv(tmp,
#                        file = paste0("website/BBS_",YYYY,"_annual_indices_for_Killdeer.csv"))
#





# SOCB extra indices ------------------------------------------------------

#
# indices_socb <- indices %>%
#   filter((for_web == FALSE & !(region %in% c("survey-wide","United States of America")))) %>%
#   inner_join(.,smooth_join,
#              by = c("species",
#                     "region",
#                     "region_type",
#                     "trend_time",
#                     "year")) %>%
#   group_by(species,region,region_type,trend_time) %>%
#   mutate(LOESS_index = loess_func(index,year),
#          results_code = "BBS",
#          season = "breeding",
#          version = YYYY,
#          area_code = ifelse(region == "survey-wide","survey-wideal",region),
#          area_code = gsub(area_code,pattern = "United States of America",
#                           replacement = "USA"),
#          area_code = ifelse(region_type == "bcr",paste0("BCR_",region),area_code),
#          trend_time = as.character(trend_time),
#          trend_time = ifelse(trend_time == "Three-generation","3Gen-Recent",trend_time)) %>%
#   ungroup() %>%
#   left_join(.,socb_areas, by = "area_code") %>%
#   #select(-area_code) %>%
#   rename(species_name = species,
#          species_code = Species_ID_core,
#          species_id = naturecounts_species_id,
#          period = trend_time,
#          upper_ci = index_q_0.95,
#          lower_ci = index_q_0.05) %>%
#   select(-c(index_q_0.025,
#             index_q_0.975,
#             region_type,
#             region)) %>%
#   relocate(results_code,
#            version,
#            area_code,
#            season,
#            period,
#            species_name,
#            species_code,
#            species_id,
#            year,
#            index,
#            upper_ci,
#            lower_ci,
#            LOESS_index,
#            smooth_index) %>%
#   select(results_code,
#          version,
#          area_code,
#          area_name,
#          season,
#          period,
#          species_name,
#          species_code,
#          species_id,
#          year,
#          index,
#          upper_ci,
#          lower_ci,
#          LOESS_index,
#          smooth_index)
#
# readr::write_excel_csv(indices_socb,
#                        file = paste0(external_dir,"/website/BBS_",YYYY,"_extra_annual_indices_for_socb.csv"))
#
#
#
#



