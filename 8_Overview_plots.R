


###  - State of Canada's Birds
YYYY <- 2024

webmaps <- FALSE # set to true if needing to create all map images for ECCC website

library(bbsBayes2)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(foreach)
library(doParallel)
#setwd("C:/GitHub/CWS_2023_BBS_Analyses")
output_dir <- "D:/BBS_Trends_CWS/output"
#output_dir <- "output"
external_dir <- "D:/BBS_Trends_CWS"
# output_dir <- "F:/CWS_2023_BBS_Analyses/output"

sp_list <- readRDS("sp_list_w_generations.rds") %>%
  filter(model == TRUE)
#
# three_gens <- read_csv("data/full_bbs_species_list_w_generation_length.csv")
#
# three_gens <- three_gens %>%
#   select(aou,GenLength)
#
# sp_list <- sp_list %>%
#   inner_join(.,three_gens,
#              by = c("aou"))
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
#   mutate(species_id = ifelse(species_id2 != species_id,species_id2,species_id)) %>%
#   select(species_id,species_code) %>%
#   mutate(aou = as.integer(species_code),
#          naturecounts_species_id = species_id) %>%
#   distinct() %>%
#   select(-c(species_id,species_code))
#
# sp_list <- sp_list %>%
#   left_join(.,naturecounts_codes,
#             by = "aou") %>%
#   arrange(Sort_Order_core)
#
# # Compile all trends and indices ------------------------------------------------------


# Compare to last year's trends -------------------------------------------


lastyear = read_csv(past0("data/All_BBS_trends_",YYYY-1,".csv"))

lastyear_inds <- read_csv(paste0("data/All_BBS_Full_indices_",YYYY-1,".csv"))
lastyear_inds_smooth <- read_csv(paste0("data/All_BBS_Smoothed_Indices_",YYYY-1,".csv"))


ly_trends <- lastyear[,c("species","bbs_num","region","trend_time",
                         "n_strata_included","n_routes",
                         "trend",
                         "trend_q_0.05","trend_q_0.95",
                         "width_of_95_percent_credible_interval")] %>%
  filter(region %in% c("continent","Canada","United States of America")) %>%
  select(-c(species))



species_to_run <- sp_list %>%
  arrange(naturecounts_sort_order)

pdf(file = paste0("Figures/BBS_High_level_summary_",YYYY,".pdf"),
    height = 9,
    width = 17)

for(jj in (1:nrow(species_to_run))){


  species <- as.character(species_to_run[jj,"english"])
  espece <- as.character(species_to_run[jj,"french"])
  aou <- as.integer(species_to_run[jj,"aou"])

  # identifying first years for selected species ----------------------------
  fy <- 1970
  if(aou %in% c(4661,4660)){ #Alder and Willow Flycatcher
    fy <- 1978 #5 years after the split
  }
  if(aou %in% c(10,11,22860)){ # Clark's and Western Grebe and EUCD
    fy <- 1990 #5 years after the split and first year EUCD observed on > 3 BBS routes
  }
  if(aou == 6121){ # CAve Swallow
    fy = 1985
  }


  if(file.exists(paste0(external_dir,"/Figures/temp_rds_storage/",aou,"_maps.RDS"))){
  species_f_bil <- gsub(paste(species,espece),pattern = "[[:space:]]|[[:punct:]]",
                        replacement = "_")

  trends_ly <- ly_trends %>%
    filter(bbs_num == aou) %>%
    mutate(version = "Last year")

  if(nrow(trends_ly) > 0){
  trends_1 <- readRDS(paste0(external_dir,"/Trends/",aou,"_trends.rds")) %>%
    filter(region %in% c("continent","Canada","United States of America")) %>%
    mutate(aou = aou,
           version = "This year") %>%
    bind_rows(.,trends_ly)
  }else{
    trends_1 <- readRDS(paste0(external_dir,"/Trends/",aou,"_trends.rds")) %>%
      filter(region %in% c("continent","Canada","United States of America")) %>%
      mutate(aou = aou,
             version = "This year")
}
tplot <- ggplot(data = trends_1)+
  geom_pointrange(aes(x = region,y = trend,ymin = trend_q_0.05,ymax = trend_q_0.95,
                      colour = version),
                  position = position_dodge(width = 0.5))+
  facet_wrap(vars(trend_time),
             scales = "free_x")+
  scale_colour_viridis_d(direction = -1,
                         name = paste("Trends",YYYY,"\n  and",YYYY-1))+
  coord_flip()+
  geom_hline(yintercept = 0)+
  theme_bw()


  tmaps <- readRDS(paste0(external_dir,"/Figures/temp_rds_storage/",aou,"_maps.RDS"))

 # trajs <- readRDS(paste0(external_dir,"/Figures/temp_rds_storage/",aou,"_highlevel_simple_trajs.RDS"))

  inds <- readRDS(paste0(external_dir,"/Indices/Inds_",aou,".rds"))

  ind <- readRDS(paste0(external_dir,"/Indices/Ind_plot_",aou,".rds"))


  trajs_all <- plot_indices(ind,
                        title = FALSE,
                        ci_width = 0.9,
                        add_observed_means = FALSE,
                        add_number_routes = FALSE,
                        min_year = fy,
                        title_size = 12,
                        axis_title_size = 10,
                        axis_text_size = 10)


  lastyear_inds_sp <- lastyear_inds %>%
    filter(bbs_num == aou,
           trend_time == "Long-term")

  lastyear_inds_smooth_sp <- lastyear_inds_smooth %>%
    filter(bbs_num == aou,
           trend_time == "Long-term")

  trajs <- vector("list",3)
  names(trajs) <- c("continent","Canada","United_States_of_America")

  for(j in names(trajs_all)[1:3]){
    if(!j %in% c("continent","Canada","United_States_of_America")){next}

    t1 <- trajs_all[[j]]
    #t2 <- trajshort[[j]]
    if(grepl("United_States_of_America",j)){
      rr <- str_replace_all(j,"_"," ")
    }else{
      rr <- str_replace_all(j,"_","-")
    }

    labl <- paste(esp,"/",sp,"-",rr)

    n1 <- inds$indices %>%
      filter(region == rr,
             year >= fy)

    strats <- str_split_1(unlist(n1[1,"strata_included"]),
                          pattern = " ; ")
    upy <- max(n1$index_q_0.95,na.rm = TRUE)/2


    t1plot <- t1 +
      geom_ribbon(data = n1, aes(x = year,y = index,ymin = index_q_0.05,ymax = index_q_0.95),
                  fill = grey(0.5),alpha = 0.2)+
      geom_line(data = n1, aes(x = year,y = index),
                colour = grey(0.5))+
      labs(subtitle = labl)






      ly_inds <- lastyear_inds_sp %>%
        filter(region == rr)
      ly_inds_smooth <- lastyear_inds_smooth_sp %>%
        filter(region == rr)

      if(j == "continent"){
        ly_inds <- lastyear_inds_sp %>%
          filter(region == "continent")
      }

      if(nrow(ly_inds) > 0){
        tmpPlot <- t1plot +
          geom_line(data = ly_inds,
                    aes(x = year,
                        y = index_q_0.05),
                    colour = "darkgreen",alpha = 0.3, linetype = 6)+
          geom_line(data = ly_inds,
                    aes(x = year,
                        y = index_q_0.95),
                    colour = "darkgreen",alpha = 0.3, linetype = 6)+
          geom_line(data = ly_inds_smooth,
                    aes(x = year,
                        y = index),
                    colour = "darkgreen",alpha = 0.4, linetype = 6)

        trajs[[j]] <- tmpPlot
      }else{
        trajs[[j]] <- t1plot
      }


  }


  design <- "
135
246
"

  tt_sel <- trends_1 %>%
    filter(trend_time == "Long-term", version == "This year",
           region == "continent")
  main_title <- paste(species,espece,"continent long-term trend has",tt_sel$reliability,"overall reliability",
                      ": coverage",tt_sel$reliab.cov*100,"% :", tt_sel$precision, "precision :",
                      tt_sel$backcast_reliab,"local data ")

  if(is.null(trajs[["Canada"]])){
    tt_sel <- trends_1 %>%
      filter(trend_time == "Long-term", version == "This year",
             region %in% c("United States of America"))
    main_caption <- paste("US long-term trend has",tt_sel$reliability,"overall reliability",
                          ": coverage",tt_sel$reliab.cov*100,"% of the species' range :", tt_sel$precision, "precision :",
                          tt_sel$backcast_reliab,"local data  :", tt_sel$n_routes,"routes total and ",
                          tt_sel$mean_n_routes,"on average, each year: ","local data available for",tt_sel$backcast_flag*100,
                          "% of the years and regions included")

    layt <- trajs[[1]] + trajs[[3]] + plot_spacer() +
      tplot + tmaps[[1]] + tmaps[[2]] +
      plot_layout(design = design,
                  guides = "collect",
                  widths = 1)+
      plot_annotation(title = main_title,
                      caption = main_caption,
                      theme = theme(plot.caption = element_text(size = 10),
                                    plot.title = element_text(size = 12)))


  }else{

    tt_sel <- trends_1 %>%
      filter(trend_time == "Long-term", version == "This year",
             region %in% c("Canada"))
    can_title <- paste("Canada long-term trend has",tt_sel$reliability,"overall reliability",
                       ": coverage",tt_sel$reliab.cov*100,"% of the species' range :", tt_sel$precision, "precision :",
                       tt_sel$backcast_reliab,"local data :", tt_sel$n_routes,"routes total and",
                       tt_sel$mean_n_routes,"on average, each year: ","local data available for",tt_sel$backcast_flag*100,
                       "% of the years and regions included")

    tt_sel <- trends_1 %>%
      filter(trend_time == "Long-term", version == "This year",
             region %in% c("United States of America"))
    us_title <- paste("US long-term trend has",tt_sel$reliability,"overall reliability",
                      ": coverage",tt_sel$reliab.cov*100,"% of the species' range :", tt_sel$precision, "precision :",
                      tt_sel$backcast_reliab,"local data  :", tt_sel$n_routes,"routes total and ",
                      tt_sel$mean_n_routes,"on average, each year: ","local data available for",tt_sel$backcast_flag*100,
                      "% of the years and regions included")
    main_caption <- paste0(can_title,"\n",us_title)

    layt <- trajs[[1]] + trajs[[3]] + trajs[[2]] +
      tplot + tmaps[[1]] + tmaps[[2]] +
      plot_layout(design = design,
                  guides = "collect",
                  widths = 1) +
      plot_annotation(title = main_title,
                      caption = main_caption,
                      theme = theme(plot.caption = element_text(size = 10),
                                    plot.title = element_text(size = 12)))
}
  print(layt)

  print(species)
  }else{
  print(paste("No estimates for",species,aou))
}
}

dev.off()








# Plotting trend maps -----------------------------------------------------

re_run <- TRUE

start_years <- c("Long-term","Short-term","Three-generation")



n_cores <- 6
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)


test <- foreach(jj = rev(c(1:nrow(species_to_run))),
                .packages = c("bbsBayes2",
                              "tidyverse",
                              "cmdstanr",
                              "patchwork"),
                .errorhandling = "pass") %dopar%
  {

  species <- as.character(species_to_run[jj,"english"])
  espece <- as.character(species_to_run[jj,"french"])
  aou <- as.integer(species_to_run[jj,"aou"])
  species_f_bil <- gsub(paste(species,espece),pattern = "[[:space:]]|[[:punct:]]",
                        replacement = "_")

  if(!file.exists(paste0(external_dir,"/Indices/Inds_",aou,".rds")) &
     (!file.exists(paste0(external_dir,"/Figures/trend_maps/",species_f_bil,"_trend_maps.pdf")) |
      re_run)){

    tmaps <- readRDS(paste0(external_dir,"/Figures/temp_rds_storage/",aou,"_maps.RDS"))
    qmaps <- readRDS(paste0(external_dir,"/Figures/temp_rds_storage/",aou,"_quart_maps.RDS"))

    pdf(paste0(external_dir,"/Figures/trend_maps/",species_f_bil,"_trend_maps.pdf"),
        width = 11,
        height = 8.5)


    for(j in (start_years)){
      tt <- tmaps[[j]]+
        labs(title = paste(j,":",espece,"/",species))

    print(tt)

    tt <- qmaps[[j]]+
      plot_annotation(title = paste(j,":",espece,"/",species))

    print(tt)


      }
  dev.off()


    }


  }
parallel::stopCluster(cluster)

