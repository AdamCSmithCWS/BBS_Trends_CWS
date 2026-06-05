


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


lastyear = read_csv(paste0("data/All_BBS_trends_",YYYY-1,".csv"))

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





use_last_year <- FALSE

if(use_last_year){
pdf(file = paste0("Figures/BBS_High_level_summary_",YYYY,".pdf"),
    height = 9,
    width = 17)
}else{
  pdf(file = paste0("Figures/BBS_trends_summary_",YYYY,".pdf"),
      height = 9,
      width = 17)
}
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

  if(nrow(trends_ly) > 0 & use_last_year){
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

  trends_1 <- trends_1 |>
    mutate(region = ifelse(region == "continent",
                           "Survey-wide",
                           region))

tplot <- ggplot(data = trends_1)+
  geom_pointrange(aes(x = region,y = trend,ymin = trend_q_0.05,ymax = trend_q_0.95,
                      colour = version),
                  position = position_dodge(width = 0.5))+
  facet_wrap(vars(trend_time),
             scales = "free_x")

if(use_last_year){
  tplot <- tplot +
    coord_flip()+
    geom_hline(yintercept = 0)+
    theme_bw()+
    ylab("Trends %/year")+
    xlab("")+
    scale_colour_viridis_d(direction = -1,
                           name = paste("Trends",YYYY,"\n  and",YYYY-1))
}else{
  tplot <- tplot +
    coord_flip()+
    geom_hline(yintercept = 0)+
    ylab("Trends %/year")+
    xlab("")+
    theme_bw()+
    theme(legend.position = "none")
}

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

if(use_last_year){
  lastyear_inds_sp <- lastyear_inds %>%
    filter(bbs_num == aou,
           trend_time == "Long-term")

  lastyear_inds_smooth_sp <- lastyear_inds_smooth %>%
    filter(bbs_num == aou,
           trend_time == "Long-term")
}else{
  lastyear_inds_sp <- NULL
  lastyear_inds_smooth_sp <- NULL

}
  trajs <- vector("list",3)
  names(trajs) <- c("continent","Canada","United_States_of_America")

  for(j in names(trajs_all)[1:3]){
    if(!j %in% c("continent","Canada","United_States_of_America")){next}

    t1 <- trajs_all[[j]]
    #t2 <- trajshort[[j]]
    if(grepl("United_States_of_America",j)){
      rr <- str_replace_all(j,"_"," ")
    }
  if(j == "continent"){
    rr <- "Survey-wide"
  }
    if(j == "Canada"){
      rr <- j
    }
    labl <- paste(species,"-",rr)

    n1 <- inds$indices |>
      mutate(region = ifelse(region == "continent",
                             "Survey-wide",
                             region)) %>%
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





if(use_last_year){
      ly_inds <- lastyear_inds_sp |>
        mutate(region = ifelse(region == "continent",
                               "Survey-wide",
                               region)) %>%
        filter(region == rr)
      ly_inds_smooth <- lastyear_inds_smooth_sp|>
        mutate(region = ifelse(region == "continent",
                               "Survey-wide",
                               region)) %>%
        filter(region == rr)


      if(j == "continent"){
        ly_inds <- lastyear_inds_sp %>%
          filter(region == "continent")|>
          mutate(region = ifelse(region == "continent",
                                 "Survey-wide",
                                 region))
      }
}else{
  ly_inds <- NULL
  ly_inds_smooth <- NULL
}
      if(use_last_year){
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
      }else{
        trajs[[j]] <- t1plot
      }


} # end regional loop

  design <- "
135
246
"

  tt_sel <- trends_1 %>%
    filter(trend_time == "Long-term", version == "This year",
           region == "Survey-wide")
  main_title <- paste(species,"Survey-wide long-term trend has",tt_sel$reliability,"overall reliability",
                      "; coverage",tt_sel$reliab.cov*100,"% ; ", tt_sel$precision, "precision ; ",
                      tt_sel$backcast_reliab,"local data ")

  if(is.null(trajs[["Canada"]])){
    tt_sel <- trends_1 %>%
      filter(trend_time == "Long-term", version == "This year",
             region %in% c("United States of America"))
    main_caption <- paste("US long-term trend has",tt_sel$reliability,"overall reliability",
                          "; coverage",tt_sel$reliab.cov*100,"% of the species' range ;", tt_sel$precision, "precision ; ",
                          tt_sel$backcast_reliab,"local data  ; ", tt_sel$n_routes,"routes total and ",
                          tt_sel$mean_n_routes,"on average, each year; ","local data available for",tt_sel$backcast_flag*100,
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
                       "; coverage",tt_sel$reliab.cov*100,"% of the species' range ; ", tt_sel$precision, "precision ; ",
                       tt_sel$backcast_reliab,"local data ; ", tt_sel$n_routes,"routes total and",
                       tt_sel$mean_n_routes,"on average, each year; ","local data available for",tt_sel$backcast_flag*100,
                       "% of the years and regions included")

    tt_sel <- trends_1 %>%
      filter(trend_time == "Long-term", version == "This year",
             region %in% c("United States of America"))
    us_title <- paste("US long-term trend has",tt_sel$reliability,"overall reliability",
                      "; coverage",tt_sel$reliab.cov*100,"% of the species' range ; ", tt_sel$precision, "precision ; ",
                      tt_sel$backcast_reliab,"local data  ; ", tt_sel$n_routes,"routes total and ",
                      tt_sel$mean_n_routes,"on average, each year; ","local data available for",tt_sel$backcast_flag*100,
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

  print(paste(species,round(jj/nrow(species_to_run),2)))
  }else{
  print(paste("No estimates for",species,aou))
}
}

dev.off()







# French overview file ----------------------------------------------------




if(use_last_year){
  pdf(file = paste0("Figures/BBS_Haut_resume_",YYYY,".pdf"),
      height = 9,
      width = 17)
}else{
  pdf(file = paste0("Figures/BBS_tendence_resume_",YYYY,".pdf"),
      height = 9,
      width = 17)
}
for(jj in c(1:nrow(species_to_run))){


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

    if(nrow(trends_ly) > 0 & use_last_year){
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

    trends_1 <- trends_1 |>
      mutate(trend_time = ifelse(trend_time == "Three-generation",
                                 "Trois générations",
                                 trend_time),
             trend_time = ifelse(trend_time == "Long-term",
                                 "Long terme",
                                 trend_time),
             trend_time = ifelse(trend_time == "Short-term",
                                 "Court terme",
                                 trend_time),
             region = ifelse(region == "United States of America",
                             "États-Unis d'Amérique",
                             region),
             region = ifelse(region == "continent",
                             "Zone complète de l'enquête",
                             region))|>
      mutate(reliability = cat_translate(reliability),
             precision = cat_translate(precision),
             backcast_reliab = cat_translate(backcast_reliab))



    tplot <- ggplot(data = trends_1)+
      geom_pointrange(aes(x = region,y = trend,ymin = trend_q_0.05,ymax = trend_q_0.95,
                          colour = version),
                      position = position_dodge(width = 0.5))+
      facet_wrap(vars(trend_time),
                 scales = "free_x")

    if(use_last_year){
      tplot <- tplot +
        coord_flip()+
        geom_hline(yintercept = 0)+
        theme_bw()+
        ylab("Tendances %/an")+
        xlab("")+
        scale_colour_viridis_d(direction = -1,
                               name = paste("Tendances",YYYY,"\n  et",YYYY-1))
    }else{
      tplot <- tplot +
        coord_flip()+
        geom_hline(yintercept = 0)+
        ylab("Tendances %/an")+
        xlab("")+
        theme_bw()+
        theme(legend.position = "none")
    }

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

    if(use_last_year){
      lastyear_inds_sp <- lastyear_inds %>%
        filter(bbs_num == aou,
               trend_time == "Long-term")

      lastyear_inds_smooth_sp <- lastyear_inds_smooth %>%
        filter(bbs_num == aou,
               trend_time == "Long-term")
    }else{
      lastyear_inds_sp <- NULL
      lastyear_inds_smooth_sp <- NULL

    }
    trajs <- vector("list",3)
    names(trajs) <- c("continent","Canada","United_States_of_America")

    for(j in names(trajs_all)[1:3]){
      if(!j %in% c("continent","Canada","United_States_of_America")){next}

      rr2 <- ifelse(grepl("United_States_of_America",j),"United States of America",j)

      rr <- ifelse(grepl("United_States_of_America",j),
                   "États-Unis d'Amérique",
                  j)

      rr <- ifelse(j == "continent",
                   "Zone complète de l'enquête",
                   rr)


      t1 <- trajs_all[[j]]



      labl <- paste(espece,"-",rr)

      n1 <- inds$indices %>%
        filter(region == rr2,
               year >= fy)

      strats <- str_split_1(unlist(n1[1,"strata_included"]),
                            pattern = " ; ")
      upy <- max(n1$index_q_0.95,na.rm = TRUE)/2


      t1plot <- t1 +
        geom_ribbon(data = n1, aes(x = year,y = index,ymin = index_q_0.05,ymax = index_q_0.95),
                    fill = grey(0.5),alpha = 0.2)+
        geom_line(data = n1, aes(x = year,y = index),
                  colour = grey(0.5))+
        labs(subtitle = labl)+
        xlab("Année")+
        ylab("Indice annuel d'abondance\n(nombre moyen)")





      if(use_last_year){
        ly_inds <- lastyear_inds_sp %>%
          filter(region == rr2)
        ly_inds_smooth <- lastyear_inds_smooth_sp %>%
          filter(region == rr2)


        if(j == "continent"){
          ly_inds <- lastyear_inds_sp %>%
            filter(region == "continent")
        }
      }else{
        ly_inds <- NULL
        ly_inds_smooth <- NULL
      }
      if(use_last_year){
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
      }else{
        trajs[[j]] <- t1plot
      }


    } # end regional loop

    design <- "
135
246
"

    tmaps_alt1 <- tmaps[[1]]+
      labs(title = "Tendance à Long terme")
    tmaps_alt2 <- tmaps[[2]]+
      labs(title = "Tendance à Court terme")

    tt_sel <- trends_1 %>%
      filter(trend_time == "Long terme", version == "This year",
             region == "Zone complète de l'enquête")
    main_title <- paste("La tendance à long terme de Zone complète de l'enquête",espece,"présente une fiabilité globale",tt_sel$reliability,
                        "; couverture de ",tt_sel$reliab.cov*100,"% ; precision ", tt_sel$precision,"; poids des données locales ",
                        tt_sel$backcast_reliab)

    if(is.null(trajs[["Canada"]])){
      tt_sel <- trends_1 %>%
        filter(trend_time == "Long terme", version == "This year",
               region %in% c("États-Unis d'Amérique"))
      main_caption <- paste("La tendance à long terme de États-Unis d'Amérique présente une fiabilité globale",tt_sel$reliability,
                            "; couverture de ",tt_sel$reliab.cov*100,"% ; precision ", tt_sel$precision,"; poids des données locales ",
                            tt_sel$backcast_reliab,"; ",tt_sel$n_routes,"parcours au total et ",
                            tt_sel$mean_n_routes,"en moyenne par an; ","et données locales de",tt_sel$backcast_flag*100,
                            "% des regions et ans.")

      layt <- trajs[[1]] + trajs[[3]] + plot_spacer() +
        tplot + tmaps_alt1 + tmaps_alt2 +
        plot_layout(design = design,
                    guides = "collect",
                    widths = 1)+
        plot_annotation(title = main_title,
                        caption = main_caption,
                        theme = theme(plot.caption = element_text(size = 10),
                                      plot.title = element_text(size = 12)))


    }else{

      tt_sel <- trends_1 %>%
        filter(trend_time == "Long terme", version == "This year",
               region %in% c("Canada"))
      can_title <- paste("La tendance à long terme de Canada présente une fiabilité globale",tt_sel$reliability,
                         "; couverture de ",tt_sel$reliab.cov*100,"% ; precision", tt_sel$precision,"; poids des données locales ",
                         tt_sel$backcast_reliab,"; ",tt_sel$n_routes,"parcours au total et ",
                         tt_sel$mean_n_routes,"en moyenne par an; ","et données locales de",tt_sel$backcast_flag*100,
                         "% des regions et ans.")

      tt_sel <- trends_1 %>%
        filter(trend_time == "Long terme", version == "This year",
               region %in% c("États-Unis d'Amérique"))

      us_title <- paste("La tendance à long terme de États-Unis d'Amérique présente une fiabilité globale",tt_sel$reliability,
                        "; couverture de ",tt_sel$reliab.cov*100,"% ; precision", tt_sel$precision,"; poids des données locales ",
                        tt_sel$backcast_reliab,"; ",tt_sel$n_routes,"parcours au total et ",
                        tt_sel$mean_n_routes,"en moyenne par an; ","et données locales de",tt_sel$backcast_flag*100,
                        "% des regions et ans.")
      main_caption <- paste0(can_title,"\n",us_title)

      layt <- trajs[[1]] + trajs[[3]] + trajs[[2]] +
        tplot + tmaps_alt1 + tmaps_alt2 +
        plot_layout(design = design,
                    guides = "collect",
                    widths = 1) +
        plot_annotation(title = main_title,
                        caption = main_caption,
                        theme = theme(plot.caption = element_text(size = 10),
                                      plot.title = element_text(size = 12)))
    }
    print(layt)

    print(paste(species,round(jj/nrow(species_to_run),2)))
  }else{
    print(paste("No estimates for",species,aou))
  }
}

dev.off()








# Plotting trend maps -----------------------------------------------------

re_run <- TRUE

start_years <- c("Long-term","Short-term","Three-generation")



n_cores <- 8
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

  if(file.exists(paste0(external_dir,"/Indices/Inds_",aou,".rds")) &
     (!file.exists(paste0(external_dir,"/Figures/trend_maps/",species_f_bil,"_trend_maps.pdf")) |
      re_run)){

    tmaps <- readRDS(paste0(external_dir,"/Figures/temp_rds_storage/",aou,"_maps.RDS"))
    qmaps <- readRDS(paste0(external_dir,"/Figures/temp_rds_storage/",aou,"_quart_maps.RDS"))

    pdf(paste0(external_dir,"/Figures/trend_maps/",species_f_bil,"_trend_maps.pdf"),
        width = 11,
        height = 8.5)


    for(j in (start_years)){
      tt <- tmaps[[j]]+
        labs(title = paste(j,"; ",espece,"/",species))

    print(tt)

    tt <- qmaps[[j]]+
      plot_annotation(title = paste(j,"; ",espece,"/",species))

    print(tt)


      }
  dev.off()


    }


  }
parallel::stopCluster(cluster)

