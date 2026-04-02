## Calculates all possible three-generation trends for every species
## graphs the trend maps, and estimated trends through time to understand
## how the three-generation trends are changing through time.
## Intended to provide some assessment of the acceleration of trends
## slowing declines, increasing declines, etc.
## also the dependency of a trend estimate on the end years

library(bbsBayes2)
library(tidyverse)
library(foreach)
library(doParallel)
library(patchwork)

YYYY <- 2024
short_time <- 10

#setwd("C:/Users/SmithAC/Documents/GitHub/CWS_2023_BBS_Analyses")
#setwd("C:/GitHub/CWS_2023_BBS_Analyses")


output_dir <- "D:/BBS_Trends_CWS/output"
#output_dir <- "output"
external_dir <- "D:/BBS_Trends_CWS"
# output_dir <- "F:/CWS_2022_BBS_Analyses/output"


# custom functions to calculate reliability categories and determine website inclusion
source("functions/web_trends.R")
source("functions/reliability.R")



#output_dir <- "output"
n_cores = 8
re_run <- TRUE #set to TRUE to overwrite any previous output from this script


sp_list <- readRDS("sp_list_w_generations.rds") %>%
  filter(model == TRUE)

#
# sp_rerun <- c("Northern Shrike","Willow Ptarmigan", "Herring Gull",
#               "Common Loon",
#               "American Pipit",
#               "Redpoll (Common/Hoary)")
# sp_list <- sp_list %>%
#   filter(english %in% sp_rerun)


regs_to_estimate <- c("continent","country","prov_state","bcr","stratum","bcr_by_country")

# load previous trend data -----------------------------------------------------------

lastyear = read_csv(paste0("data/All_BBS_trends_",YYYY-1,".csv"))




# CV_threshold <- function(m,ci,thresh = 100){
#   y <- ifelse(ci/m > thresh,TRUE,FALSE)
#   return(y)
# }
#



# reliability category definitions ----------------------------------------

prec_cuts = c(abs(2*((0.7^(1/20))-1)),
              abs(2*((0.5^(1/20))-1)))*100
names(prec_cuts) <- c("High","Medium")

cov_cuts = c(0.5,0.25)
names(cov_cuts) <- c("High","Medium")

pool_cuts = c(0.33,0.1)
names(pool_cuts) <- c("High","Medium")

backcast_cuts = c(0.90,0.75)
names(backcast_cuts) <- c("High","Medium")




# build cluster -----------------------------------------------------------

post_diff_mean <- function(x,y){
  x = as.numeric(unlist(x))
  y = as.numeric(unlist(y))
  d <- y-x
  return(mean(d))
}
post_diff_q <- function(x,y,q = 0.5){
  x = as.numeric(unlist(x))
  y = as.numeric(unlist(y))
  d <- y-x
  return(quantile(d,q))
}
p_pos <- function(x,y){
  x = as.numeric(unlist(x))
  y = as.numeric(unlist(y))
  d <- y-x
  return(length(which(d > 0))/length(d))
}


### consider alternate approaches of comparing
### the recent three gen to long-term (including recent period)
### the recent three gen to all previous period
###
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)


test <- foreach(i = rev(c(1:nrow(sp_list))),
                .packages = c("bbsBayes2",
                              "tidyverse",
                              "cmdstanr",
                              "patchwork"),
                .errorhandling = "pass") %dopar%
  {

   #for(i in 1:nrow(sp_list)){
    sp <- as.character(sp_list[i,"english"])
    esp <- as.character(sp_list[i,"french"])
    aou <- as.integer(sp_list[i,"aou"])
    species_f_bil <- gsub(paste(esp,sp),pattern = "[[:space:]]|[[:punct:]]",
                          replacement = "_")


    if(file.exists(paste0(external_dir,"/Indices/Inds_",aou,".rds")) &
       (!file.exists(paste0(external_dir,"/Trends/Accelerating_trends/",aou,"_trend_acceleration.rds")) | re_run)){



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

      ## set three generations
      ## unless < 10, then 10   or  unless > number of years available, then n-years
      gen3 <- min((YYYY-fy),max(10,round(as.numeric(sp_list[i,"GenLength"])*3)))

    gen3_time <- max(10,round(as.numeric(sp_list[i,"GenLength"])*3))

      inds <- readRDS(paste0(external_dir,"/Indices/Inds_",aou,".rds"))


# Annual trend calculations  -------------------------

      ends <- c(seq((YYYY),(fy+gen3),by = -1*gen3))


      if(length(ends) < 2){
        gen3 <- floor((YYYY-fy)/2)
        ends <- c(seq((YYYY),(fy+gen3),by = -1*gen3))

      }
      trend_accel_out <- NULL

      # ind <- readRDS(paste0(external_dir,"/Indices/Ind_plot_",aou,".rds"))
      #
      #
      # trajs_all <- plot_indices(ind,
      #                           title = FALSE,
      #                           ci_width = 0.9,
      #                           add_observed_means = FALSE,
      #                           add_number_routes = FALSE,
      #                           min_year = fy,
      #                           title_size = 12,
      #                           axis_title_size = 10,
      #                           axis_text_size = 10)
      #
      #
      # # trajs <- readRDS(paste0(external_dir,"/Figures/temp_rds_storage/",aou,"_highlevel_simple_trajs.RDS"))
      # traj1 <- trajs_all[["continent"]]
      # n1 <- inds$indices %>%
      #   filter(region == "continent",
      #          year >= fy)
      #
      # strats <- str_split_1(unlist(n1[1,"strata_included"]),
      #                       pattern = " ; ")
      # upy <- max(n1$index_q_0.95,na.rm = TRUE)/2
      #
      #
      # trajs <- traj1 +
      #   geom_ribbon(data = n1, aes(x = year,y = index,ymin = index_q_0.05,ymax = index_q_0.95),
      #               fill = grey(0.5),alpha = 0.2)+
      #   geom_line(data = n1, aes(x = year,y = index),
      #             colour = grey(0.5))
      #

      # pdf(paste0(external_dir,"/trends/rolling_trend_maps/",species_f_bil,"_rolling_trend_map.pdf"),
      #     height = 11,
      #     width = 8.5)

      for(dd in ends[-length(ends)]){
        trends_10temp <- generate_trends(inds,
                                         min_year = dd-gen3,
                                         max_year = dd,
                                         export_full_posterior = TRUE)

        trends_10temp2 <- generate_trends(inds,
                                         min_year = dd-(gen3*2),
                                         max_year = dd-gen3,
                                         export_full_posterior = TRUE)



        t_dif <- trends_10temp$trends_full_posterior %>%
          select(-c(start_year,end_year)) %>%
          rename(trend_full_posterior_1 = trend_full_posterior) %>%
          left_join(trends_10temp2$trends_full_posterior,
                    by = c("region","region_type")) %>%
          select(-c(start_year,end_year)) %>%
          rowwise() %>%
          mutate(trend_dif = post_diff_mean(trend_full_posterior,trend_full_posterior_1),
                 trend_dif_median = post_diff_q(trend_full_posterior,trend_full_posterior_1,0.5),
                 trend_dif_q5 = post_diff_q(trend_full_posterior,trend_full_posterior_1,0.05),
                 trend_dif_q10 = post_diff_q(trend_full_posterior,trend_full_posterior_1,0.1),
                 trend_dif_q90 = post_diff_q(trend_full_posterior,trend_full_posterior_1,0.9),
                 trend_dif_q95 = post_diff_q(trend_full_posterior,trend_full_posterior_1,0.95),
                 p_pos = p_pos(trend_full_posterior,trend_full_posterior_1),
                 trend_years = paste0(dd-gen3,":",dd," - ",dd-(gen3*2),":",dd-gen3)) %>%
          select(-c(trend_full_posterior,trend_full_posterior_1))

        trend1 <- trends_10temp2$trends %>%
          select(region,region_type,
                 starts_with("trend")) %>%
          rename_with(.fn = ~gsub("trend","previous_trend",.x))

        trend2 <- trends_10temp$trends %>%
          left_join(trend1,by = c("region","region_type"))

        t_dif_1 <- t_dif %>%
          left_join(trend2,
                    by = c("region","region_type"))



        trend_accel_out <- bind_rows(trend_accel_out,t_dif_1)
}
trend_accel_out <- trend_accel_out%>%
  mutate(species = sp,
         espece = esp,
         bbs_num = aou,
         three_gen = gen3_time,
         three_gen_used = gen3) %>%
  mutate(across(where(is.double) & !contains("year") &
                  !starts_with("n_") & !starts_with("bbs_num"),~signif(.,3))) %>%
  relocate(species,espece,bbs_num,three_gen,three_gen_used,
           region,region_type,trend,previous_trend,
           starts_with("trend_dif"))







      saveRDS(trend_accel_out, file = paste0(external_dir,"/Trends/Accelerating_trends/",aou,"_trend_acceleration.rds"))

  }


}

parallel::stopCluster(cluster)




avian_core <- readxl::read_xlsx("data/Avian_Core_20251124.xlsx") %>%
  filter(Full_Species__Espèce_complète == "Yes - Oui") %>%
  select(Species_ID_Espèce, BBS_Number__Numéro_BBS, Sort_Order__Ordre_de_tri) %>%
  rename_with(.,.fn = ~paste0(.x,"_core")) %>%
  distinct() %>%
  mutate(aou = as.integer(BBS_Number__Numéro_BBS_core))

# Compile and explore acceleration ----------------------------------------

trend_accell_sum <- NULL

for(i in 1:nrow(sp_list)){

aou <- as.integer(sp_list[i,"aou"])

if(!file.exists(paste0(external_dir,"/Trends/Accelerating_trends/",aou,"_trend_acceleration.rds"))){
  next
  }

trend_accel_tmp <- readRDS(file = paste0(external_dir,"/Trends/Accelerating_trends/",aou,"_trend_acceleration.rds"))


trend_accell_sum <- bind_rows(trend_accell_sum,trend_accel_tmp)

}


trend_accell_high <- trend_accell_sum %>%
  filter(region_type %in% c("continent","country"),
         end_year == 2024) %>%
  mutate(change_cat = ifelse(p_pos >= 0.5 & (previous_trend > 0 & trend > 0),
                             "Accelerating Increase",
                             ""),
         change_cat = ifelse(p_pos < 0.5 & (previous_trend > 0 & trend > 0),
                             "Decelerating Increase",
                             change_cat),
         change_cat = ifelse(p_pos < 0.5 & (previous_trend < 0 & trend < 0),
                             "Accelerating Decrease",
                             change_cat),
         change_cat = ifelse(p_pos >= 0.5 & (previous_trend < 0 & trend < 0),
                             "Decelerating Decrease",
                             change_cat),
         change_cat = ifelse(p_pos < 0.5 & (previous_trend > 0 & trend < 0),
                             "Increase to Decrease",
                             change_cat),
         change_cat = ifelse(p_pos >= 0.5 & (previous_trend < 0 & trend > 0),
                             "Decrease to Increase",
                             change_cat),
         change_cat_full = ifelse(p_pos > 0.85 | p_pos < 0.15,
                                  change_cat,
                                  paste("uncertain",change_cat)),
         change_cat_alpha = ifelse(p_pos > 0.95 | p_pos < 0.05,
                                   1,0.2))


write_excel_csv(trend_accell_high,
                "Three_generation_trend_change_analysis_broad_scale.csv")

three_gen_plot <- ggplot(data = trend_accell_high,
                         aes(x = previous_trend,
                             y = trend_dif_median,
                             colour = change_cat,
                             alpha = change_cat_alpha))+
  geom_abline(slope = c(1,-1),
              intercept = c(0,0))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point()+
  scale_colour_viridis_d(name = "Category of\ntrend difference",
                         option = "turbo")+
  scale_alpha_continuous(name = "Probability of\ndifference",
                         breaks = c(0.21,1),
                         labels = c("< 0.95","> 0.95"))+
  facet_wrap(vars(region))+
  coord_cartesian(xlim = c(-20,20),
                  ylim = c(-20,20))+
  theme_bw()



three_gen_plot

trend_accell_high_can <- trend_accell_high %>%
  filter(region == "Canada")

trend_accell_high_can_labs <- trend_accell_high_can %>%
  filter(change_cat_alpha == 1) %>%
  inner_join(avian_core, by = c("bbs_num" = "aou")) %>%
  mutate(labl = Species_ID_Espèce_core)

three_gen_plot_can <- ggplot(data = trend_accell_high_can,
                         aes(x = previous_trend,
                             y = trend,
                             colour = change_cat))+
  geom_abline(slope = c(1),
              intercept = c(0),
              alpha = 0.5)+
  geom_hline(yintercept = 0,
             alpha = 0.5)+
  geom_vline(xintercept = 0,
             alpha = 0.5)+
  geom_point(aes(alpha = change_cat_alpha))+
  ggrepel::geom_text_repel(data = trend_accell_high_can_labs,
                  aes(label = labl,
                      x = previous_trend,
                      y = trend,
                      colour = change_cat),
                  size = 2.5,
                  min.segment.length = 0)+
  scale_colour_viridis_d(name = "Category of\ntrend difference",
                         option = "turbo")+
  scale_alpha_continuous(name = "Probability of difference\nbetween recent and previous",
                         breaks = c(0.21,1),
                         labels = c("< 0.95","> 0.95"))+
 coord_cartesian(xlim = c(-20,20),
                  ylim = c(-20,20))+
  ylab("Recent three-generation trend")+
  xlab("Previous three-generation trend")+
  theme_bw()


pdf("Canada_trend_change_three_generation.pdf",
    width = 11, height = 8.5)
print(three_gen_plot_can)
dev.off()
