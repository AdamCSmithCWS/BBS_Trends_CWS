## creating the data file for each species so that HPC does not require the
## bbsBayes::fetch_bbs_data()
## run locally then upload full prepared_data folder to HPC
##

library(bbsBayes2)
library(tidyverse)

output_dir <- "e://BBS_Trends_CWS/output"


sp_list <- readRDS("species_list.rds")
strat <- "bbs"

strat_alt <- load_map(strat)


# Northern Strata that are not worth including
strats_3 <- c("CA-MB-70", "CA-NL-69", "CA-NT-68", "CA-NT-70",
              "CA-NU-69", "CA-NU-68", "CA-NU-70",
              "CA-QC-69", "CA-QC-68", "CA-QC-70", "CA-YT-70",
              "US-AK-69")


 strat_alt <- strat_alt |>
  filter(!strata_name %in% strats_3)

for(i in 1:nrow(sp_list)){



  sp <- as.character(sp_list[i,"english"])
  aou <- as.integer(sp_list[i,"aou"])


  fy <- NULL
  if(aou %in% c(4661,4660)){ #Alder and Willow Flycatcher
    fy <- 1978 #5 years after the split
  }
  if(aou %in% c(10,11,22860)){ # Clark's and Western Grebe and EUCD
    fy <- 1990 #5 years after the split and first year EUCD observed on > 3 BBS routes
  }
  if(aou == 6121){ # CAve Swallow
    fy = 1985
  }



  s <- try(stratify(by = strat,
                release = 2026,
                species = sp,
                quiet = TRUE,
                distance_to_strata = 4000) |>
  prepare_data(min_max_route_years = 2,
               quiet = TRUE,
               min_year = fy),
  silent = TRUE)

  if(class(s) == "try-error"){
    print(paste("not enough routes with",sp))
    next}

  if(any(s$meta_strata$strata_name %in% strats_3)){

    strat_alt <- load_map(strat) |>
      filter(!strata_name %in% strats_3)

    s <- stratify(by = strat,
                  strata_custom = strat_alt,
                  release = 2026,
                  species = sp,
                  quiet = TRUE,
                  distance_to_strata = 4000)  |>
      prepare_data(min_max_route_years = 2,
                   quiet = TRUE,
                   min_year = fy)

  }
  ## bbsBayes2 models do not currently work unless n_strata > 1
  if(nrow(s$meta_strata) == 1){
    warning(paste("Only 1 stratum for",sp,"skipping to next species"))
    next
  }
# data suitable for running models in the GPSC
  saveRDS(s,file = paste0("prepared_data/",aou,"_data.rds"))

  if(nrow(s$meta_strata) > 2){ #spatial models are irrelevant with < 3 strata
    bbs_dat_sp <- prepare_spatial(s,
                                  strata_map = load_map(strat),
                                  queen = TRUE)

    #print(bbs_dat_sp$spatial_data$map)

    saveRDS(bbs_dat_sp,paste0("raw_data/spatial_neighbours_",aou,".rds"))
    bbs_dat <- bbs_dat_sp |>
      prepare_model(.,
                    model = "gamye",
                    model_variant = "spatial")

  }else{
    bbs_dat <- prepare_model(s,
                             model = "gamye",
                             model_variant = "hier")
  }


  saveRDS(bbs_dat,file = paste0("prepared_data_local/",aou,"_data.rds"))

  print(round(i/nrow(sp_list),2))
}


 completed <- list.files(output_dir) |>
   str_extract(pattern = "(?<=_)[[:digit:]]{2,6}")
 #completed <- as.integer(completed[-which(is.na(completed))])

 ## write text file with just aous for bash script in gpsc
 # aous <- sp_list[c((nrow(sp_list)-202):(nrow(sp_list)-101)),"aou"]
 for(i in 1:nrow(sp_list)){



   sp <- as.character(sp_list[i,"english"])
   aou <- as.integer(sp_list[i,"aou"])
if(file.exists(paste0("prepared_data/",aou,"_data.rds"))){
s <- readRDS(paste0("prepared_data/",aou,"_data.rds"))

   if(nrow(s$meta_strata) == 2){
     sp_list[i,"rerun"] <- TRUE
   }
}

 }


 #aous <- sp_list[c((nrow(sp_list)-302):(nrow(sp_list)-203)),"aou"]
 aous <- sp_list[which(sp_list$rerun),"aou"]

 aous <- aous[-which(aous$aou %in% completed),] |>
   distinct()

 write_tsv(aous, "aou_list.txt",
           col_names = FALSE)


