## id species with BCRs that have too much missing data.


library(bbsBayes2)
library(tidyverse)

sp_list <- readRDS("species_list.rds") %>%
  filter(model == TRUE)

strats_3 <- c("CA-NT-6N","CA-NT-7W","CA-AB-6N")


for(i in (1:nrow(sp_list))){  # tmp_clr){ #
  sp <- as.character(sp_list[i,"english"])
  aou <- as.integer(sp_list[i,"aou"])


  # identifying first years for selected species ----------------------------
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



  strat <- "bbs"

  s <- try(stratify(by = strat,
                release = 2025,
                species = sp,
                quiet = TRUE,
                distance_to_strata = 4000) %>%
    prepare_data(min_max_route_years = 2,
                 quiet = TRUE,
                 min_year = fy),
    silent = TRUE)
  if(class(s) == "try-error"){next}

  if(any(s$meta_strata$strata_name %in% strats_3)){
    sp_list[i,"northern"] <- TRUE
  }else{
    sp_list[i,"northern"] <- FALSE
  }

}

write_excel_csv(sp_list,"northern_strata_species.csv")

table(sp_list$vm,sp_list$northern)
