## script to run species analysis on GPSC


library(bbsBayes2)
strat <- "bbs"


aou <- Sys.getenv("aou") # aou index for each species supplied by bash script

# load prepared species data
s <- readRDS(paste0("/home/acs001/BBS_Trends_CWS/prepared_data/",aou,"_data.rds"))

# if there are enough strata to make the spatial model useful
if(nrow(s$meta_strata) > 2){ #spatial models are irrelevant with < 3 strata
  bbs_dat_sp <- prepare_spatial(s,
                                strata_map = load_map(strat),
                                queen = TRUE)

  #print(bbs_dat_sp$spatial_data$map)

  bbs_dat <- bbs_dat_sp %>%
    prepare_model(.,
                  model = "gamye",
                  model_variant = "spatial")

}else{ # else just fit the hierarchical version of the model
  bbs_dat <- prepare_model(s,
                           model = "gamye",
                           model_variant = "hier")
}

# file name for saved output
out_name <- paste0("fit_",aou)
# location for saved output
out_loc <- "/home/acs001/BBS_Trends_CWS/output"

# run model using 4 chains, each sampled in parallel = requires 4 cores
m <- run_model(bbs_dat,
               parallel_chains = 4, # default
               chains = 4, # default
               refresh = 0,
               adapt_delta = 0.9,
               output_basename = out_name,
               output_dir = out_loc)






