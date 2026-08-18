## script to run species analysis on GPSC


library(bbsBayes2)
strat <- "bbs"


aou <- Sys.getenv("aou") # aou index for each species supplied by bash script
# location for saved output
out_loc <- "/gpfs/fs7/eccc/esrp/cws/acs001/BBS_Trends_CWS/output"

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

# run model using 4 chains, each sampled in parallel = requires 4 cores
m <- run_model(bbs_dat,
               parallel_chains = 4, # default
               chains = 4, # default
               refresh = 0,
               adapt_delta = 0.9,
               iter_warmup = 2000,
               iter_sampling = 6000,
               thin = 6,
               output_basename = out_name,
               output_dir = out_loc,
               save_model = FALSE)

               #output_dir = out_loc)


save_model_run <- function(model_output,
                           retain_csv = TRUE, path = NULL, quiet = FALSE,
                           save_file_path = NULL) {


  model_fit <- model_output$model_fit

  if(is.null(path)) {
    #if(!retain_csv){
    csv_path <- model_fit$output_files()


    if(any(!file.exists(csv_path))) {


      stop("Cannot find original model file location, please specify `path`",
           call. = FALSE)
    }
    #}

    path <- csv_path %>%
      normalizePath() %>%
      stringr::str_remove("-[0-9]{1,3}.csv$") %>%
      unique() %>%
      paste0(".rds")

    if(is.null(save_file_path)){
      save_file_path <- path
    }else{
      #check_dir(dirname(save_file_path))
      if(ext(save_file_path) != "rds") {
        stop("save_file_path must have a .rds extension", call. = FALSE)
      }
    }

    if(!quiet) message("Saving model output to ", save_file_path)
  } else {
    csv_path <- model_fit$output_files()
    if(!is.null(save_file_path)){
      #check_dir(dirname(save_file_path))
    }else{
      save_file_path <- path
    }
    if(ext(save_file_path) != "rds") {
      stop("save_file_path must have a .rds extension", call. = FALSE)
    }
  }

  # Ensure all lazy data loaded (see ?cmdstanr::save_object)
  model_fit$draws()
  try(model_fit$sampler_diagnostics(), silent = TRUE)
  try(model_fit$init(), silent = TRUE)
  try(model_fit$profiles(), silent = TRUE)

  # Update entire model output object and save
  model_output[["model_fit"]] <- model_fit
  readr::write_rds(model_output, save_file_path)

  if(!retain_csv){

    unlink(csv_path) # deleting the csv files
  }
  invisible(model_output)
}

ext <- function(file) {
  stringr::str_extract(file, "(?<=\\.)[[:alnum:]]+$")
}

save_model_run(m,
               retain_csv = FALSE,
               path = paste0(out_loc,
                             "/",
                             out_name,
                             ".rds"))



summ <- get_summary(m)
saveRDS(summ,paste0(out_loc,"/Convergence/summ_",aou,".rds"))

# saving raw data locally
raw_data <- m$raw_data

saveRDS(raw_data,paste0(out_loc,"/Raw_data/Raw_",aou,".rds"))




