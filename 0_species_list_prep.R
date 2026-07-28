### setup up annual bbs analyses
## generate list of species with simple data summaries
## identify which species to run and sort them into parallel-run groups (column "vm")

if(packageVersion("bbsBayes2") != "1.2026.2"){
install.packages("bbsBayes2",
                 repos = c(bbsbayes = 'https://bbsbayes.r-universe.dev',
                           CRAN = 'https://cloud.r-project.org'))



  pckgs <- c("cmdstanr",
             "sf",
             "tidyverse",
             "doParallel",
             "foreach")


  install.packages(pckgs)

  }



library(bbsBayes2)
library(tidyverse)

all <- load_bbs_data(release = 2026)
species_list <- all$species %>%
  filter(unid_combined == TRUE,
         !grepl("unid",english),
         !grepl("Unid",english))

bird <- all$birds %>%
  select(aou,route_data_id,species_total) %>%
  inner_join(.,species_list,
            by = "aou")

route <- all$routes

sp_sum <- route %>%
  select(country_num,state_num,route,route_name,bcr,
         year,state,country,route_data_id) %>%
  inner_join(.,bird,
            by = "route_data_id")


sp_list <- sp_sum %>%
  group_by(aou,english,french,species) %>%
  summarise(n_obs = n(),
            n_birds = sum(species_total),
            n_routes = length(unique(route_name)),
            n_years = length(unique(year))) %>%
  arrange(n_routes,n_obs,n_years) %>%
  mutate(model = ifelse((n_obs > 30 &
                        n_years > 30 &
                        n_routes > 2),
                        TRUE,FALSE))

#split the species list into groups to spread across many processors
#
# sp_list_mod <- sp_list %>%
#   filter(model)
# sp_list_mod[,"vm"] <- rep(1:10,length.out = nrow(sp_list_mod)) # setting a permanent list of which species go to which vms


sp_list_mod <- sp_list |>
  filter(model) #|>
  #arrange(n_routes) |>
  #rowwise() |>

tmp <- quantile(sp_list_mod$n_obs,seq(0,1,by = 0.1))
tmp[1] <- tmp[1]-1
tmp[length(tmp)] <- tmp[length(tmp)]+1
sp_list_mod <- sp_list_mod |>
  mutate(vm = as.integer(cut(n_obs,breaks = tmp)))





saveRDS(sp_list_mod, "species_list.rds")






