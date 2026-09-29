# Script to predict the seroprevalence of unsampled communities using binomial
# spatial models (via INLA), for each of five pathogens (CHIKV, DENV, JEV, CHOL, HEV).


# Input:  data/deidentified_data.csv - long-format serosurvey data, one row per
#         individual per pathogen, with binary infection status (`value`). 
#         data/spatial_objects.RDS - spatial objects for INLA. 
# Output: outputs/cv_spatial/cv_results_list_hold_out_size_%d.RDS
#         outputs/spatial_performance_plot.jpg
# 
# Method:
#
# Between 1 and 65 communities are heldout, and the model is trained on the 
# remaining communities. At each hold-out size the communities to withhold are 
# drawn at random and a spatial-only model, is fitted to the remaining communities
# for each pathogen. Seroprevalence is predicted in the withheld communities. 
# This is repeated over n_iterations random draws of the hold-out set. 
# Predicted and observed community seroprevalence are compared using mean
# absolute error, the correlation across communities, and the Brier skill score. 
# The three metrics are plotted against the number of communities used for training.

################################################################################
# Note: results will differ from those in the manuscript because the coordinates
# released here are rounded, which blurs the finer distance bands, and the 
# script is set up to only run 5 iterations. 
################################################################################


# set up -----------------------------------------------------------------------


library(INLA)
library(tidyverse)
library(purrr)
library(patchwork)

# Set up simulation parameters
set.seed(14)
n_iterations = 5 # 100 to replicate results in manuscript

source("R/functions.R")

cols = c("#CD3572", "#E0BBC7",  "#FD8D3C",
         "#1E3D8A", "#8CA1CC")

# read in data 
model_data = read.csv("data/deidentified_data.csv") |>  select(-X) |> 
  mutate(age_group = factor(age_group, levels = c("<10","11-20", "21-30", "31-40", "41-50", "51-60", ">60"))) 

spatial_objects = readRDS(file = "data/spatial_objects.RDS")

df  = model_data %>% 
  select(community_id, sample_id, pathogen, value, lon, lat, logpopden)

communities  = unique(df$community_id)

pathogens = unique(df$pathogen)

# stablilise chikungunya model fitting 
spde_prior = inla.spde2.pcmatern(
  mesh = spatial_objects$mesh,
  prior.range = c(20, 0.5),     
  prior.sigma = c(0.5, 0.05)      
)

folder = "cv_spatial"

dir.create(paste0("outputs/", folder), recursive = T, showWarnings = F)

sizes <- c(1, 10, 30, 50, 65) # number of communities to holdout 

# Run cross-validation ---------------------------------------------------------
spatial_formula <- y ~ f(spatial, model = spatial_objects$spde)

chik_spatial_formula <- y ~ f(spatial, model = spde_prior)

cv_results <- lapply(sizes, function(h) {
  run_seroprev_cv(model_data      = df,
                  pathogens       = pathogens,
                  n_iterations    = n_iterations,
                  hold_out_size   = h,
                  spatial_objects = spatial_objects,
                  main_formula    = spatial_formula,
                  chik_formula    = chik_spatial_formula,
                  communities     = communities,
                  folder          = folder,
                  spde_prior      = spde_prior)
})

names(cv_results) <- paste0("holdout_", sizes)

# Plot cross-validation performance  -------------------------------------------

# read in model results  

spatial <- sizes |> 
  sprintf(fmt = paste0("outputs/", folder, "/cv_results_list_hold_out_size_%d.RDS")) |> 
  setNames(paste0("spatial_", sizes)) |> 
  lapply(readRDS)

# spatial only 
spatial_performance = lapply(spatial, calculate_performance_metrics) %>% 
  bind_rows(.id = "holdout") %>% 
  pivot_longer(cols = MAE:Correlation) %>% 
  mutate(communities = 70 - as.numeric(str_remove(holdout, "spatial_")))

# Create individual plots
p1 <- make_sub_plot("BSS") + geom_hline(yintercept = 0, linetype = 2) + ylab("Brier skill score")
p2 <- make_sub_plot("MAE")  + ylab("Mean absolute error")
p3 <- make_sub_plot("Correlation") 

# Combine side-by-side with a shared legend at the top
spatial_performance_plot <- (p1 | p2 | p3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top", legend.title = element_blank())

ggsave(spatial_performance_plot, 
       filename = "outputs/spatial_performance_plot.jpg",
       units = "cm",
       height = 6,
       width = 12)

