# Script to fit spatial binomial regression models using R INLA.
# Takes as input the data frame 'deidentified_data.csv'.

# set up -----------------------------------------------------------------------

source("R/functions.R")

library(tidyverse)
library(INLA)
library(sf)
library(rnaturalearth) # for map
library(fmesher)

# read in data
model_data = read.csv("data/deidentified_data.csv")[,-1]

cols = c("#CD3572", "#E0BBC7", "#FD8D3C", "#1E3D8A", "#8CA1CC")

# format data

model_data  = model_data %>%
  mutate(travel = factor(travel, levels = c(4, 1, 2, 3)),
         age_group = factor(
           age_group,
           levels = c("<10", "11-20", "21-30", "31-40", "41-50", "51-60", ">60")
         ))

length(unique(model_data$sample_id)) # 2878

# labels
pathogens = unique(model_data$pathogen)

fixed_labels = c(
  "Intercept",
  "11-20",
  "21-30",
  "31-40",
  "41-50",
  "51-60",
  ">60",
  "Male",
  "7 days",
  "7 days - 1 month",
  "1-6 months",
  "Pigs",
  "Log population density",
  "Aeg",
  "Alb"
)

# check data
model_data %>%
  summarise(across(everything(), ~ sum(is.na(.))))

# Prepare spatial data ---------------------------------------------------------

# Create spatial data (not long pathogen format)
spatial_data = model_data %>%
  pivot_wider(names_from = pathogen, values_from = value) %>%
  dplyr::select(hhid, community_id, lat, lon, sample_id) %>%
  arrange(sample_id)

# Convert to UTM
utm_zone = 46
coords_sf = st_as_sf(spatial_data, coords = c("lon", "lat"), crs = 4326)
utm_crs = paste0("+proj=utm +zone=", utm_zone, " +north +datum=WGS84")
coords_utm = st_transform(coords_sf, crs = utm_crs)
utm_coords = st_coordinates(coords_utm)

# Add UTM coordinates
spatial_data$utm_easting = utm_coords[, 1]
spatial_data$utm_northing = utm_coords[, 2]

# Convert to km for mesh creation
spatial_data$utm_easting_km = spatial_data$utm_easting / 1000
spatial_data$utm_northing_km = spatial_data$utm_northing / 1000

coords_km = as.matrix(spatial_data[, c("utm_easting_km", "utm_northing_km")])

# boundary
bnd = inla.nonconvex.hull(coords_km, convex = -0.05)  # Negative value = tighter fit

mesh1 = inla.mesh.2d(
  loc = coords_km,
  boundary = bnd,
  max.edge = c(30, 50), # 30km inner, 50km outer
  cutoff = 2  # if households are less than 2km apart, builds only a single vertex
  ) 


# Define SPDE model (Matérn covariance)
spde = inla.spde2.matern(mesh = mesh1,
                         alpha = 2,
                         constr = TRUE)

# Create projection matrix for household locations
A.spde = inla.spde.make.A(mesh = mesh1, loc = coords_km)

# Create spatial index
spatial.index = inla.spde.make.index("spatial", n.spde = spde$n.spde)

# combine all spatial objects
spatial_objects = list(
  mesh = mesh1,
  spde = spde,
  A.spde = A.spde,
  spatial.index = spatial.index
)

# priors for chikungunya to stabilise estimates 
spde_prior = inla.spde2.pcmatern(
  mesh = mesh1,
  prior.range = c(20, 0.5),
  prior.sigma = c(0.5, 0.05)
)

prior_hh = "hyper = list(prec = list(prior = 'pc.prec', param = c(0.5, 0.02)))"
prior_comm = "hyper = list(prec = list(prior = 'pc.prec', param = c(.8, 0.02)))"

# Fit models with spatial and both random effects  -----------------------------

formula_full = y ~  age_group + sex + travel  +
  pigs + logpopden + aeg + alb +
  # Use priors for all
  f(hhid, model = "iid") +
  f(community_id, model = "iid") +
  f(spatial, model = spatial_objects$spde)

formula_chik = y ~  age_group + sex + travel  +
  pigs + logpopden + aeg + alb +
  f(hhid, model = "iid", hyper = prior_hh) +
  f(community_id, model = "iid", hyper = prior_comm) +
  f(spatial, model = spde_prior)

results_full = fit_spatial_models(
  formula_full,
  model_data,
  spatial_objects,
  pathogens,
  "spatial and random effects",
  variable_labels = fixed_labels,
  warnings = F,
  spatial = T,
  formula_chik = formula_chik
)

multivariate_fixed_effects_df = results_full$fixed_effects

# Fit univariate model with spatial and random effects -------------------------


predictor_vars = c("age_group", "sex", "travel", "pigs", "logpopden", "aeg", "alb")

# data frame to store summary statistics
univariate_summary = data.frame(
  variable = character(),
  model_name = character(),
  stringsAsFactors = FALSE
)

univariate_results = list()

for (var in predictor_vars) {
  formula_univariate = as.formula(
    paste0(
      "y ~ ",
      var,
      " + ",
      "f(hhid, model = 'iid') + ",
      "f(community_id, model = 'iid') +",
      "f(spatial,  model = spatial_objects$spde)"
    )
  )
  
  
  formula_univariate_chik = as.formula(
    paste0(
      "y ~ ",
      var,
      " + ",
      "f(hhid, model = 'iid', ",
      prior_hh,
      ") + ",
      "f(community_id, model = 'iid', ",
      prior_comm,
      ") +",
      "f(spatial,  model = spde_prior)"
    )
  )
  
  univariate_results[[var]] = fit_spatial_models(
    formula_univariate,
    model_data,
    spatial_objects,
    pathogens,
    model_name = paste("Univariate:", var),
    warnings = F,
    spatial = T,
    formula_chik = formula_univariate_chik
  )
  
  print(paste0(var, " completed"))
}


# Extract fixed effects for all variables 
fixed_effects_list  = lapply(univariate_results, function(model) {
  model$fixed_effects
})

univariate_fixed_effects_df = fixed_effects_list %>%
  bind_rows()
  

# save estimates ---------------------------------------------------------------
write.csv(univariate_fixed_effects_df, "outputs/binomial_univariate_results.csv")
write.csv(multivariate_fixed_effects_df,"outputs/binomial_multivariate_results.csv")
