# Script to calculate the degree of spatial clustering at different distance 
# scales for each pathogen. Takes as input the data frame 'deidentified_data.csv'.  

# set up -----------------------------------------------------------------------
library(tidyverse)
library(sf) 
library(units)


my_theme = 
  theme(
    legend.text = element_text(size = 7),
    legend.spacing.y = unit(0, "cm"),  # Reduced to zero
    legend.spacing.x = unit(0, "cm"),  # Added to control horizontal spacing
    legend.box.spacing = unit(0, "pt") ,  # space between plot panel and legend box
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.1, "lines"),
    legend.title = element_text(size = 7),
    legend.key.size = unit(5, "pt"),
    strip.text = element_text(size = 7))

cols = c("#CD3572", "#E0BBC7",  "#FD8D3C",
         "#1E3D8A", "#8CA1CC")

pathogens =  c("CHIKV", "DENV", "JEV", "CHOL", "HEV")

# Set seed for reproducibility
set.seed(123)

# Number of bootstrap iterations
n_boot = 100  # 100 to replicate results in manuscript

# get inf status ---------------------------------------------------------------
data_long = read.csv("data/deidentified_data.csv")

data_long = data_long %>% 
  mutate(pathogen = factor(pathogen, levels = pathogens))

inf = data_long %>% 
  select(pathogen, sample_id, value) %>% 
  pivot_wider(names_from = pathogen) %>% 
  arrange(sample_id) %>% 
  select(-sample_id)

# calculate distance between all individuals -----------------------------------
coords = data_long %>%
  select(lat, lon, sample_id, hhid, community_id) %>%
  distinct() %>%
  arrange(sample_id)

ind_sf = st_as_sf(coords, coords = c("lon", "lat"),  crs = 4326)

# calculate matrix in KM and drop units
dist_matrix = drop_units(st_distance(ind_sf) / 1000) # distance matrix of pairs of people 

# mask diagonal 
diag(dist_matrix) = NA

# Create household matrix 
hh_matrix = outer(coords$hhid, coords$hhid, "==")

# Define distance bands including household
dmax = c(0.01, 0.1, 1, 10, 50, 100, 1e5, NA)  # NA for household
dmin = c(0, 0.01, 0.1, 1, 10, 50, 100, NA)     # NA for household

band_names = c("0-0.01", "0.01-0.1", "0.1-1", "1-10", "10-50", 
               "50-100", ">100", "Household")

# calculate infection concordance by pathogen and distance ---------------------

list_results = list()

for(j in pathogens){
  
  infvec = inf[, which(colnames(inf) == j), drop = TRUE]
  inf_ind = which(infvec == 1)
  n_infected = length(inf_ind)
  n_total = length(infvec)
  
  # results dataframe
  results_df = data.frame(
    pathogen = j, 
    dmin = dmin, 
    dmax = dmax,
    band = band_names,
    concordance = NA, 
    rr = NA,
    rr_mean = NA,
    rr_lower = NA,
    rr_upper = NA,
    n_both_pos = NA 
  )

  # max distance concordance
  distmax = dist_matrix >= 100
  
  max_concordance = sum(distmax[inf_ind, inf_ind], na.rm = T) / 
    sum(distmax[inf_ind, ], na.rm = T)

  # Bootstrap array to store RR values
  boot_rr = matrix(NA, nrow = n_boot, ncol = length(dmin))
  
  # bootstap loop
  
  for(b in 1:n_boot){
    
    # Resample all individuals with replacement
    boot_indices = sample(1:n_total, n_total, replace = TRUE)
    
    # Get infection status for bootstrapped sample
    boot_infvec = infvec[boot_indices]
    boot_inf_ind = which(boot_infvec == 1)
    
    # Resample distance and household matrices
    boot_dist_matrix = dist_matrix[boot_indices, boot_indices]
    boot_hh_matrix = hh_matrix[boot_indices, boot_indices]
    
    # Bootstrap reference concordance
    boot_distmax = boot_dist_matrix >= 100
    boot_max_concordance = sum(boot_distmax[boot_inf_ind, boot_inf_ind], na.rm = TRUE) / 
      sum(boot_distmax[boot_inf_ind, ], na.rm = TRUE)
    
    # Calculate concordance for each distance band
    for(i in 1:(length(dmin) - 1)){
      
      dmat2 = (boot_dist_matrix >= dmin[i]) * (boot_dist_matrix < dmax[i])
    
      boot_concordance = sum(dmat2[boot_inf_ind, boot_inf_ind], na.rm = TRUE) / sum(dmat2[boot_inf_ind, ], na.rm = TRUE)
      
      boot_rr[b, i] = boot_concordance / boot_max_concordance
    }
  
    # Household concordance
      boot_hh_concordance = sum(boot_hh_matrix[boot_inf_ind, boot_inf_ind], na.rm = TRUE) / sum(boot_hh_matrix[boot_inf_ind, ], na.rm = TRUE)
      
      boot_rr[b, length(dmin)] = boot_hh_concordance / boot_max_concordance
  
  }
  
  
  # Calculate overall estimates
  for(i in 1:(length(dmin) - 1)){
    dmat2 = (dist_matrix >= dmin[i]) * (dist_matrix < dmax[i])
    
    results_df$concordance[i] = sum(dmat2[inf_ind, inf_ind], na.rm = TRUE) / 
      sum(dmat2[inf_ind, ], na.rm = TRUE)
    
    results_df$n_both_pos[i] = sum(dmat2[inf_ind, inf_ind], na.rm = TRUE)
    results_df$rr[i] = results_df$concordance[i] / max_concordance
    
    # Calculate 95% CI from bootstrap
    results_df$rr_mean[i] = mean(boot_rr[, i], na.rm = TRUE)
    results_df$rr_lower[i] = quantile(boot_rr[, i], 0.025, na.rm = TRUE)
    results_df$rr_upper[i] = quantile(boot_rr[, i], 0.975, na.rm = TRUE)
  }
  
  # Household concordance
  results_df$concordance[length(dmin)] = sum(hh_matrix[inf_ind, inf_ind], na.rm = TRUE) / 
    sum(hh_matrix[inf_ind, ], na.rm = TRUE) 
  
  results_df$n_both_pos[length(dmin)] = sum(hh_matrix[inf_ind, inf_ind], na.rm = TRUE)
  results_df$rr[length(dmin)] = results_df$concordance[length(dmin)] / max_concordance
  
  # Household CI
  results_df$rr_mean[length(dmin)] = mean(boot_rr[,  length(dmin)], na.rm = TRUE)
  results_df$rr_lower[length(dmin)] = quantile(boot_rr[, length(dmin)], 0.025, na.rm = TRUE)
  results_df$rr_upper[length(dmin)] = quantile(boot_rr[, length(dmin)], 0.975, na.rm = TRUE)
  
  list_results[[j]] = results_df
}

# combine results, format and plot ---------------------------------------------

results_all = list_results %>% bind_rows(.id = "pathogen")

band_levels = c("HH", "0-0.01", "0.01-0.1", "0.1-1", "1-10", "10-50", 
                "50-100", ">100")

band_labels = c("Household", "0-10m", "10-100m", "100m-1km", "1-10km", "10-50km", 
                "50-100km", "100km+")

out = results_all %>% 
  filter(band != "0-0.01") %>% 
  mutate(band = ifelse(band == "Household", "HH", band)) %>% 
  mutate(band = factor(band, levels = band_levels, labels = band_labels )) %>%  
  filter(n_both_pos > 2) %>% 
  ggplot(aes(x = band, y = rr_mean, color = pathogen, group = pathogen)) +
  geom_point(size = 0.8, position = position_dodge(.4)) +
  geom_line(position = position_dodge(.4)) +
  geom_errorbar(aes(ymin = rr_lower, ymax = rr_upper),position = position_dodge(.4), 
               width = 0.1) +
  theme_classic() +
  my_theme + theme(legend.position = "none",
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = cols)  +
  labs(
    x = "Spatial scale",
    y = "Relative risk of infection compared to >100km",
  ) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_log10(breaks = c(1, 2, 5, 20)) +
  facet_wrap(~pathogen, ncol = 1)


ggsave(out, 
       filename = "outputs/spatial_clustering.jpg",
       units = "cm",
       height = 12,
       width = 9)   

