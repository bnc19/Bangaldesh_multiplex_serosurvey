
# Fit inla model with spatial effects to multiple pathogens 

fit_spatial_models = function(formula,
                              model_data,
                              spatial_objects,
                              pathogens,
                              model_name,
                              variable_labels = NULL,
                              warnings = F,
                              spatial = T,
                              formula_chik) {
  

  my_theme = 
    theme_minimal() +
    theme(
      legend.text = element_text(size = 7),
      legend.spacing.y = unit(0, "cm"),  # Reduced to zero
      legend.spacing.x = unit(0, "cm"),  # Added to control horizontal spacing
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 7),
      legend.title = element_blank(),
      plot.title = element_text(size = 8),
      strip.text = element_text(size = 7)) 
  
  # Initialize results
  results = vector("list", length(pathogens))
  
  # Fit models for each pathogen
  for (i in 1:length(pathogens)) {
   
     p = pathogens[i]

     # Filter pathogen-specific data
     data = model_data %>%
       filter(pathogen == !!p) %>% 
       arrange(sample_id)
     
     cov = data %>% 
       select(-sample_id, -pathogen, -lon, -lat, -value) 
    
      
    # Create data stack
    stack = inla.stack(
      data = list(y = data$value),
      A = list(spatial_objects$A.spde, 1), # project spatial field to observed locations 
      effects = list( # fixed and random effects 
        spatial_objects$spatial.index, 
        cov
      ),
      tag = "data"
    )
    
    # Fit model  
    
    if(p == "CHIKV"){ # include regularising priors for CHIKV
    model = inla(
      formula_chik,
      data = inla.stack.data(stack),
      family = "binomial",
      control.predictor = list(A = inla.stack.A(stack)),
      control.compute = list(dic = TRUE, waic = TRUE),
      control.fixed = list(
        prec = list(default = 0.44) 
      )
      )
    } else {
      model = inla(
        formula,
        data = inla.stack.data(stack),
        family = "binomial",
        control.predictor = list(A = inla.stack.A(stack)),
        control.compute = list(dic = TRUE, waic = TRUE),
        control.fixed = list(
          prec = list(default = 0.44)  
        )
      ) 
    }

   if(spatial == T) {
      # Extract spatial parameters
      spde_result = inla.spde2.result(model, "spatial", spatial_objects$spde)
      spatial_summary = list(
        range = spde_result$summary.log.range.nominal, # How far spatial correlation extends (km)
        variance = spde_result$summary.log.variance.nominal # How much spatial variation exists
      )
 
   }     else {
     spatial_summary = NULL
   }
    
    data_idx = inla.stack.index(stack, tag = "data")$data 
    
    
    # Store results
    results[[i]] = list(
      fixed = model$summary.fixed,
      random = model$summary.hyperpar,
      waic = model$waic$waic,
      dic = model$dic$dic,
      spatial_params = spatial_summary,
      pos = model$summary.fitted.values$mean[data_idx],
      model_object = model
    )
    
    print(paste(p, "completed"))
  }
  
  # Process results
  fixed = map(results, "fixed")
  random = map(results, "random")
  waic = map(results, "waic")
  dic = map(results, "dic")
  if(spatial == T) {
  spatial_params = map(results, "spatial_params")
  } else {
    spatial_params = NULL
  }
  pos = map(results, "pos")
  
  # Clean fixed effects
  fixed = lapply(fixed, rownames_to_column)
  fixed_clean = bind_rows(fixed, .id = "pathogen") %>%
    mutate(pathogen = factor(pathogen, labels = pathogens))
  
  # Clean random effects
  random_clean = NULL
  if (!is.null(random) &&
      length(random) > 0 && !is.null(random[[1]])) {
    random = lapply(random, rownames_to_column)
    random_clean = bind_rows(random, .id = "pathogen") %>%
      mutate(pathogen = factor(pathogen, labels = pathogens)) %>%
      rename(lower = `0.025quant`, upper = `0.975quant`) %>%
      select(pathogen, rowname, mean, lower, upper) %>%
      mutate(across(c(mean, lower, upper), ~ round(.x, 2)))
  }
  
  # Clean WAIC and DIC
  waic_clean = data.frame(
    pathogen = factor(1:length(waic), labels = pathogens),
    waic = unlist(waic),
    model = model_name
  )
  
  dic_clean = data.frame(
    pathogen = factor(1:length(dic), labels = pathogens),
    dic = unlist(dic),
    model = model_name
  )
  
  
  if(spatial == T) {
  
  # Clean spatial parameters
  spatial_clean = NULL
  if (!is.null(spatial_params[[1]])) {
    spatial_clean = map_dfr(spatial_params, ~ {
      if (!is.null(.x)) {
        data.frame(
          range_mean_km = exp(.x$range$mean),
          range_lower_km = exp(.x$range$`0.025quant`),
          range_upper_km = exp(.x$range$`0.975quant`),
          variance_mean = exp(.x$variance$mean),
          variance_lower = exp(.x$variance$`0.025quant`),
          variance_upper = exp(.x$variance$`0.975quant`)
        )
      } else {
        data.frame(
          range_mean_km = NA,
          range_lower_km = NA,
          range_upper_km = NA,
          variance_mean = NA,
          variance_lower = NA,
          variance_upper = NA
        )
      }
    }, .id = "pathogen") %>%
      mutate(pathogen = factor(pathogen, labels = pathogens)) %>%
      mutate(across(starts_with("range"), ~ round(.x, 1))) %>%
      mutate(across(starts_with("variance"), ~ round(.x, 3)))
  }
  
  } else{
    spatial_clean = NULL 
    
  }
  # Clean p(pos)
  
pos_clean = pos %>% 
    bind_cols() %>% 
  mutate(id = data$sample_id)

colnames(pos_clean) = c(pathogens, "sample_id")

  
  return(
    list(
      fixed_effects = fixed_clean,
      random_effects = random_clean,
      spatial_parameters = spatial_clean,
      waic = waic_clean,
      dic = dic_clean,
      prob_pos = pos_clean,
      raw_results = results
    )
  )
}



