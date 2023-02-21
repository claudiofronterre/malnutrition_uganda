# LOAD REQUIRED PACKAGE AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("dplyr", "PrevMap", "sf") # package names
pacman::p_load(pkgs, character.only = T)
  
source("R/functions.R")
  
# LOAD DATA and PROCESS --------------------------------------------------------
  
scores <- c("HAZ", "WAZ", "WHZ")
  
# Districts in Uganda
uga_districts <- st_read("data/original/geodata/gadm36_UGA.gpkg", 
                         layer = "gadm36_UGA_1") %>% 
  st_transform(epsgKM(32635))

country_prev <- list()
  
for (i in scores) {
  
  # Predictions
  if(i == "WAZ") {
    predictions <- readRDS(paste0("output/mix_pred_", i, ".rds"))
  } else {
    predictions <- readRDS(paste0("output/geo_pred_", i, ".rds"))
  }
  
  # Fitted model 
  fit <- readRDS(paste0("output/fit_geo_", i, ".rds"))
    
  # Pop raster
  pop_raster <- readRDS("data/processed/pred_pop_grid.rds")$pop
    
  # DATA PREPARATION -------------------------------------------------------------
    
  # Save prediction locations
  # and convert to sf object
  pred_coords <- predictions$coords %>% 
    as.data.frame() %>% 
    st_as_sf(coords = 1:2, crs = epsgKM(32635))
    
  # Calculate mean_score
  mean_score <- rowMeans(predictions$samples)
    
  # Convert the predicted HAZ to predicted prevalence
  if(i == "WAZ") {
    st_dev <- predictions$st_dev
    predictions <- apply(predictions$mu, 1, 
                         function(x) pnorm((-2 - x) / st_dev))
    
  } else {
    predictions <- apply(predictions$samples, 1, 
                         function(x) pnorm((-2 - x) / sqrt(coef(fit)[["omega^2"]])))
  }
    
  
  predictions <- t(predictions)
    
  # Convert the prevalence to number of malnourished children 
  pop <- extract(pop_raster, pred_coords)
  burden_pop <- predictions * pop
  
  # Summary at country level
  weights <- pop / sum(pop)
  country_prev[[paste0(i, "_unweighted")]] <- colMeans(predictions)
  country_prev[[paste0(i, "_weighted")]] <- colSums(weights * predictions)
  
  # Probability of exceeding country level prevalence
  ex_prob <- rowMeans(predictions > mean(country_prev[[paste0(i, "_unweighted")]]))
  ex_prob_w <- rowMeans(predictions > mean(country_prev[[paste0(i, "_weighted")]]))
  
  # Raster with mean prevalence and mean burden
  pred_rast <- rasterFromXYZ(data.frame(st_coordinates(pred_coords),
                                        mean_burden = rowMeans(burden_pop),
                                        mean_prev = rowMeans(predictions) * 100,
                                        mean_score,
                                        ex_prob, 
                                        ex_prob_w))
    
  saveRDS(pred_rast, paste0("output/mapping/rast_", i, ".rds"))
    
  # Calculate summary at district level
  districts <- as.character(st_join(pred_coords, uga_districts)$NAME_1)
  burden_pop <- split(as.data.frame(burden_pop), f = districts)  
  burden_district <- burden_pop %>% 
    lapply(function(x) data.frame(burden = round(mean(colSums(x))))) %>% 
    bind_rows(.id = "NAME_1")
    
  map_summary <- left_join(uga_districts, burden_district)
  map_summary$burden[grep("Lake", map_summary$NAME_1)] <- NA
    
  saveRDS(map_summary, paste0("output/mapping/district_", i, ".rds"))
  print(i)
  
}

saveRDS(country_prev, "output/country_prev.rds")