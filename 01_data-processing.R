# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("sf", "dplyr", "raster") # package names
pacman::p_load(pkgs, character.only = T)

if(all(c("dplyr", "raster") %in% .packages())) {
  select <- dplyr::select
}

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------


# Malnutrition data
df <- readr::read_csv("data/original/Anthropometry_2016_UGDHS_coords.csv")

# Shapefile for Uganda
uga <- st_read("data/original/geodata/gadm36_UGA.gpkg", layer = "gadm36_UGA_0")

# DATA PREPARATION -------------------------------------------------------------

# Select essential columns 
# Set to long format
# Kepp only nutrition scorse between -6 and +6
# Remove NA
df_clean <- df %>% 
  select(ID = Cluster_ID, bmonth = Birth_month, byear = Birth_year,
         age = Age_months, gender = Sex_Male_is_1,
         long = LONGNUM, lat = LATNUM, HAZ, WAZ, WHZ) %>% 
  filter(gender >= 1) %>%
  mutate(gender = as.factor(gender)) %>% 
  tidyr::gather(key = "metric", value = "score", 
                -ID, -bmonth, -byear, -long, -lat, -age, -gender) %>% 
  mutate(score = score / 100) %>% 
  filter(score >= -6, score <= 6, long != 0) %>% 
  na.omit()

# Replace cluster ID with ordered numeric ID
df_clean$ID <- as.numeric(as.factor(df_clean$ID))
df_clean <- df_clean[order(df_clean$ID), ]


# Convert the LAT LONG coordinates to UTM (km) EPSSG: 32365
crs_utm_km <- epsgKM(32635)
df_sp <- df_clean %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  st_transform(crs = crs_utm_km)

# Add two new columns with the coordinates in UTM (km)
df_clean[, c("utm_x", "utm_y")] <- st_coordinates(df_sp)

df_clusters <- df_clean %>% 
  select(ID, long, lat, utm_x, utm_y) %>% 
  distinct()

# CREATE PREDICTION GRID -------------------------------------------------------

# Convert to the reference system of the points
uga <- st_transform(uga, crs = crs_utm_km)

# Load population raster
pop <- readRDS("data/original/geodata/U51KM_Total.rds")

# Create prediction grid
pred <- create_grid(resolution = 1, study_area = uga, pop = pop, 
                    cutoff = 0, filename = "data/processed/pred_pop_grid.rds")

# Load raster covariates and align them to the prediction grid
elevation <- raster("data/original/geodata/elevation2000_100m.tif")
slope <- raster("data/original/geodata/slope2000_100m.tif")
nle <- raster("data/original/geodata/uga_viirs_100m_2016.tif")
rainfall <- raster("data/original/geodata/Sum_Rainfall7.tif")
aridity <- raster("data/original/geodata/Annual_Aridity_index.tif")
rain_anomaly <- raster("data/original/geodata/Normalized_Precipitation_Anomaly_2010_.tif")
poverty <- raster("data/original/geodata/uga11povmpi.tif")

ttime <- raster("../rasters/2015_accessibility_to_cities_v1.0.tif") %>% 
  crop(projectRaster(pred$raster, crs = crs(elevation))) 

lst <- stack("../rasters/cru_ts4.04.1901.2019.tmp.dat.nc") %>% 
  crop(projectRaster(pred$raster, crs = crs(elevation)))

id <- which(as.numeric(substr(names(lst), start = 2, stop = 5)) %in%  2011:2016)

lst <- lst[[id]]
lst <- calc(lst, mean)

# Assign CRS to raster with missing CRS
crs(rainfall) <- crs(st_crs(4326)$proj4string)
crs(aridity) <- crs(st_crs(4326)$proj4string)
crs(rain_anomaly) <- crs(st_crs(4326)$proj4string)

elevation <- align_raster(pred = pred$raster, cov = elevation)
slope <- align_raster(pred = pred$raster, cov = slope)
nle <- align_raster(pred = pred$raster, cov = nle)
rainfall <- align_raster(pred = pred$raster, cov = rainfall)
lst <- align_raster(pred = pred$raster, cov = lst)
ttime <- align_raster(pred = pred$raster, cov = ttime)
aridity <- align_raster(pred = pred$raster, cov = aridity)
rain_anomaly <- align_raster(pred = pred$raster, cov = rain_anomaly)
poverty <- align_raster(pred = pred$raster, cov = poverty)

# Stack all the raster covariates together
covariates <- stack(elevation, slope, nle, rainfall, lst, aridity, ttime,
                    pred$pop, rain_anomaly, poverty)
names(covariates) <- c("elevation", "slope", "nle", "rainfall", "lst",
                       "aridity", "ttime", "pop", "rain_anomaly", "poverty")
saveRDS(covariates, "output/covariates/stack_covariates.rds")

# Extract them at the observed locations
cov_obs  <- extract(covariates, df_clusters[, c("utm_x", "utm_y")])

# Fill NA
idna <- unique(which(is.na(cov_obs), arr.ind = T)[,1])
cov_obs[idna, ] <- extract(covariates, df_clusters[idna, c("utm_x", "utm_y")],
                           small = T, buffer = 4, fun = mean)

colnames(cov_obs) <- c("elevation", "slope", "nle", "rainfall", "lst",
                       "aridity", "ttime", "pop", "rain_anomaly", "poverty")

# cov_obs <- cbind(cov_obs, urbanicity = extract(urbanicity, df_clusters[c("utm_x", "utm_y")]))

cov_obs <- cov_obs %>% 
  cbind(df_clusters[c("utm_x", "utm_y")]) %>% 
  as_tibble() 

# cov_obs$urbanicity <- ifelse(cov_obs$urbanicity == 1, "urban", "rural")

df_clean <- df_clean %>% 
  inner_join(cov_obs)

# Assign SPEI EVI and Rainfall to each child based on month and year of birth
# for up to 4 time lags

# SPEI lagged at 1 month 
spei <- stack("data/original/geodata/SPEI1_2010_2018_Monthly.tif")
names(spei) <- paste0(rep(1:12, length(2010:2018)), rep(2010:2018, each = 12))

# EVI
evi <- stack("data/original/geodata/EVI_2010-2018_monthly.tif")
names(evi) <- paste0(rep(1:12, length(2010:2018)), rep(2010:2018, each = 12))

# Rainfall
rain <- readr::read_csv("data/original/RFE2_JAN172020_JAN012001.csv", 
                        col_names = F)

names(rain) <- as.Date(as.Date("2001-01-01"):as.Date("2020-01-17"), 
                       origin = "1970/01/01") %>% 
  as.character()

# Subset to 2010-01-01 to 2018-12-31
rain <- rain[,match("2010-01-01", names(rain)):match("2018-12-31", names(rain))]


rain_coords <- readr::read_csv("data/original/RFE_Rainfall_Coordinates.csv")
rain$long <- rain_coords$X
rain$lat <- rain_coords$Y

# Calculate monthly averages
# Convert to rasters
rain_rast <- rain %>% 
  tidyr::gather(key = "date", value = "rainfall", -long, -lat) %>% 
  mutate(date = as.Date(date), month = lubridate::month(date), 
         year = lubridate::year(date)) %>% 
  group_by(month, year, long, lat) %>% 
  summarise(rainfall = mean(rainfall)) %>% 
  ungroup() %>% 
  mutate(col_name = paste0("X", month, year)) %>% 
  select(-month, -year) %>% 
  tidyr::spread(key = col_name, value = rainfall) %>% 
  rasterFromXYZ(crs = crs(st_crs(4326)$proj4string))

# Fix order of rasters
rain_rast <- rain_rast[[names(spei)]]
writeRaster(rain_rast, "data/processed/geodata/rainfall.tif")

time_vars <- df_clean %>% 
  distinct(bmonth, byear, long, lat)

match_id <- paste0("X", time_vars$bmonth, time_vars$byear)
lags <- -12:-1
for (i in 1:nrow(time_vars)) {
  
  for (j in lags) {
    t_id <- match(match_id[i], names(spei))
    col_name <- paste0("spei", ifelse(j < 0, "m", "p"), abs(j))
    time_vars[i, col_name] <- extract(spei[[t_id + j + 1]], 
                                      time_vars[i, c("long", "lat")])
    
    t_id <- match(match_id[i], names(evi))
    col_name <- paste0("evi", ifelse(j < 0, "m", "p"), abs(j))
    time_vars[i, col_name] <- extract(evi[[t_id + j]], 
                                      time_vars[i, c("long", "lat")])
    
    t_id <- match(match_id[i], names(rain_rast))
    col_name <- paste0("rain", ifelse(j < 0, "m", "p"), abs(j))
    time_vars[i, col_name] <- extract(rain_rast[[t_id + j]], 
                                      time_vars[i, c("long", "lat")])
  }
  
  print(i)

  }

# Merge with main dataset
df_clean <- df_clean %>% 
  left_join(time_vars)

# Extract predictors
spei <- calc(spei[[grep(2016, names(spei))]], mean)
spei <- align_raster(pred = pred$raster, cov = spei)

evi <- calc(evi[[grep(2016, names(evi))]], mean)
evi <- align_raster(pred = pred$raster, cov = evi)

rainfall <- calc(rain_rast[[grep(2016, names(rain_rast))]], mean)
rainfall <- align_raster(pred = pred$raster, cov = rainfall)


predictors <- stack(elevation, slope, nle, rainfall, lst, aridity, ttime,
                    evi, pred$pop, rain_anomaly, poverty, spei)
names(predictors) <- c("elevation", "slope", "nle", "rainfall", "lst",
                       "aridity", "ttime", "evi", "pop", "rain_anomaly", 
                       "poverty", "spei")
saveRDS(predictors, "output/covariates/stack_predictors.rds")

# SAVE DATASETS ----------------------------------------------------------------
readr::write_csv(df_clean, "data/processed/malnutrition.csv")
readr::write_csv(df_clusters, "data/processed/malnutrition_clusters.csv")
