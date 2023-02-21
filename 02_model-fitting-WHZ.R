# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("lme4", "PrevMap", "dplyr", "ggplot2", "sf", "splines") # package names
pacman::p_load(pkgs, character.only = T)

if(all(c("dplyr", "raster") %in% .packages())) {
  select <- dplyr::select
}

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Malnutrition data
df <- readr::read_csv("data/processed/malnutrition.csv")

# Selected covariates (according to HAZ)
covariates <- readRDS("output/covariates/covariates_list.rds")

# EDA --------------------------------------------------------------------------

whz <- df[df$metric == "WHZ", ]

# Remove mean rainfall (we have it lagged)
whz$rainfall <- NULL

# Log-transform elevation
whz$elevation <- log(whz$elevation)

# Standardise covariates
whz$age <- (whz$age - mean(whz$age)) / sd(whz$age)

covs <- covariates

centered_cov <- scale(whz[covs], center = T, scale = T)
mean_cov <- attr(centered_cov, "scaled:center")
sd_cov <- attr(centered_cov, "scaled:scale")
whz[covs] <- centered_cov 

# Correlation plot 
ggcorrplot::ggcorrplot(
  corr = cor(dplyr::select(whz, covs)),
  type = "lower",
  ggtheme = ggplot2::theme_minimal,
  hc.order = T,
  show.diag = F,
  #p.mat = ggcorrplot::cor_pmat(sth[, cov_names]),
  outline.col = "white",
  lab = T,
  legend.title = "Correlation",
  tl.cex = 11,
  tl.srt = 55
) 

ggsave("figs/cov_corr.pdf", width = 8.5, height = 7.5)

# Relationship with covariates
whz %>% 
  dplyr::select(score, covs) %>% 
  tidyr::gather(key = "variable", value = "value", -score) %>%
  ggplot(aes(x = value, y = score)) +
  geom_point(shape = 19, alpha = .05, size = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"),size = .5, aes(col = "GAM", fill = "GAM"),
              method.args = list(method = "REML", gamma = 1.4)) +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), size = .5, fill = "orange", col = "orange") +
  geom_smooth(method = "lm", formula = y ~ x, size = .5, aes(col = "LM", fill = "LM")) +
  scale_color_manual("Method", values = c("blue", "orange")) +
  scale_fill_manual("Method", values = c("blue", "orange")) +
  facet_wrap(~ variable, scales = "free") +
  labs(y = "WHZ", x = "") +
  theme_bw() +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"))

ggsave("figs/cov_rel_WHZ.pdf", width = 8.5, height = 9)

# NON - SPATIAL MODEL ----------------------------------------------------------

# Before fitting this type of models it is alway a good idea to rescale the
# numeric variables (substract the mean and divide by the standard devation)
# to avoid problem of convergence

# Create new ID for clusters
whz <- as.data.frame(whz)
whz$gender <- whz$gender - 1

# Formula
f <- as.formula(paste("score ~", paste(covariates, collapse = "+")))

# Formula with cluster level random effects
fcluster <- update(f, ~ . + (1|ID))

# Fit the model in equation (1)
fit <- lmer(formula = fcluster, data = whz)

# Generate summary of the model
summary(fit)

# Extract the random effects at cluster level U_i
reff <- ranef(fit)$ID$`(Intercept)`

# Extract coordinates for each cluster
coords <- whz %>%
  distinct(ID, utm_x, utm_y) %>% 
  dplyr::select(utm_x, utm_y)

# Compute and plot variogram
library(geoR)
ggvario(coords = coords, data = reff, nsim = 1000, show_nbins = F, maxdist = 50)  
ggsave("figs/variogram_WHZ.pdf", width = 8, height = 6)

# GEOSTATISTICAL MODEL ---------------------------------------------------------
fit_geo <- linear.model.MLE(formula = f,
                            coords = ~ utm_x + utm_y, 
                            ID.coords = create.ID.coords(whz, ~ utm_x + utm_y),
                            data = whz,
                            start.cov.pars = c(10, 1),
                            fixed.rel.nugget = 0,
                            kappa = 0.5, method = "nlminb", 
                            messages = T)
summary(fit_geo, l = F)

# Save the results
saveRDS(fit_geo, file = "output/fit_geo_WHZ.rds")

# PREDICTIONS -----------------------------------------------------------------

# Prediction locations and pop 
pred <- readRDS("data/processed/pred_pop_grid.rds")

# Predictors
predictors <- readRDS("output/covariates/stack_predictors.rds")
rain_lag <- grep("rain", covariates, value = T)
evi_lag <- grep("evi", covariates, value = T)
spei_lag <- grep("spei", covariates, value = T)
names(predictors)[grep("rainfall", names(predictors))] <- rain_lag
names(predictors)[grep("evi", names(predictors))] <- evi_lag
names(predictors)[grep("spei", names(predictors))] <- spei_lag


cov_pred <- extract(predictors, pred$coords) %>% 
  as.data.frame() %>% 
  select(covs)

# Scale predictors
cov_pred <- as.data.frame(scale(cov_pred, 
                                center = mean_cov[names(cov_pred)], 
                                scale = sd_cov[names(cov_pred)]))

# Remove missing covariates
id_rm <- unique(which(is.na(cov_pred), arr.ind = T)[, 1])

cov_pred <- cov_pred[-id_rm, ] 
pred$coords <- pred$coords[-id_rm, ]

# Obtain joint predictions at district level

# Load shapefile for districts in Uganda
districts <- st_read("data/processed/geodata/districts.shp") %>% 
  st_transform(epsgKM(32635))

# Assign district to each prediction point
id_district <- pred$coords %>% 
  as.data.frame() %>% 
  st_as_sf(coords = c(1, 2), crs = epsgKM(32635)) %>% 
  st_join(districts["NAME_1"]) %>% 
  pull("NAME_1") %>% 
  as.character()

# Split the coordinates in a list according to the districts and calculate the
# joint predictions
input <- data.frame(pred$coords, cov_pred, id_district)

# Remove the pred points that fall outside the IUs boundaries
# do this also for the predictors
keep <- !is.na(input$id_district)
input <- input[keep, ]

# Be careful with split because it messes up the order of the chunks of data
# according to the order of f

# Indicate the levels in the factor so split will respect the original order
input_list <- split(input, f = input$id_district)


# Predictions
pred_list <- lapply(input_list,
                    function(x) {
                      
                      out <- spatial.pred.linear.MLE(object = fit_geo,
                                                     grid.pred = x[c("x", "y")],
                                                     predictors = x,
                                                     type = "joint",
                                                     n.sim.prev = 1000)
                      return(out$samples)
                    })


samples <- do.call(rbind, pred_list)
predictions <- list(samples = samples, 
                    coords = do.call(rbind, input_list)[c("x", "y")])

saveRDS(predictions, "output/geo_pred_WHZ.rds")