# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("sf", "raster", "rasterVis", "ggplot2", "PrevMap", "knitr",
         "kableExtra") # package names
pacman::p_load(pkgs, character.only = T)

if(all(c("dplyr", "raster") %in% .packages())) {
  select <- dplyr::select
}

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Shapefile for Uganda
uga <- st_read("data/original/geodata/gadm36_UGA.gpkg", layer = "gadm36_UGA_0")

# SPEI
spei <- stack("data/original/geodata/SPEI1_2010_2018_Monthly.tif")

# EVI
evi <- stack("data/original/geodata/EVI_2010-2018_monthly.tif")

# Rainfall
rain <- stack("data/processed/geodata/rainfall.tif")

# Full malnutrition dataset
df <- readr::read_csv("data/processed/malnutrition.csv")

# Fitted models
fit_haz <- readRDS("output/fit_geo_HAZ.rds")
fit_waz <- readRDS("output/fit_geo_WAZ.rds")
fit_whz <- readRDS("output/fit_geo_WHZ.rds")

# Rasters of covariates
covariates <- readRDS("output/covariates/stack_covariates.rds")

# DATA PREPARATION -------------------------------------------------------------

# Change names
months <- as.character(lubridate::month(1:12, label = T))
names_months <- paste(rep(months, length(2010:2018)), rep(2010:2018, each = 12))
names(spei) <- names(evi) <- names(rain) <- names_months

# Select only from 2011 to 2016
ids <- match(c("Jan 2011", "Dec 2016"), names_months)
spei <- spei[[ids[1]:ids[2]]]
evi <- evi[[ids[1]:ids[2]]]
rain <- rain[[ids[1]:ids[2]]]


# Crop to Uganda
spei <- spei %>% 
  crop(uga) %>% 
  mask(uga)


evi <- evi %>% 
  crop(uga) %>% 
  mask(uga)

rain <- rain %>% 
  crop(uga) %>% 
  mask(uga)

# Convert Uganda boundaries CRS to utm
uga_utm <- uga %>% 
  st_transform(epsgKM(32635))

# MAPS OF RASTER COVARIATES ----------------------------------------------------

# Time varying

# SPEI
ggspei <- gplot(spei) +
  geom_tile(aes(fill = value)) +
  geom_polygon(aes(x = long, y = lat, group = group),
               data = as(uga, "Spatial"), color = "black",
               fill = NA, size = .3) +
  facet_wrap(~ variable, 
             labeller = labeller(variable = function(x) gsub("\\.", " ", x)), 
             ncol = 12) +
  scale_fill_viridis_c(na.value = NA, direction = -1) +
  coord_equal() +
  theme_bw(base_size = 13) +
  labs(x = "", y = "", fill = "SPEI") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  theme(line = element_blank(), axis.text = element_blank(), 
        legend.position = "top", legend.title = element_text(face = "bold")) 
ggsave(plot = ggspei, filename = "figs/spei_ts.png", width = 9.5, height = 8)
ggsave(plot = ggspei, filename = "figs/spei_ts.pdf", width = 9.5, height = 8)

# EVI
ggevi <- gplot(evi) +
  geom_tile(aes(fill = value)) +
  geom_polygon(aes(x = long, y = lat, group = group),
               data = as(uga, "Spatial"), color = "black",
               fill = NA, size = .3) +
  facet_wrap(~ variable, 
             labeller = labeller(variable = function(x) gsub("\\.", " ", x)), 
             ncol = 12) +
  scale_fill_gradient(na.value = NA, low = "saddlebrown", high = "green3") +
  coord_equal() +
  theme_bw(base_size = 13) +
  labs(x = "", y = "", fill = "EVI") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  theme(line = element_blank(), axis.text = element_blank(), 
        legend.position = "top", legend.title = element_text(face = "bold"),
        legend.key.width = unit(x = 2, units = "cm"))  

ggsave(plot = ggevi, filename = "figs/evi_ts.png", width = 9.5, height = 8)
ggsave(plot = ggevi, filename = "figs/evi_ts.pdf", width = 9.5, height = 8)

# Rain
ggrain <- gplot(rain) +
  geom_tile(aes(fill = value)) +
  geom_polygon(aes(x = long, y = lat, group = group),
               data = as(uga, "Spatial"), color = "black",
               fill = NA, size = .3) +
  facet_wrap(~ variable, 
             labeller = labeller(variable = function(x) gsub("\\.", " ", x)), 
             ncol = 12) +
  scale_fill_viridis_c(na.value = NA) +
  coord_equal() +
  theme_bw(base_size = 13) +
  labs(x = "", y = "", fill = "Rainfall") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  theme(line = element_blank(), axis.text = element_blank(), 
        legend.position = "top", legend.title = element_text(face = "bold"))  

ggsave(plot = ggrain, filename = "figs/rain_ts.png", width = 9.5, height = 8)
ggsave(plot = ggrain, filename = "figs/rain_ts.pdf", width = 9.5, height = 8)



# Time-constant
covs_in_model <- c("slope", "lst", "aridity", "ttime", "nle")
covariates <- covariates[[covs_in_model]]
cov_names <- c("Slope angle", "Mean temperature", "Aridity index",
               "Travel time cities", "NLE  (log-scale)")
covariates[["nle"]] <- log(covariates[["nle"]])
gglist <- list()

for (i in 1:length(covs_in_model)) {
  gg <- gplot(covariates[[i]]) +
    geom_tile(aes(fill = value)) +
    geom_polygon(aes(x = long, y = lat, group = group),
                 data = as(uga_utm, "Spatial"), color = "black",
                 fill = NA, size = .8) +
    scale_fill_viridis_c(na.value = NA) +
    coord_equal() +
    theme_bw(base_size = 8) +
    labs(x = "", y = "", fill = cov_names[[i]]) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
    theme(line = element_blank(), axis.text = element_blank(), 
          legend.position = "top", legend.title = element_text(face = "bold"), 
          rect = element_blank()) 
  ggsave(plot = gg, filename = paste0("figs/", names(covariates[[i]]), ".png"),
         width = 5, height = 5)
  ggsave(plot = gg, filename = paste0("figs/", names(covariates[[i]]), ".pdf"),
         width = 5, height = 5)
  gglist[[i]] <- gg
}
 
cowplot::plot_grid(plotlist = gglist)
ggsave(filename = "figs/covs_static.png", width = 6)
ggsave(filename = "figs/covs_static.pdf", width = 6)

# RELATIONSHIP WITH AGE and GENDER ---------------------------------------------

agerel <- ggplot(df, aes(x = age, y = score, col = metric, fill = metric)) +
  geom_point(alpha = .2, shape = 19, size = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"),
              method.args = list(method = "REML", gamma = 1.4)) +
  facet_wrap(~ metric, ncol = 3) +
  labs(x = "Age (months)", y = "Z-score") +
  scale_x_continuous(breaks = seq(0, 60, by = 10)) +
  scale_y_continuous(breaks = seq(-6, 6, by = 1)) +
  scale_color_brewer(type = "q", palette = 7) +
  scale_fill_brewer(type = "q", palette = 7) +
  theme_bw(16) +
  theme(legend.position = "none")

ggsave(plot = agerel, "figs/age_rel.pdf", width = 9, height = 6)

genderel <- ggplot(df, aes(x = as.factor(gender), y = score, fill = metric)) +
  geom_boxplot(aes(col = metric), alpha = .2) +
  facet_wrap(~ metric, ncol = 3) +
  labs(y = "Z-score", x = "") +
  scale_x_discrete(labels = c("Male", "Female")) +
  scale_y_continuous(breaks = seq(-6, 6, by = 1)) +
  scale_color_brewer(type = "q", palette = 7) +
  scale_fill_brewer(type = "q", palette = 7) +
  theme_bw(16) +
  theme(legend.position = "none")

ggpubr::ggarrange(agerel, genderel, ncol = 1, labels = c("A", "B"))
ggsave("figs/age_gender_rel.pdf", width = 9, height = 12)

# TABLES -----------------------------------------------------------------------

haz_tab <- create_tab(fit_haz)
waz_tab <- create_tab(fit_waz)
whz_tab <- create_tab(fit_whz)

all_tab <- haz_tab %>% 
  full_join(waz_tab, by = "Parameter") %>% 
  full_join(whz_tab, by = "Parameter") 

par_latex <- c(1, (nrow(all_tab)-2):nrow(all_tab))
all_cov <- all_tab[-par_latex, ]

all_nocov <- all_tab[par_latex, ]
all_nocov$Parameter <- c("$\\alpha$",
                         "$\\sigma^2$",
                         "$\\phi$",
                         "$\\omega^2$")

all_tab <- bind_rows(all_nocov[1, ],
                     all_cov, 
                     all_nocov[2:4, ])

all_tab %>% 
  readr::write_csv("output/summary_table.csv")

tab <- kable(all_tab, booktabs = T, linesep = "", align = "l", escape = F,
             caption = "Monte Carlo maximum likelihood estimates and corresponding 95\\% confidence intervals.",
             format = "latex",
             col.names = c("Parameter", rep(c("Estimate", "95\\% CI"), times = 3))) %>%
  kable_styling(latex_options = c("HOLD_position", "scale_down")) %>%
  add_header_above(c(" " = 1, "HAZ" = 2, "WAZ" = 2,
                     "WHZ" = 2), bold = T) %>%
  row_spec(0, bold = T)

saveRDS(tab, "output/model_table.rds")
