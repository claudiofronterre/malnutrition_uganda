# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("dplyr", "sf", "ggplot2", "raster", "tmap") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Uganda border
uga <- rgeoboundaries::geoboundaries("uganda") %>% 
  st_transform(epsgKM(32635))

# Districts
uga_districts <- st_read("data/original/geodata/gadm36_UGA.gpkg", 
                         layer = "gadm36_UGA_1") %>% 
  st_transform(epsgKM(32635))

haz_rast <- readRDS("output/mapping/rast_HAZ.rds")
waz_rast <- readRDS("output/mapping/rast_WAZ.rds")
whz_rast <- readRDS("output/mapping/rast_WHZ.rds")

haz_burden <- readRDS("output/mapping/district_HAZ.rds")
waz_burden <- readRDS("output/mapping/district_WAZ.rds")
whz_burden <- readRDS("output/mapping/district_WHZ.rds")

# HAZ --------------------------------------------------------------------------

haz_df <- haz_rast %>% 
  rasterToPoints() %>% 
  as_tibble() %>% 
  tidyr::gather(key = "summary", value = "value", -x, -y)

# Mean Score
df <- haz_df[haz_df$summary == "mean_score", ]
lims <- c(round(min(df$value), 1), round(max(df$value), 1))
brks <- seq(lims[1], lims[2], l = 5)
g1 <- ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = value)) +
  geom_sf(data = uga_districts, inherit.aes = F, fill = NA, col = "white") +
  geom_sf(data = uga, inherit.aes = F, fill = NA, col = "black") +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_fill_distiller("HAZ-scores", type = "seq", direction = 1, palette = 3,
                       breaks = brks, limits = lims) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", legend.key.width = unit(1.5, "cm"), 
        legend.title = element_text(face = "bold")) 
g1
ggsave("figs/haz_score.pdf", width = 4, height = 4.5)


# Mean prevalence
df <- haz_df[haz_df$summary == "mean_prev", ]
lims <- c(round(min(df$value)), round(max(df$value))) 
brks <- round(seq(lims[1], lims[2], l = 6)) 
g2 <- ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = value/100)) +
  geom_sf(data = uga_districts, inherit.aes = F, fill = NA, col = "white") +
  geom_sf(data = uga, inherit.aes = F, fill = NA, col = "black") +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_fill_distiller("Stunting Prevalence", type = "seq", direction = 1,
                       breaks = brks / 100, limits = lims / 100,
                       labels = scales::percent, palette = 3) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", legend.key.width = unit(1.5 , "cm"), 
        legend.title = element_text(face = "bold")) 
g2
ggsave("figs/haz_prev.pdf", width = 4, height = 4.5)


# Mean burden
library(ggspatial)
df <- haz_burden
df$burden <- df$burden / 1000
lims <- c(floor(min(df$burden, na.rm = T)), 
          ceiling(max(df$burden, na.rm = T)))  
brks <- round(seq(lims[1], lims[2], l = 6)) 
g3 <- ggplot(df) +
  geom_sf(aes(fill = burden), col = "white") +
  geom_sf(data = uga, inherit.aes = F, fill = NA, col = "black") +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_fill_distiller("Stunting Burden (thousands)", type = "seq", direction = 1,
                       breaks = brks, limits = lims, na.value = NA,
                       palette = 3) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", legend.key.width = unit(1.5 , "cm"), 
        legend.title = element_text(face = "bold")) +
  annotation_scale(location = "br", pad_y = unit(.5, "cm")) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         height = unit(.7, "cm"), width = unit(.7, "cm"), 
                         pad_y =  unit(1.2, "cm"))
g3
ggsave("figs/haz_burden.pdf", width = 4, height = 4.5)

# All together
gall_haz <- cowplot::plot_grid(plotlist = list(g1, g2, g3), nrow = 1)
gall_haz
ggsave("figs/haz_all.pdf", width = 4 * 3, height = 4.5)

# WAZ --------------------------------------------------------------------------

waz_df <- waz_rast %>% 
  rasterToPoints() %>% 
  as_tibble() %>% 
  tidyr::gather(key = "summary", value = "value", -x, -y)

# Mean Score
df <- waz_df[waz_df$summary == "mean_score", ]
lims <- c(round(min(df$value), 1), round(max(df$value), 1))
brks <- seq(lims[1], lims[2], l = 5)
g1 <- ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = value)) +
  geom_sf(data = uga_districts, inherit.aes = F, fill = NA, col = "white") +
  geom_sf(data = uga, inherit.aes = F, fill = NA, col = "black") +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_fill_distiller("WAZ-scores", type = "seq", palette = 16, direction = 1, 
                       breaks = brks, limits = lims) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", legend.key.width = unit(1.5, "cm"), 
        legend.title = element_text(face = "bold")) 
g1
ggsave("figs/waz_score.pdf", width = 4, height = 4.5)


# Mean prevalence
df <- waz_df[waz_df$summary == "mean_prev", ]
lims <- c(round(min(df$value)), round(max(df$value))) 
brks <- round(seq(lims[1], lims[2], l = 6)) 
g2 <- ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = value/100)) +
  geom_sf(data = uga_districts, inherit.aes = F, fill = NA, col = "white") +
  geom_sf(data = uga, inherit.aes = F, fill = NA, col = "black") +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_fill_distiller("Underweight Prevalence", type = "seq", direction = 1,
                       breaks = brks / 100, limits = lims / 100,
                       labels = scales::percent, palette = 16) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", legend.key.width = unit(1.5 , "cm"), 
        legend.title = element_text(face = "bold")) 
g2
ggsave("figs/waz_prev.pdf", width = 4, height = 4.5)


# Mean burden
library(ggspatial)
df <- waz_burden
df$burden <- df$burden / 1000
lims <- c(floor(min(df$burden, na.rm = T)), 
          ceiling(max(df$burden, na.rm = T)))  
brks <- round(seq(lims[1], lims[2], l = 6)) 
g3 <- ggplot(df) +
  geom_sf(aes(fill = burden), col = "white") +
  geom_sf(data = uga, inherit.aes = F, fill = NA, col = "black") +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_fill_distiller("Underweight Burden (thousands)", type = "seq", direction = 1,
                       breaks = brks, limits = lims, na.value = NA,
                       palette = 16) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", legend.key.width = unit(1.5 , "cm"), 
        legend.title = element_text(face = "bold")) +
  annotation_scale(location = "br", pad_y = unit(.5, "cm")) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         height = unit(.7, "cm"), width = unit(.7, "cm"), 
                         pad_y =  unit(1.2, "cm"))
g3
ggsave("figs/waz_burden.pdf", width = 4, height = 4.5)

# All together
gall_waz <- cowplot::plot_grid(plotlist = list(g1, g2, g3), nrow = 1)
gall_waz
ggsave("figs/waz_all.pdf", width = 4 * 3, height = 4.5)

# WHZ --------------------------------------------------------------------------

whz_df <- whz_rast %>% 
  rasterToPoints() %>% 
  as_tibble() %>% 
  tidyr::gather(key = "summary", value = "value", -x, -y)

# Mean Score
df <- whz_df[whz_df$summary == "mean_score", ]
lims <- c(round(min(df$value), 1), round(max(df$value), 1))
brks <- seq(lims[1], lims[2], l = 5)
g1 <- ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = value)) +
  geom_sf(data = uga_districts, inherit.aes = F, fill = NA, col = "white") +
  geom_sf(data = uga, inherit.aes = F, fill = NA, col = "black") +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_fill_distiller("WHZ-scores", type = "seq", palette = 8, direction = 1, 
                       breaks = brks, limits = lims) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", legend.key.width = unit(1.5, "cm"), 
        legend.title = element_text(face = "bold")) 
g1
ggsave("figs/whz_score.pdf", width = 4, height = 4.5)


# Mean prevalence
df <- whz_df[whz_df$summary == "mean_prev", ]
lims <- c(round(min(df$value)), round(max(df$value))) 
brks <- round(seq(lims[1], lims[2], l = 6)) 
g2 <- ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = value/100)) +
  geom_sf(data = uga_districts, inherit.aes = F, fill = NA, col = "white") +
  geom_sf(data = uga, inherit.aes = F, fill = NA, col = "black") +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_fill_distiller("Wasting Prevalence", type = "seq", direction = 1,
                       breaks = brks / 100, limits = lims / 100,
                       labels = scales::percent, palette = 8) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", legend.key.width = unit(1.5 , "cm"), 
        legend.title = element_text(face = "bold")) 
g2
ggsave("figs/whz_prev.pdf", width = 4, height = 4.5)


# Mean burden
library(ggspatial)
df <- whz_burden
df$burden <- df$burden / 1000
lims <- c(floor(min(df$burden, na.rm = T)), 
          ceiling(max(df$burden, na.rm = T)))  
brks <- round(seq(lims[1], lims[2], l = 6)) 
g3 <- ggplot(df) +
  geom_sf(aes(fill = burden), col = "white") +
  geom_sf(data = uga, inherit.aes = F, fill = NA, col = "black") +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_fill_distiller("Underweight Burden (thousands)", type = "seq", direction = 1,
                       breaks = brks, limits = lims, na.value = NA,
                       palette = 8) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", legend.key.width = unit(1.5 , "cm"), 
        legend.title = element_text(face = "bold")) +
  annotation_scale(location = "br", pad_y = unit(.5, "cm")) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         height = unit(.7, "cm"), width = unit(.7, "cm"), 
                         pad_y =  unit(1.2, "cm"))
g3
ggsave("figs/whz_burden.pdf", width = 4, height = 4.5)

# All together
gall_whz <- cowplot::plot_grid(plotlist = list(g1, g2, g3), nrow = 1)
gall_whz
ggsave("figs/whz_all.pdf", width = 4 * 3, height = 4.5)

# Final map
gall <- cowplot::plot_grid(plotlist = list(gall_haz, gall_waz, gall_whz), nrow = 3)
gall
ggsave("figs/malnutrition_all.pdf", width = 4 * 3, height = 4.5 * 3)

# SCORES -----------------------------------------------------------------------

uga <- st_union(uga_districts)

scores_rast <- stack(haz_rast$mean_score, 
                     waz_rast$mean_score, 
                     whz_rast$mean_score)
names(scores_rast) <- c("HAZ", "WAZ", "WHZ")

scores_df <- scores_rast %>% 
  rasterToPoints() %>% 
  as_tibble() %>% 
  tidyr::gather(key = "indicator", value = "value", -x, -y)

# Mean Score
df <- scores_df
brks <- tapply(df$value, df$indicator, FUN = quantile, na.rm = T,
               prob = seq(0, 1, l = 7)) %>% 
  lapply(round, digit = 1)
brks <- seq(-2.5, 1.5, by = 0.3)

pal <- tmaptools::get_brewer_pal("RdYlBu", n =length(brks) - 1, 
                                 contrast = c(0, 1))
crs(scores_rast) <- crs(uga)
tm_shape(uga) +
  tm_fill() +
tm_shape(scores_rast) +
  tm_raster(title = "", palette = pal, style = "fixed",
            breaks = brks, legend.is.portrait = F) +
tm_shape(uga_districts) +
  tm_borders() +
tm_shape(uga) +
  tm_borders(col = "black") +
tm_layout(outer.margins = 0, asp = 0, 
          legend.outside.position = "top", legend.outside.size = .20,
          legend.title.fontface = "bold", legend.position = c("center", .6),
          title = "Predicted mean Z-scores", 
          title.snap.to.legend = T, title.size = .95,
          title.position = c("center", "TOP"), title.fontface = "bold",
          legend.outside = T, legend.just = c("center", "top"))
tmap_save(filename = "figs/zscores.pdf", width = 9, heigh = 4)
# PREVALENCES ------------------------------------------------------------------

prev_rast <- stack(haz_rast$mean_prev, 
                     waz_rast$mean_prev, 
                     whz_rast$mean_prev)
names(prev_rast) <- c("Stunting", "Underweight", "Wasting")

prev_df <- prev_rast %>% 
  rasterToPoints() %>% 
  as_tibble() %>% 
  tidyr::gather(key = "indicator", value = "value", -x, -y)

# Mean prevalence
df <- prev_df
brks <- tapply(df$value, df$indicator, FUN = quantile, na.rm = T,
               prob = seq(0, 1, l = 7)) 
brks <- c(0, 2, 5, 7, 12, 15, 25, 28, 31, 34, 58)

pal <- tmaptools::get_brewer_pal("-RdYlBu", n =length(brks) - 1, 
                                 contrast = c(0, 1))
crs(prev_rast) <- crs(uga)
tm_shape(uga) +
  tm_fill() +
tm_shape(prev_rast) +
  tm_raster(title = "", palette = pal, style = "fixed",
            breaks = brks, legend.is.portrait = F) +
tm_shape(uga_districts) +
  tm_borders() +
tm_shape(uga) +
  tm_borders(col = "black") +
tm_layout(outer.margins = 0, asp = 0, 
          legend.outside.position = "top", legend.outside.size = .20,
          legend.title.fontface = "bold", legend.position = c("center", .6),
          title = "Predicted mean prevalence (%)", 
          title.snap.to.legend = T, title.size = .95,
          title.position = c("center", "TOP"), title.fontface = "bold",
          legend.outside = T, legend.just = c("center", "top"))
tmap_save(filename = "figs/prevalence.pdf", width = 9, heigh = 4)

# BURDEN -----------------------------------------------------------------------

haz_burden$indicator <- "Stunting"
waz_burden$indicator <- "Underweight"
whz_burden$indicator <- "Wasting"

burden <-rbind(haz_burden, waz_burden, whz_burden)
burden$burden <- burden$burden / 1000
 
brks <- round(quantile(burden$burden, na.rm = T, probs = seq(0, 1, l = 11)))
brks <- c(0, 2, 4, 6, 8, 10, 15, 20, 25, 45, 65, 85, 105)
pal <- tmaptools::get_brewer_pal("Reds", n =length(brks) - 1, 
                                 contrast = c(0, 1))
tm_shape(uga) +
  tm_fill() +
tm_shape(burden) +
  tm_polygons(col = "burden", legend.is.portrait = F,
              title = "", pal = pal, style = "fixed", breaks = brks,
              showNA = F) +
tm_facets(by = "indicator") +
tm_shape(uga) +
  tm_borders(col = "black") +
tm_layout(outer.margins = 0, asp = 0, 
          legend.outside.position = "top", legend.outside.size = .20,
          legend.title.fontface = "bold", legend.position = c("center", .6),
          title = "Predicted mean burden (thousands)", 
          title.snap.to.legend = T, title.size = .95,
          title.position = c("center", "TOP"), title.fontface = "bold",
          legend.outside = T, legend.just = c("center", "top"))
tmap_save(filename = "figs/burden.pdf", width = 9, heigh = 4)
