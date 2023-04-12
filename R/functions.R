# Function to create grid of prediction and 
# project the population raster
create_grid <- function(resolution, study_area, pop, cutoff = 0, filename = NULL) {
  crs <- st_crs(study_area)$proj4string
  pred_grid <- raster(res = resolution, ext = extent(study_area), crs = crs, val = 1)
  res_pred <- res(pred_grid)[1]
  res_pop <- res(pop)[1] * 110.567
  fct <- floor(res_pred / res_pop)
  if(fct >= 2) pop <- aggregate(pop, fun = sum, fact = fct) 
  if(crs != projection(pop)) pop <- projectRaster(pop, pred_grid)
  pop[pop < 0] <- 0
  compare_cond <- compareRaster(pred_grid, pop, res = T, orig = T, 
                                stopiffalse = F, showwarning = T)
  stopifnot(compare_cond)    
  pop_mask <- pop
  pop_mask[pop_mask <= cutoff] <- NA
  # writeRaster(pop, filename)
  pred_raster = mask(pred_grid, pop_mask)
  out <- list(pop = pop, 
              raster = pred_raster,
              coords = rasterToPoints(pred_raster)[, 1:2])
  if(!is.null(filename)) saveRDS(out, file = filename)
  return(out) 
}

align_raster <- function(pred, cov, filename = NULL) {
  crs_pred <- projection(pred)
  crs_cov <- projection(cov)
  res_pred <- res(pred)[1]
  if(length(grep("\\+proj=longlat", crs_cov)) == 1) multiplier <- 110.567
  if(length(grep("\\+units=km", crs_cov)) == 1) multiplier <- 1
  if(length(grep("\\+units=m", crs_cov)) == 1) multiplier <- 1 / 1000
  res_cov <- min(res(cov)) * multiplier
  fct <- floor(res_pred / res_cov)
  if(fct >= 2) cov <- aggregate(cov, fun = mean, fact = fct) 
  if(crs_pred != crs_cov) cov <- projectRaster(cov, pred)
  
  compare_cond <- compareRaster(pred, cov, res = T, orig = T, 
                                stopiffalse = F, showwarning = T)
  stopifnot(compare_cond)    
  cov <- mask(cov, pred)
  cov_name <- sub("_.*", "", x = names(cov))
  names(cov) <- paste0(cov_name, "_", res_pred, "km")
  if(!is.null(filename)) writeRaster(cov, filename = filename)
  cov
}


plotRaster <- function(x, ...) rasterVis::levelplot(x, margin = F, ...)

# Convert epsg to epsg KM
epsgKM <- function(x) {
  crs <- st_crs(x)
  proj4KM <- gsub(pattern = "+.units=m", replacement = "+units=km", 
                  crs$proj4string)
  return(proj4KM)
}

elogit <- function(num, den) {
  log((num + 0.5) / (den - num + 0.5))
}

# Envelope for variogram
variog_envelope <- function (geodata, coords = geodata$coords, data = geodata$data, 
                              obj.variog, nsim = 99, save.sim = FALSE, messages) 
{
  call.fc <- match.call()
  if (missing(geodata)) 
    geodata <- list(coords = coords, data = data)
  if (missing(messages)) 
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                         TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  obj.variog$v <- NULL
  if ((is.matrix(data) | is.data.frame(data))) 
    if (ncol(data) > 1) 
      stop("envelops can be computed for only one data set at once")
  if (!is.null(obj.variog$estimator.type)) 
    estimator.type <- obj.variog$estimator.type
  else estimator.type <- "classical"
  if (abs(obj.variog$lambda - 1) > 1e-04) {
    if (abs(obj.variog$lambda) < 1e-04) 
      data <- log(data)
    else data <- ((data^obj.variog$lambda) - 1)/obj.variog$lambda
  }
  xmat <- unclass(trend.spatial(trend = obj.variog$trend, geodata = geodata))
  if (obj.variog$trend != "cte") {
    if (is.vector(data)) {
      data <- lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    }
    else {
      only.res <- function(y, x) {
        lm(y ~ xmat + 0)$residuals
      }
      data <- apply(data, 2, only.res, x = xmat)
    }
  }
  if (messages.screen) 
    cat(paste("variog.env: generating", nsim, "simulations by permutating data values\n"))
  simula <- list(coords = coords)
  n.data <- length(data)
  perm.f <- function(i, data, n.data) {
    return(data[sample(1:n.data)])
  }
  simula$data <- apply(as.matrix(1:nsim), 1, perm.f, data = data, 
                       n.data = n.data)
  if (messages.screen) 
    cat(paste("variog.env: computing the empirical variogram for the", 
              nsim, "simulations\n"))
  nbins <- length(obj.variog$bins.lim) - 1
  if (obj.variog$direction == "omnidirectional") {
    bin.f <- function(sim) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      temp <- .C("binit", as.integer(obj.variog$n.data), 
                 as.double(as.vector(coords[, 1])), as.double(as.vector(coords[, 
                                                                               2])), as.double(as.vector(sim)), as.integer(nbins), 
                 as.double(as.vector(obj.variog$bins.lim)), as.integer(estimator.type == 
                                                                         "modulus"), as.double(max(obj.variog$u)), as.double(cbin), 
                 vbin = as.double(vbin), as.integer(FALSE), as.double(sdbin), 
                 PACKAGE = "geoR")$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f)
  }
  else {
    variog.vbin <- function(x, ...) {
      variog(geodata = geodata, 
             data = x, uvec = obj.variog$uvec, estimator.type = obj.variog$estimator.type, 
             nugget.tolerance = obj.variog$nugget.tolerance, max.dist = obj.variog$max.dist, 
             pairs.min = obj.variog$pairs.min, direction = obj.variog$direction, 
             tolerance = obj.variog$tolerance, messages.screen = FALSE,...)$v
    }
    simula.bins <- apply(simula$data, 2, variog.vbin)
  }
  simula.bins <- simula.bins[obj.variog$ind.bin, ]
  if (save.sim == FALSE) 
    simula$data <- NULL
  if (messages.screen) 
    cat("variog.env: computing the envelops\n")
  limits <- apply(simula.bins, 1, quantile, prob = c(0.025, 0.975))
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ], 
                  v.upper = limits[2, ])
  if (save.sim) 
    res.env$simulations <- simula$data
  res.env$call <- call.fc
  oldClass(res.env) <- "variogram.envelope"
  return(res.env)
}

# Calculate and plot the variogram
ggvario <- function(coords, 
                    data, 
                    bins = 15, 
                    maxdist = max(dist(coords))/3, 
                    uvec = NULL, 
                    nsim = 999,
                    color = "royalblue1", 
                    xlab = "distance", 
                    show_nbins = T) {
  coords <- as.matrix(coords)
  min_dist <- min(dist(coords))
  if(is.null(uvec)) uvec <- seq(min_dist, maxdist, l = bins)
  empvario <- variog(coords = coords, data = data, uvec = uvec, messages = F)
  envmc <- variog_envelope(coords = coords, data = data, 
                           obj.variog = empvario, nsim = nsim, messages = F)
  dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                        lowemp = envmc$v.lower, upemp = envmc$v.upper, 
                        nbins = empvario$n)
  p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
    geom_ribbon(aes(ymin = lowemp, ymax = upemp), fill = color, alpha = .3) +
    geom_point(aes(y = empirical), col = "black", fill = color, shape = 21, size = 3) +
    scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                       breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
    scale_y_continuous(name = "semivariance", 
                       #breaks = round(seq(0, max(dfvario$upemp, dfvario$empirical), l = 6), 1), 
                       limits = c(0, max(dfvario$upemp, dfvario$empirical))) +
    ggtitle("Empirical semivariogram") +
    theme_classic()
  p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
  if(show_nbins) p2 else p1
}

create_labels <- function(x, greater = F, smaller = F, text_sep = " - ") {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = text_sep)
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }
  
  return(labs)
}

floor_dec <- function(x, digits=1) round(x - 5*10^(-digits-1), digits)
ceiling_dec <- function(x, digits=1) round(x + 5*10^(-digits-1), digits)

raster_map <- function(x, raster, shape, boundaries, palette, breaks = NULL, n, 
                       digits = 0, legend_title, panel_title = NULL, labg = F, 
                       text_sep = " - ", labels = NULL, ...) {
  
  # Generate breaks for color legend
  if (is.null(breaks)) {
    breaks <- as.numeric(quantile(raster[[x]], probs = seq(0, 1, l = n + 1)))
    breaks[2:n] <- round(breaks[2:n], digits)
    breaks[1] <- floor_dec(breaks[1], digits)
    breaks[n + 1] <- ceiling_dec(breaks[n + 1], digits) 
  }
  
  # Generate labels for color legend
  if(is.null(labels)) {
    labs <- create_labels(breaks, greater = labg, text_sep = text_sep)
  }
  
  
  # Generate color palette
  if (is.character(palette)) {
    pal <- tmaptools::get_brewer_pal(palette, n = length(labs), contrast = c(0, 1), plot = F)
  } else {
    pal <- palette(length(labs))
  }
  
  # Produce map
  tm_shape(boundaries) +
    tm_polygons(col = "grey", border.col = "black") + 
    tm_shape(raster) +
    tm_raster(col = x, palette = pal, breaks = breaks, style = "fixed", 
              legend.show = T, labels = labs, title = legend_title,
              midpoint = "NA") +
    tm_shape(shape, is.master = T) +
    tm_borders(col = "black") +
    tm_compass(position = c("right", "top")) +
    tm_scale_bar(position = c("right", "bottom"), text.size = .8) +
    tm_layout(outer.margins = 0, asp = 0,
              # =legend.bg.color = "white", legend.frame = "black",
              legend.text.size = 0.95, legend.title.size = 1.35,
              panel.show = T, panel.labels = panel_title, 
              panel.label.size = 1.4,
              # panel.label.fontfamily = "bold", 
              ...)
}


district_map <- function(x, shape, palette, breaks = NULL, n, digits = 0, 
                         legend_title, panel_title = NULL, labg = F, ...) {
  
  # Generate breaks for color legend
  if (is.null(breaks)) {
    breaks <- as.numeric(quantile(shape[[x]], probs = seq(0, 1, l = n + 1), 
                                  na.rm = T))
    breaks[2:n] <- round(breaks[2:n], digits)
    breaks[1] <- floor_dec(breaks[1], digits)
    breaks[n + 1] <- ceiling_dec(breaks[n + 1], digits) 
  }
  
  # Generate labels for color legend
  labs <- create_labels(breaks, greater = labg)
  
  # Generate color palette
  if (is.character(palette)) {
    pal <- tmaptools::get_brewer_pal(palette, n = length(labs), contrast = c(0, 1), plot = F)
  } else {
    pal <- palette(length(labs))
  }
  
  # Produce map
  tm_shape(shape) +
    tm_polygons(col = x, palette = pal, breaks = breaks, style = "fixed", 
                legend.show = T, labels = labs, title = legend_title,
                border.col = "black", showNA = F) +
    tm_compass(position = c("right", "top")) +
    tm_scale_bar(position = c("right", "bottom"), text.size = .8) +
    tm_layout(outer.margins = 0, asp = 0,
              # =legend.bg.color = "white", legend.frame = "black",
              legend.text.size = 0.95, legend.title.size = 1.35,
              panel.show = T, panel.labels = panel_title, 
              panel.label.size = 1.4,
              # panel.label.fontfamily = "bold", 
              ...)
}

create_tab <- function(x) {
  require(dplyr)
  
  # Parameters 
  pars <- as.numeric(x$estimate)
  
  # Standard errors
  se <- sqrt(diag(x$covariance))
  
  # Parameters on their support
  estimates <- PrevMap:::coef.PrevMap(x) 
  p <- length(estimates) - 3
  estimates <- round(estimates, 5)
  
  # CIs
  ci_up <- pars + qnorm(0.975) * se
  ci_low <- pars - qnorm(0.975) * se
  ci <- cbind(ci_low, ci_up) 
  ci[p + 3, ] <- ci[p + 3, ] + ci[p + 1, ]
  ci[(p + 1):length(estimates), ] <- exp(ci[(p + 1):length(estimates), ])
  ci <- round(ci, 5)
  ci <- apply(ci, 1, function(x) paste0("(", x[1], ", ", x[2], ")"))
  
  # P-values
  t <- abs(pars / se)
  p_values <- 2 * (1 - pnorm(t))
  round(p_values, digits = 6)
  p_values[c(1, (p+1):(p+3))] <- NA
  
  # Put everything together
  tab <- tibble(Parameter = names(estimates), Estimate = estimates, CI = ci,
                pvalue = p_values)
  return(tab)
}
