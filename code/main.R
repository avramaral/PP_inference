##############################
# Load packages
##############################

set.seed(123)

library("spatstat") # Point process
library("ggplot2")  # Plotting 
library("tmap")     # Plotting
library("raster")   # Create `raster` object
library("rgeos")    # Plotting
library("INLA")     # Fit Gaussian latent models
library("sf")       # Manipulate spatial objects

##############################

# E.g., simulate a point pattern from an intensity function

f <- function (x, y) { (x ** 2 + y ** 2) }
x <- seq(-1, 1, 0.05)
y <- seq(-1, 1, 0.05)
z <- outer(X = x, Y = y, FUN = f)
w <- owin(xrange = c(-1, 1), yrange = c(-1, 1)) 

pp <- rpoint(n = 200, f = f, win = w)

persp(x, y, z, theta = 30)
plot(pp, main = "")

##############################
# Poisson process
##############################

# Simulating homogeneous Poisson process
sim.HPP <- function (lambda, max = 1, min = 0, ...) {
  m <- lambda * (max - min)
  N <- rpois(1, lambda = m)
  sort(runif(N, min = min, max = max))
}

hpp <- sim.HPP(lambda = 0.5, max = 100, min = 0)
print(hpp)
print(length(hpp))

# Plot 1D PP
plot1D.pp <- function (pp, xlim = NA, ...) {
  if (is.na(xlim)) { xlim <- c(floor(min(pp)), ceiling(max(pp))) }
  plot(x = NA, xlab = "", ylab = "", xlim = xlim, ylim = c(0, 1), axes = FALSE, frame = FALSE, asp = diff(xlim) * 0.25) 
  arrows(x0 = pp, y0 = 1, x1 = pp, y1 = 0)
  axis(1, pos = 0)
}

plot1D.pp(hpp)

# As we will see next, the MLE for lambda is computed as follows
(lambda_hat <- length(hpp) / 100)

# Simulating a non-homogeneous Poisson process

# Define the intensity function (depending on `beta_0` and `beta_1`)
l <- function(x, par, log = FALSE, ...) { 
  eta <- par[1] + par[2] * x
  if (log) { res <- eta } else { res <- exp(eta) }
  res
}

pts <- seq(0, 100, 0.1)
par <- c(-1, 0.015)
plot(x = pts, y = l(pts, par), type = "l", xlab = "x", ylab = "y")

integrate(l, lower = 0, upper = 100, par = par)

# Simulate a NHPP based on the accept-reject approach
sim.NHPP <- function (par, FUN, max = 1, min = 0, ...) {
  max_lb <- optimize(FUN, c(min, max), par = par, maximum = TRUE)$objective
  lambda <- max_lb * (max - min)
  # As in `sim.IPP()`
  N <- rpois(1, lambda = lambda)
  pts <- sort(runif(N, min = min, max = max))
  # Accept or reject
  p.accept <- FUN(pts, par) / max_lb
  pp <- pts[runif(length(pts)) <= p.accept]
  attributes(pp) <- list(simulated = length(pts), accepted = length(pp), rate = length(pp) / length(pts))
  pp
}

nhpp <- sim.NHPP(par = par, FUN = l, max = 100, min = 0)
print(unlist(attributes(nhpp)))

plot1D.pp(nhpp)

# Parameter estimation based on the likelihood function
lik.NHPP <- function (par, FUN, pp, max = 1, min = 0, ...) {
  int.l <- integrate(FUN, low = min, upp = max, par = par)$value
  sum.t <- sum(FUN(x = pp, par = par, log = T))
  -(sum.t - int.l)
}

initial_values <- c(0, 0)
(theta_hat <- optim(par = initial_values, fn = lik.NHPP, FUN = l, pp = nhpp, min = 0, max = 100)[1:2])

# Hypothesis testing
L  <- psp(x0 = 0, y0 = 0, x1 = 100, y1 = 0, owin(c(0, 100), c(-5, 5)))
pp <- as.ppp(cbind(c(nhpp), 0), W = L)

q <- quadratcount(X = pp, nx = 8, ny = 1)
plot(q, main = "")

ts <- quadrat.test(X = pp, nx = 8, ny = 1)
ts

##############################
# Cox process
##############################

# Example 1

print(chorley)
plot(chorley, cols = c("red", rgb(0, 1, 0, 0.5)), pch = c(19, 4), cex = 0.75, main = "Cancer cases")

# Filter `lung cancer` cases
lung <- chorley[chorley$marks == "lung"]
lung <- ppp(x = lung$x, y = lung$y, window = lung$window)

plot(lung, cols = "green", pch = 4, cex = 0.75, main = "Lung-cancer cases")

# Map for the analyzed area and create a `raster` object
resolution <- 0.5
map <- as(st_as_sf(lung$window), "Spatial") # Convert it to a `SpatialPolygonsDataFrame` object
map$cancer <- "lung"
plot(map)

r <- raster(map, resolution = resolution) # Create a `raster` object based on the map and resolution
(n_row <- nrow(r))
(n_col <- ncol(r))

r[] <- 0 # Set all `NA` to `0`
dpts <- SpatialPoints(cbind(rev(lung$x), rev(lung$y))) # Convert the locations to a `SpatialPoints` object

# Get the cell number from a raster where `dpts` are located
(tab <- table(cellFromXY(r, dpts))) 

r[as.numeric(names(tab))] <- tab # Assign the number of observed events to the `raster` object
plot(r)
plot(map, add = T)

# Create a `grid` variable based on a `raster` object
grid <- rasterToPolygons(r) # Convert it to a `SpatialPolygonsDataFrame` object

grid <- grid[as.vector(matrix(1:nrow(grid), nrow = n_row, ncol = n_col, byrow = T)), ] # Rearrange the indices numbering

# Create and rename variables
grid$id <- 1:nrow(grid)
grid$Y <- grid$layer
grid$cellarea <- resolution * resolution
plot(grid)

# Obtain intersection only
gridmap <- raster::intersect(x = grid, y = map) # Compute the intersection between `grid` and `map`
grid <- grid[grid$id %in% gridmap$id, ]

plot(grid)
plot(map, border = "red", lwd = 1, add = T)

summary(grid)

################
# Fit model
################

# Change prior of the precision parameter from (1 / sigma^2) ~ Gamma(1, 0.0005) for PC prior for sigma
# Prob(sigma > u) = alpha
prior.list <- list(prec = list(prior = "pc.prec", param = c(0.25, 0.01))) # c(u, alpha)

formula <- Y ~ 1 + f(id, model = "matern2d", nrow = n_row, ncol = n_col, nu = 1, hyper = prior.list) # Intercept + Matérn spatial random effects

res <- inla(formula,
            family = "poisson",
            data = grid@data,
            E = cellarea) # Acts like an offset

summary(res)

# Plotting random effects
grid$R.E. <- res$summary.random$id[grid$id, "mean"]

gridborder <- gUnaryUnion(grid) # Plot the random effects using `tmap` package
tm_shape(grid) + 
  tm_polygons(col = c("R.E."), style = "cont", border.col = "transparent", midpoint = NA) +
  tm_shape(gridborder) + 
  tm_borders() +
  tm_facets(ncol = 1) + 
  tm_legend(legend.position = c("left", "bottom")) + tm_layout(fontfamily = "LM Roman 10", frame = FALSE)

# Plotting intensity
cellarea   <- resolution * resolution
grid$Mean  <- res$summary.fitted.values[, "mean"] * cellarea
grid$Lower <- res$summary.fitted.values[, "0.025quant"] * cellarea
grid$Upper <- res$summary.fitted.values[, "0.975quant"] * cellarea

# Changing the margin size to accommodate the plot caption
bbox_new <- st_bbox(grid)
xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values

bbox_new[1] <- bbox_new[1] - (0.25 * xrange)
bbox_new    <- bbox_new %>% st_as_sfc()

# Main plot for the estimated intensity (along with a 95% equal-tail credible interval)
tm_shape(grid, bbox = bbox_new) +
  tm_polygons(col = c("Lower", "Mean", "Upper"),
              style = 'fixed', border.col = "transparent",
              breaks = c(0, 0.25, 0.5, 0.75, 1.0, 2, 5, 10, ceiling(max(grid$Upper)) + 1)) +
  tm_shape(gridborder) + 
  tm_borders() +
  tm_facets(nrow = 3) + 
  tm_legend(legend.position = c("left", "bottom")) + 
  tm_layout(fontfamily = "LM Roman 10", frame = FALSE, inner.margins  = c(0, 0, 0, 0))


# Example 2 (Spatio-temporal)
terror_country <- readRDS(file = "terror_country.rds")
table(terror_country$country)

# Select Iraq (as it a place with many observations)
country_code <- "IRQ"

terror_country <- terror_country[terror_country$country == country_code, ]
terror_country <- terror_country[(terror_country$iyear >= 2010) & (terror_country$iyear <= 2015), ] 
coordinates(terror_country) <- c("longitude", "latitude")
proj4string(terror_country) <- "+proj=longlat"

area_country <- readRDS(file = "area_country.rds")
area_country <- area_country[area_country$sov_a3 == country_code, ]
area_country <- spTransform(x = area_country, CRSobj = CRS("+proj=longlat"))

plot(area_country, main = "")
plot(terror_country, add = T, col = "green")

# Do as before
resolution <- 0.5
r <- raster(area_country, resolution = resolution) 
(n_row <- nrow(r))
(n_col <- ncol(r))

# Create raster object with counting for each year (from 2010 to 2015)
terror_country$year <- terror_country$iyear - min(terror_country$iyear) + 1 
n_years <- length(unique(terror_country$year))

tab <- list()
ras <- list()
grids <- list()
grids_map <- list()

for (y in 1:n_years) {
  tab[[y]] <- table(cellFromXY(r, terror_country[terror_country$year == y, ]))
  ras[[y]] <- r
  ras[[y]][as.numeric(names(tab[[y]]))] <- tab[[y]]
  values(ras[[y]])[is.na(values(ras[[y]]))] <- 0
  grids[[y]] <- rasterToPolygons(ras[[y]]) 
  grids[[y]] <- grids[[y]][as.vector(matrix(1:nrow(grids[[y]]), nrow = n_row, ncol = n_col, byrow = T)), ] 
  grids[[y]]$id <- 1:nrow(grids[[y]])
  grids[[y]]$Y <- grids[[y]]$layer
  grids[[y]]$cellarea <- resolution * resolution
  grids_map[[y]] <- raster::intersect(x = grids[[y]], y = area_country) 
  grids[[y]] <- grids[[y]][grids[[y]]$id %in% grids_map[[y]]$id, ]
  
  plot(ras[[y]], main = y)
  plot(area_country, add = T)
}

# Create `data_inla` object
for (y in 1:n_years) {
  if (y == 1) {
    data_inla <- grids[[y]]@data
  } else {
    data_inla <- rbind(data_inla, grids[[y]]@data)
  }
}
data_inla <- cbind(data_inla, id_time = rep(x = 1:n_years, each = nrow(grids[[1]])))
data_inla[c(1:3, ((nrow(data_inla) - 2):nrow(data_inla))), ]

################
# Fit model
################

prior.list <- list(prec = list(prior = "pc.prec", param = c(0.25, 0.01))) # Priors

# Separable space-temporal model: eta_ij = beta_0 + u_i + v_j
formula <- Y ~ 1 + f(id, 
                     model = "matern2d", 
                     nrow = n_row, 
                     ncol = n_col, 
                     nu = 1, 
                     hyper = prior.list) + f(id_time, model = "ar1")

# Model with interaction: eta_ij = beta_0 + u_ij
# formula <- Y ~ 1 + f(id, 
#                      model = "matern2d", 
#                      nrow = n_row, 
#                      ncol = n_col, 
#                      nu = 1, 
#                      group = id_time,
#                      control.group = list(model = "ar1"),
#                      hyper = prior.list)

res <- inla(formula,
            family = "poisson",
            data = data_inla,
            E = resolution)

summary(res)

# Plot results
grid       <- grids[[1]]
cells_grid <- nrow(grids[[1]])
cellarea   <- resolution * resolution

grid$M1 <- res$summary.fitted.values[, "mean"][1:cells_grid] * cellarea
grid$M2 <- res$summary.fitted.values[, "mean"][(cells_grid + 1):(cells_grid * 2)] * cellarea
grid$M3 <- res$summary.fitted.values[, "mean"][(cells_grid * 2 + 1):(cells_grid * 3)] * cellarea
grid$M4 <- res$summary.fitted.values[, "mean"][(cells_grid * 3 + 1):(cells_grid * 4)] * cellarea
grid$M5 <- res$summary.fitted.values[, "mean"][(cells_grid * 4 + 1):(cells_grid * 5)] * cellarea
grid$M6 <- res$summary.fitted.values[, "mean"][(cells_grid * 5 + 1):(cells_grid * 6)] * cellarea

max_int <- ceiling(max(c(grid$M1, grid$M2, grid$M3, grid$M4, grid$M5, grid$M6)))
max_int <- ceiling(max_int / 10) * 10

bbox_new <- st_bbox(grid)
xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values

bbox_new[1] <- bbox_new[1] - (0.25 * xrange)
bbox_new <- bbox_new %>% st_as_sfc()

gridborder <- gUnaryUnion(grid) 
tm_shape(grid, bbox = bbox_new) +
  tm_polygons(col = c("M1", "M2", "M3", "M4", "M5", "M6"),
              style = 'fixed', border.col = "transparent",
              breaks = c(0, 1, 3, 5, 10, 20, 30, 50, max_int)) + 
  tm_shape(gridborder) + tm_borders() +
  tm_facets(ncol = 3) + tm_legend(legend.position = c("left", "bottom")) +
  tm_layout(fontfamily = "LM Roman 10", frame = FALSE)
