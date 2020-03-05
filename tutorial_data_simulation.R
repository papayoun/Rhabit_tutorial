rm(list = ls())
get_distance_square <- function(lim, resol, center = c(0, 0),
                                raster_like = T){
  xgrid <- seq(lim[1], lim[2], by = resol)
  ygrid <- seq(lim[3], lim[4], by = resol)
  coords <- as.matrix(expand.grid(xgrid, ygrid))
  distance_to_center <- apply(coords, 1, function(z) sum((z - center)^2))
  if (raster_like) {
    return(list(x = xgrid, y = ygrid, 
                z = matrix(distance_to_center, nrow = length(xgrid))))
  }
  return(data.frame(x = coords[, 1], y = coords[, 2], val = distance_to_center))
}

set.seed(1)# repeatability
lim <- c(-10, 10, -10, 10) # limits of map
resol <- 0.1 # grid resolution
rho <- 4; nu <- 1; sigma2 <- 5# Matern covariance parameters
mean_function <- function(z){# mean function
  -log(3 + sum(z^2))
}
J <- 3 # number of covariates

 
my_covariates <- c(replicate(J, # Number of environmental covariates 
                          Rhabit::simSpatialCov(lim, nu, rho, sigma2, resol = resol,
                                                mean_function = mean_function,
                                                raster_like = T),
                          simplify = F),
                list(get_distance_square(lim, resol)))
environment_covariates <- my_covariates[1:J]
distance_covariate <- my_covariates[J + 1]
save(environment_covariates, 
     file = "environment_covariates.RData")
save(distance_covariate, 
     file = "distance_covariate.RData")

beta_true <- c(5, 10, 0, -1.5)
plot_ud <- Rhabit::getUD(my_covariates, beta_true, log = F) %>% 
  Rhabit::rasterToGGplot() %>% 
  ggplot(aes(x = x, y = y)) +
  geom_raster(aes(fill = val)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c() +
  labs(fill = expression(pi(z)), x = expression(z[1]), y = expression(z[2]))


t_final <- 1
time_step <- 1 / 1000
times <- seq(0, t_final, by = time_step)
set.seed(123) #repeatability
tracks <- purrr::rerun(20,
                       Rhabit::simLangevinMM(beta = beta_true, gamma2 = 5,
                        times = times, loc0 = runif(2, -9, 9),
                        cov_list = my_covariates, keep_grad = F)) %>% 
  bind_rows(.id = "ID")
write.table(tracks, file = "tutorial_simulated_data.txt",
            row.names = FALSE, sep = "\t")

plot_ud + 
  geom_path(data = tracks, aes(color = id)) +
  geom_point(data = tracks %>% group_by(id) %>% 
               summarise(x = first(x), y = first(y)),
             col = "blue")
locs <- tracks %>% 
  select(x, y) %>% 
  as.matrix()
times <- pull(tracks, t)
ID <- pull(tracks, id)
grad_array <- Rhabit::bilinearGradArray(locs, my_covariates)
langevin_result <- Rhabit::langevinUD(locs = locs, times = times, ID = ID, 
                                      leverage = F,
                                      grad_array = grad_array[,,1:2, drop = F])

norm_res <- apply(langevin_result$residuals,1,function(x) sum(x^2))

residuals <- langevin_result$residuals
head(residuals)
angle <- complex(real = residuals[, 1], 
                 imaginary = residuals[, 2]) %>% 
  Arg() %>% 
  {(. *180 / pi)}

tracks %>% rowid_to_column(var = "num") %>% 
  filter((num %% 100) !=0) %>% 
  mutate(angle = angle) %>% 
  ggplot() + 
  geom_point(aes(x = x, y = y, color = angle)) +
  scale_color_viridis_c()
