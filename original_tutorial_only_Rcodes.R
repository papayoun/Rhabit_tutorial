## ----setup, include = FALSE----------------------------------------------
rm(list = ls())
knitr::opts_chunk$set(echo = TRUE, comment = NA, cache = TRUE)


## ----install_rhabit, eval = FALSE----------------------------------------
## devtools::install_github("papayoun/Rhabit") # Installing (run only once)


## ----load_packages, message = FALSE, cache = FALSE-----------------------
library(Rhabit) # Package dedicated to the Langevin movement model
library(tidyverse) # Set of useful packages, including ggplot2
library(raster)


## ----trajectory_data, echo = -2------------------------------------------
trajectory_data <- read.table(file = "tutorial_simulated_data.txt",
                              sep = "\t", header = TRUE,
                              colClasses = c(ID = "factor")) # ID column is a factor
knitr::kable(head(trajectory_data))


## ----summary_trajectory_data, comment = NA-------------------------------
str(trajectory_data)


## ----plotting_simulated_data---------------------------------------------
ggplot(trajectory_data) + # Data set which will be plotted
  aes(x = x, y = y) + # Axis
  geom_path(aes(color = ID)) + # Plotting trajectories, One color per trajectory
  # We add starting points
  geom_point(data = dplyr::filter(trajectory_data, t == 0), col = "blue") + 
  theme(legend.position = "none") + # No need for a legend
  labs(x = "X-axis", y = "Y-axis", title = "Trajectory data set")


## ----environment_covariates----------------------------------------------
environment_covariates <- get(load("environment_covariates.RData"))


## ----properties_env_cov--------------------------------------------------
is(environment_covariates) # A list
(J <- length(environment_covariates)) # J is the number of covariates


## ----str_environment_covariates------------------------------------------
str(environment_covariates)


## ----plot_covariates-----------------------------------------------------
plotCovariates(environment_covariates)


## ----plot_covariates_plus_traj, echo = FALSE-----------------------------
# Same plot, adding trajectory
plotCovariates(environment_covariates, trajectory_data = trajectory_data)


## ----get_covariates, eval = FALSE----------------------------------------
## # Gathering all covariates in data.frame
## covariates_df <- map_dfr(environment_covariates, # To each element of the list, apply
##                          Rhabit::rasterToDataFrame, # This function
##                          .id = "Covariate")


## ----gradient_array------------------------------------------------------
locations <- trajectory_data %>% # From trajectory_data, extract 
  dplyr::select(x, y) %>% # the  sampled locations  (the x and y columns)
  as.matrix() # and convert it to a matrix 
gradient_array <- Rhabit::bilinearGradArray(locs = locations,
                                            cov_list = environment_covariates)


## ----langevin_fitting----------------------------------------------------
langevin_fit <- Rhabit::langevinUD(locs = locations, # Matrix of locations 
                                   times = trajectory_data$t, # Vector of times
                                   ID = as.character(trajectory_data$ID), # Vector of IDs
                                   grad_array = gradient_array) # Gradient at locations


## ----langevin_fit_elements-----------------------------------------------
str(langevin_fit)


## ----betaHat_first_model-------------------------------------------------
langevin_fit$betaHat


## ----CI_first_model------------------------------------------------------
langevin_fit$CI


## ----AIC_first_model-----------------------------------------------------
langevin_fit$AIC


## ----residuals-----------------------------------------------------------
langevin_fit$residuals


## ----predicted-----------------------------------------------------------
# Note that the first value of each trajectory cannot be predicted
langevin_fit$predicted


## ----plot_predicted------------------------------------------------------
langevin_fit$predicted %>% 
  rename(x_pred = x, y_pred = y) %>%
  select(x_pred, y_pred) %>% 
  bind_cols(trajectory_data) %>% 
  filter(ID == "1") %>% 
  slice(1:30) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_point() +
  geom_point(data = slice(trajectory_data, 1), # Mark the first point
             col = "purple", size = 5)+ 
  geom_path(aes(group = ID), color = "blue", linetype = 2) + # Plot the path
  geom_segment(aes(xend = lead(x_pred), # Plot the expected displacement
                   yend = lead(y_pred)), # With arrows
               color = "red", 
               arrow = arrow(length = unit(0.01, "npc")))


## ----plot_ud-------------------------------------------------------------
Rhabit::plotUD(environment_covariates, langevin_fit$betaHat) +
  # geom_path(data = trajectory_data, aes(group = ID)) + # If you want to add trajectories
  labs(x = "X-axis", y = "Y-axis",
       title = "Estimated UD")


## ----get_all_covariates--------------------------------------------------
distance_covariate <- get(load("distance_covariate.RData"))
all_covariates <- c(environment_covariates,
                    distance_covariate)


## ----seal_tracks---------------------------------------------------------
seal_tracks <- read.table(file = "seal_data.csv", header = T, sep = ",",
                          colClasses = c(ID = "factor", # Specifying column classes
                                         time = "POSIXct")) %>% 
  mutate(x = x / 1000,
         y = y / 1000)


## ----plotting_ssl_data---------------------------------------------------
ggplot(seal_tracks) +
  aes(x = x, y = y, col = ID) +
  geom_path() +
  labs(x = "Easting (meters)", 
       y = "Northing (meters)")


## ----head_seal_tracks----------------------------------------------------
seal_tracks_Rhabit <- seal_tracks %>% # 
  group_by(ID) %>% # Applying same treatment to each group
  # Changing the time column, in an increasing numeric vector (here in hours)
  mutate(time = c(0, cumsum(diff(as.numeric(time))))) %>% 
  ungroup()


## ----loading_ssl_covariates----------------------------------------------
habitat_covariates <- raster::brick("habitat_seal_data/aleut_habitat.grd", 
                                    values = TRUE) %>%
  # Raster brick, in some projection
  raster::as.list() %>% # Turns it into a list
  purrr::map(Rhabit::rasterToRhabit) %>%  # To each list element, we apply the rasterToRhabit function
  purrr::map(function(covariable){
    covariable$x <- covariable$x / 1000
    covariable$y <- covariable$y / 1000
    return(covariable)
  })


## ----plot_ssl_covariates-------------------------------------------------
Rhabit::plotCovariates(habitat_covariates)


## ----fitting_langevin_sll------------------------------------------------
# First, get locations as a matrix
locations_ssl <- seal_tracks_Rhabit %>% # From trajectory_data, extract 
  dplyr::select(x, y) %>% # the  sampled locations  (the x and y columns)
  as.matrix() # and convert it to a matrix 
# Then, get gradients at each gradient array
gradient_array_ssl <- Rhabit::bilinearGradArray(locs = locations_ssl,
                                            cov_list = habitat_covariates)
# Then, fit it!
langevin_fit_ssl <- langevinUD(locs = locations_ssl,
                               times = seal_tracks_Rhabit$time,
                               ID = seal_tracks_Rhabit$ID, 
                               grad_array = gradient_array_ssl)


## ----estimates_ssl-------------------------------------------------------
langevin_fit_ssl$betaHat


## ----CI_ssl--------------------------------------------------------------
langevin_fit_ssl$CI


## ----cor_estimates-------------------------------------------------------
# The variance-covariance matrix is transformed to a correlation matrix
cov2cor(langevin_fit_ssl$betaHatVariance)


## ----plot_UD-------------------------------------------------------------
plotUD(covariates = habitat_covariates, beta = langevin_fit_ssl$betaHat)


## ----checking_predicted_ssl----------------------------------------------
langevin_fit_ssl$predicted %>% 
  rename(x_pred = x, y_pred = y) %>%
  select(x_pred, y_pred) %>% 
  bind_cols(seal_tracks) %>% 
  filter(ID == "14809") %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_point() +
  geom_path(aes(group = ID), color = "blue", linetype = 2) + # Plot the path
  geom_segment(aes(xend = lead(x_pred), # Plot the expected displacement
                   yend = lead(y_pred)), # With arrows
               color = "red", 
               arrow = arrow(length = unit(0.01, "npc"))) +
  coord_cartesian(xlim = c(-1780, -1760), ylim = c(560, 590))

