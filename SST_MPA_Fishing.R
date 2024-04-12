##Thermal effects x Mgt strategies

library(tidyverse)
library(ggpubr)
library(RColorBrewer)

source('dispersal_fn.R')

projections <- read.csv("anomaly_df.csv")
SST_dev <- projections[['anomaly']]

timesteps <- length(SST_dev)
r = 0.3
K = c(101.3, 101.3)
S = 0.5


patch_area_sequences <- list(seq(0, 1, by = 0.05))
patch_area_grid <- do.call(expand.grid, patch_area_sequences)
patch_area_grid$Var2 <- 1 - patch_area_grid$Var1
patch_area_list <- split(patch_area_grid, 1:nrow(patch_area_grid))
number_patches <- ncol(patch_area_grid)

fishing_effort_sequences <- list(seq(0, 0.5, by = 0.01), 0)
fishing_effort_grid <- do.call(expand.grid, fishing_effort_sequences)
fishing_effort_list <- split(fishing_effort_grid, 1:nrow(fishing_effort_grid))
catchability <- 1

parameter_grid <- expand.grid(patch_area = patch_area_list,
                              fishing_effort = fishing_effort_list)

#saving outputs
outcome_population <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_notemp <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

outcome_population_r1 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_r1 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

outcome_population_r2 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_r2 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

outcome_population_r3 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_r3 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

outcome_population_K1 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_K1 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

outcome_population_K2 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_K2 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))




population <- array(NA, dim = c(timesteps, number_patches))
harvest <- array(NA, dim = c(timesteps, number_patches))

population_r1 <- array(NA, dim = c(timesteps, number_patches))
harvest_r1 <- array(NA, dim = c(timesteps, number_patches))

population_r2 <- array(NA, dim = c(timesteps, number_patches))
harvest_r2 <- array(NA, dim = c(timesteps, number_patches))

population_r3 <- array(NA, dim = c(timesteps, number_patches))
harvest_r3 <- array(NA, dim = c(timesteps, number_patches))

population_K1 <- array(NA, dim = c(timesteps, number_patches))
harvest_K1 <- array(NA, dim = c(timesteps, number_patches))

population_K2 <- array(NA, dim = c(timesteps, number_patches))
harvest_K2 <- array(NA, dim = c(timesteps, number_patches))


r_temp_1 <- array(NA, dim = c(timesteps, 1))
r_temp_2 <- array(NA, dim = c(timesteps, 1))
r_temp_3 <- array(NA, dim = c(timesteps, 1))

K_temp_1 <- array(NA, dim = c(timesteps, 1))
K_temp_2 <- array(NA, dim = c(timesteps, 1))



for (iter in 1:nrow(parameter_grid)){
  
  #dispersal matrix
  tmp <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S)
  
  #starting numbers
  population[1,] <- c(10,10)
  population_r1[1,] <- c(10,10)
  population_r2[1,] <- c(10,10)
  population_r3[1,] <- c(10,10)
  population_K1[1,] <- c(10,10)
  population_K2[1,] <- c(10,10)
  
  for (t in 2:timesteps) {
    
    #no temp
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) # find which patch is above K
      population[t-1, patch_zero]     <- 0 # force patch above K to equal K
    }
    #harvest and population growth
    harvest[t,] <- population[t-1,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population[t,] <- (population[t-1,] + r * population[t-1,] * (1 - population[t-1,] / K) - harvest[t,]) %*% tmp
    
    #if its above K after dispersal, set it back down to K
    if(population[t, 1] > K[1] & population[t, 2] > K[2]){
      population[t, ] <- K
    }
    else if(population[t, 1] > K[1] | population[t, 2] > K[2]){
      patch_above     <- which(population[t, ] > K) # find which patch is above K
      patch_not_above <- which(!(population[t, ] > K)) # find patch not above K
      spillover       <- population[t, patch_above] -  K[patch_above] # set spillover to the difference between population and K
      
      population[t, patch_above]     <- K[patch_above] # force patch above K to equal K
      population[t, patch_not_above] <- population[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population[t, patch_not_above] > K[patch_not_above]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population[t, patch_not_above] <- K[patch_not_above] # if it does then set that patch to carrying capacity after spillover
      }
    } 
    
    
    #temp dependent r - "current" temp as optimal
    r_temp_1[t] <- 0.3 + 0*SST_dev[t] + -0.0037*SST_dev[t]^2
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) # find which patch is above K
      population_r1[t-1, patch_zero]     <- 0 # force patch above K to equal K
    }
    
    #harvest and population growth
    harvest_r1[t,] <- population_r1[t-1,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_r1[t,] <- (population_r1[t-1,] + r_temp_1[t] * population_r1[t-1,] * 
                            (1 - population_r1[t-1,] / K) - harvest_r1[t,]) %*% tmp
    
    if(population_r1[t, 1] > K[1] & population_r1[t, 2] > K[2]){
      population_r1[t, ] <- K
    }
    else if(population_r1[t, 1] > K[1] | population_r1[t, 2] > K[2]){
      patch_above     <- which(population_r1[t, ] > K) # find which patch is above K
      patch_not_above <- which(!(population_r1[t, ] > K)) # find patch not above K
      spillover       <- population_r1[t, patch_above] -  K[patch_above] # set spillover to the difference between population and K
      
      population_r1[t, patch_above]     <- K[patch_above] # force patch above K to equal K
      population_r1[t, patch_not_above] <- population_r1[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population_r1[t, patch_not_above] > K[patch_not_above]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population_r1[t, patch_not_above] <- K[patch_not_above] # if it does then set that patch to carrying capacity after spillover
      }
    } 
    
    
    #temp dependent r - higher optimal temp
    r_temp_2[t] <- 0.3 + 0*(SST_dev[t] - 1) + -0.0037*(SST_dev[t] - 1 )^2
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) # find which patch is above K
      population_r2[t-1, patch_zero]     <- 0 # force patch above K to equal K
    }
    
    harvest_r2[t,] <- population_r2[t-1,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_r2[t,] <- (population_r2[t-1,] + r_temp_2[t] * population_r2[t-1,] * 
                            (1 - population_r2[t-1,] / K) - harvest_r2[t,]) %*% tmp
    
    if(population_r2[t, 1] > K[1] & population_r2[t, 2] > K[2]){
      population_r2[t, ] <- K
    }
    else if(population_r2[t, 1] > K[1] | population_r2[t, 2] > K[2]){
      patch_above     <- which(population_r2[t, ] > K) # find which patch is above K
      patch_not_above <- which(!(population_r2[t, ] > K)) # find patch not above K
      spillover       <- population_r2[t, patch_above] -  K[patch_above] # set spillover to the difference between population and K
      
      population_r2[t, patch_above]     <- K[patch_above] # force patch above K to equal K
      population_r2[t, patch_not_above] <- population_r2[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population_r2[t, patch_not_above] > K[patch_not_above]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population_r2[t, patch_not_above] <- K[patch_not_above] # if it does then set that patch to carrying capacity after spillover
      }
    } 
    
    
    #temp dependent r - lower optimal temp
    r_temp_3[t] <- 0.3 + 0*(SST_dev[t] + 1) + -0.0037*(SST_dev[t] + 1 )^2
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) 
      population_r3[t-1, patch_zero]     <- 0
    }
    
    harvest_r3[t,] <- population_r3[t-1,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_r3[t,] <- (population_r3[t-1,] + r_temp_3[t] * population_r3[t-1,] * 
                            (1 - population_r3[t-1,] / K) - harvest_r3[t,]) %*% tmp
    
    if(population_r3[t, 1] > K[1] & population_r3[t, 2] > K[2]){
      population_r3[t, ] <- K
    }
    else if(population_r3[t, 1] > K[1] | population_r3[t, 2] > K[2]){
      patch_above     <- which(population_r3[t, ] > K) # find which patch is above K
      patch_not_above <- which(!(population_r3[t, ] > K)) # find patch not above K
      spillover       <- population_r3[t, patch_above] -  K[patch_above] # set spillover to the difference between population and K
      
      population_r3[t, patch_above]     <- K[patch_above] # force patch above K to equal K
      population_r3[t, patch_not_above] <- population_r3[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population_r3[t, patch_not_above] > K[patch_not_above]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population_r3[t, patch_not_above] <- K[patch_not_above] # if it does then set that patch to carrying capacity after spillover
      }
    } 
    
    
    #temp dependent K - piece wise linear function
    K_temp_1[t] <- -4.95243768 * SST_dev[t] + 101.3
    if (K_temp_1[t] < 10) {
      K_temp_1[t] <- 10
    } else if (K_temp_1[t] > 101.3) {
      K_temp_1[t] <- 101.3
    }
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) 
      population_K1[t-1, patch_zero]     <- 0
    }
    
    harvest_K1[t,] <- population_K1[t-1,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_K1[t,] <- (population_K1[t-1,] + r * population_K1[t-1,] * 
                            (1 - population_K1[t-1,] / K_temp_1[t]) - harvest_K1[t,]) %*% tmp
    
    
    
    if(population_K1[t, 1] > K_temp_1[t] & population_K1[t, 2] > K_temp_1[t]){
      population_K1[t, ] <- K_temp_1[t]
    }
    else if(population_K1[t, 1] > K_temp_1[t] | population_K1[t, 2] > K_temp_1[t]){
      patch_above     <- which(population_K1[t, ] > K_temp_1[t]) # find which patch is above K
      patch_not_above <- which(!(population_K1[t, ] > K_temp_1[t])) # find patch not above K
      spillover       <- population_K1[t, patch_above] - K_temp_1[t] # set spillover to the difference between population and K
      
      population_K1[t, patch_above]     <- K_temp_1[t] # force patch above K to equal K
      population_K1[t, patch_not_above] <- population_K1[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population_K1[t, patch_not_above] > K_temp_1[t]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population_K1[t, patch_not_above] <- K_temp_1[t] # if it does then set that patch to carrying capacity after spillover
      }
    }
    if(as.numeric(parameter_grid[['patch_area']][[iter]][,1]) == 0){
      population_K1[t,1] <- 0
    } else if (as.numeric(parameter_grid[['patch_area']][[iter]][,2]) == 0) {
      population_K1[t,2] <- 0
    }
    
    
    #temp dependent K - quadratic function
    K_temp_2[t] <- 101.3 + 0*SST_dev[t] + -0.7*SST_dev[t]^2
    if (K_temp_2[t] < 10) {
      K_temp_2[t] <- 10 }
    else if (K_temp_2[t] > 101.3) {
      K_temp_2[t] <- 101.3
    }
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0)
      population_K2[t-1, patch_zero]     <- 0 
    }
    
    harvest_K2[t,] <- population_K2[t-1,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_K2[t,] <- (population_K2[t-1,] + r * population_K2[t-1,] * 
                            (1 - population_K2[t-1,] / K_temp_2[t]) - harvest_K2[t,]) %*% tmp
    
    if(population_K2[t, 1] > K_temp_2[t] & population_K2[t, 2] > K_temp_2[t]){
      population_K2[t, ] <- K_temp_2[t]
    }
    else if(population_K2[t, 1] > K_temp_2[t] | population_K2[t, 2] > K_temp_2[t]){
      patch_above     <- which(population_K2[t, ] > K_temp_2[t]) # find which patch is above K
      patch_not_above <- which(!(population_K2[t, ] > K_temp_2[t])) # find patch not above K
      spillover       <- population_K2[t, patch_above] - K_temp_2[t] # set spillover to the difference between population and K
      
      population_K2[t, patch_above]     <- K_temp_2[t] # force patch above K to equal K
      population_K2[t, patch_not_above] <- population_K2[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population_K2[t, patch_not_above] > K_temp_2[t]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population_K2[t, patch_not_above] <- K_temp_2[t] # if it does then set that patch to carrying capacity after spillover
      }
    }
    if(as.numeric(parameter_grid[['patch_area']][[iter]][,1]) == 0){
      population_K2[t,1] <- 0
    } else if (as.numeric(parameter_grid[['patch_area']][[iter]][,2]) == 0) {
      population_K2[t,2] <- 0
    }
  }
  
  outcome_population[iter, ] <- colMeans(population[(t-19):t,])
  outcome_harvest_notemp[iter, ] <- colMeans(harvest[(t-19):t,])
  
  outcome_population_r1[iter, ] <- colMeans(population_r1[(t-19):t,])
  outcome_harvest_r1[iter, ] <- colMeans(harvest_r1[(t-19):t,])
  
  outcome_population_r2[iter, ] <- colMeans(population_r2[(t-19):t,])
  outcome_harvest_r2[iter, ] <- colMeans(harvest_r2[(t-19):t,])
  
  outcome_population_r3[iter, ] <- colMeans(population_r3[(t-19):t,])
  outcome_harvest_r3[iter, ] <- colMeans(harvest_r3[(t-19):t,])
  
  outcome_population_K1[iter, ] <- colMeans(population_K1[(t-19):t,])
  outcome_harvest_K1[iter, ] <- colMeans(harvest_K1[(t-19):t,])
  
  outcome_population_K2[iter, ] <- colMeans(population_K2[(t-19):t,])
  outcome_harvest_K2[iter, ] <- colMeans(harvest_K2[(t-19):t,])
  
}


colnames(outcome_population) <- c("open_fish_notemp", "mpa_fish_notemp")
outcome_population <- as.data.frame(outcome_population)
colnames(outcome_harvest_notemp) <- c("open_harvest_notemp", "mpa_harvest_notemp")

colnames(outcome_population_r1) <- c("open_fish_r1", "mpa_fish_r1")
outcome_population_r1 <- as.data.frame(outcome_population_r1)
colnames(outcome_harvest_r1) <- c("open_harvest_r1", "mpa_harvest_r1")

colnames(outcome_population_r2) <- c("open_fish_r2", "mpa_fish_r2")
outcome_population_r2 <- as.data.frame(outcome_population_r2)
colnames(outcome_harvest_r2) <- c("open_harvest_r2", "mpa_harvest_r2")

colnames(outcome_population_r3) <- c("open_fish_r3", "mpa_fish_r3")
outcome_population_r3 <- as.data.frame(outcome_population_r3)
colnames(outcome_harvest_r3) <- c("open_harvest_r3", "mpa_harvest_r3")

colnames(outcome_population_K1) <- c("open_fish_K1", "mpa_fish_K1")
outcome_population_K1 <- as.data.frame(outcome_population_K1)
colnames(outcome_harvest_K1) <- c("open_harvest_K1", "mpa_harvest_K1")

colnames(outcome_population_K2) <- c("open_fish_K2", "mpa_fish_K2")
outcome_population_K2 <- as.data.frame(outcome_population_K2)
colnames(outcome_harvest_K2) <- c("open_harvest_K2", "mpa_harvest_K2")

#Fish population dataframe
outcome <- cbind(parameter_grid, outcome_population, outcome_population_r1, outcome_population_r2, outcome_population_r3,
                 outcome_population_K1, outcome_population_K2)

outcome$area_open <- map_dbl(outcome$patch_area, 1)
outcome$area_mpa <- map_dbl(outcome$patch_area, 2)
outcome$fishing_p1 <- map_dbl(outcome$fishing_effort, 1)
outcome$fishing_p2 <- map_dbl(outcome$fishing_effort, 2)

##Harvest
outcome_harvest <- cbind(parameter_grid, outcome_harvest_notemp, outcome_harvest_r1, outcome_harvest_r2, outcome_harvest_r3,
                         outcome_harvest_K1, outcome_harvest_K2)

outcome_harvest$area_open <- map_dbl(outcome_harvest$patch_area, 1)
outcome_harvest$area_mpa <- map_dbl(outcome_harvest$patch_area, 2)
outcome_harvest$fishing_p1 <- map_dbl(outcome_harvest$fishing_effort, 1)
outcome_harvest$fishing_p2 <- map_dbl(outcome_harvest$fishing_effort, 2)




#RELATIVE TO OPTIMAL VALUES
#optimal nonspatial would be MPA = 0 and fishing at MSY

#from testingmsy.R
optimal_fishing_effort = 0.14
optimal_biomass = 5.047030e+01
optimal_harvest = 7.584363e+00
optimal_CPUE = 7.584363e+00 / 0.14

#######################. BIOMASS. ##################################
#optimal_biomass <- outcome %>% filter(area_mpa == 0 & fishing_p1 == optimal_fishing_effort) %>%
# summarize(open_fish_notemp, open_fish_r1, open_fish_r2, open_fish_r3, open_fish_K1, open_fish_K2)

biomass_mapping <- data.frame(
  fishing_p1 = unique(outcome$fishing_p1)
)
biomass_mapping$optimal_biomass <- optimal_biomass

outcome_index <- merge(outcome, biomass_mapping, by = "fishing_p1")

outcome_index$combined_fish_notemp <- (outcome_index$open_fish_notemp * outcome_index$area_open) + 
  (outcome_index$mpa_fish_notemp * outcome_index$area_mpa)
outcome_index$combined_fish_r1 <- (outcome_index$open_fish_r1 * outcome_index$area_open) + 
  (outcome_index$mpa_fish_r1 * outcome_index$area_mpa)
outcome_index$combined_fish_r2 <- (outcome_index$open_fish_r2 * outcome_index$area_open) + 
  (outcome_index$mpa_fish_r2 * outcome_index$area_mpa)
outcome_index$combined_fish_r3 <- (outcome_index$open_fish_r3 * outcome_index$area_open) + 
  (outcome_index$mpa_fish_r3 * outcome_index$area_mpa)
outcome_index$combined_fish_K1 <- (outcome_index$open_fish_K1 * outcome_index$area_open) + 
  (outcome_index$mpa_fish_K1 * outcome_index$area_mpa)
outcome_index$combined_fish_K2 <- (outcome_index$open_fish_K2 * outcome_index$area_open) + 
  (outcome_index$mpa_fish_K2 * outcome_index$area_mpa)


outcome_combined <- outcome_index %>%
  select(combined_fish_notemp, combined_fish_r1, combined_fish_r2, combined_fish_r3, combined_fish_K1, combined_fish_K2,
         optimal_biomass, fishing_p1, area_mpa)

outcome_combined$rel_biomass_notemp <- outcome_combined$combined_fish_notemp / outcome_combined$optimal_biomass
outcome_combined$rel_biomass_r1 <- outcome_combined$combined_fish_r1 / outcome_combined$optimal_biomass
outcome_combined$rel_biomass_r2 <- outcome_combined$combined_fish_r2 / outcome_combined$optimal_biomass
outcome_combined$rel_biomass_r3 <- outcome_combined$combined_fish_r3 / outcome_combined$optimal_biomass
outcome_combined$rel_biomass_K1 <- outcome_combined$combined_fish_K1 / outcome_combined$optimal_biomass
outcome_combined$rel_biomass_K2 <- outcome_combined$combined_fish_K2 / outcome_combined$optimal_biomass

outcome_combined <- outcome_combined %>%
  select(rel_biomass_notemp, rel_biomass_r1, rel_biomass_r2, rel_biomass_r3, rel_biomass_K1,
         rel_biomass_K2, fishing_p1, area_mpa)


outcome_combined_long <- pivot_longer(outcome_combined, rel_biomass_notemp:rel_biomass_K2, names_to = "model_version",
                                      values_to = "rel_biomass")

outcome_combined_long$model_version <- factor(outcome_combined_long$model_version, 
                                              levels = c("rel_biomass_notemp", "rel_biomass_r1", "rel_biomass_r2", 
                                                         "rel_biomass_r3", "rel_biomass_K1", "rel_biomass_K2"),
                                              labels = c("Baseline", "r1", "r2", 
                                                         "r3", "K1", "K2"))

##RELATIVE HARVEST

harvest_mapping <- data.frame(
  fishing_p1 = unique(outcome$fishing_p1)
)

harvest_mapping$optimal_harvest <- optimal_harvest

outcome_harvest_index <- merge(outcome_harvest, harvest_mapping, by = "fishing_p1")

outcome_harvest_2 <- outcome_harvest_index %>%
  select(open_harvest_notemp, open_harvest_r1, open_harvest_r2, open_harvest_r3, 
         open_harvest_K1, open_harvest_K2,
         optimal_harvest,
         fishing_p1, area_mpa)

outcome_harvest_2$rel_open_harvest_notemp <- outcome_harvest_2$open_harvest_notemp / outcome_harvest_2$optimal_harvest
outcome_harvest_2$rel_open_harvest_r1 <- outcome_harvest_2$open_harvest_r1 / outcome_harvest_2$optimal_harvest
outcome_harvest_2$rel_open_harvest_r2 <- outcome_harvest_2$open_harvest_r2 / outcome_harvest_2$optimal_harvest
outcome_harvest_2$rel_open_harvest_r3 <- outcome_harvest_2$open_harvest_r3 / outcome_harvest_2$optimal_harvest
outcome_harvest_2$rel_open_harvest_K1 <- outcome_harvest_2$open_harvest_K1 / outcome_harvest_2$optimal_harvest
outcome_harvest_2$rel_open_harvest_K2 <- outcome_harvest_2$open_harvest_K2 / outcome_harvest_2$optimal_harvest

outcome_harvest_2 <- outcome_harvest_2 %>%
  select(rel_open_harvest_notemp, rel_open_harvest_r1, rel_open_harvest_r2, rel_open_harvest_r3, rel_open_harvest_K1,
         rel_open_harvest_K2, fishing_p1, area_mpa)


outcome_harvest_long <- pivot_longer(outcome_harvest_2, rel_open_harvest_notemp:rel_open_harvest_K2, names_to = "model_version",
                                     values_to = "rel_harvest")

outcome_harvest_long$model_version <- factor(outcome_harvest_long$model_version, 
                                             levels = c("rel_open_harvest_notemp", "rel_open_harvest_r1", "rel_open_harvest_r2", 
                                                        "rel_open_harvest_r3", "rel_open_harvest_K1", "rel_open_harvest_K2"),
                                             labels = c("Baseline", "r1", "r2", 
                                                        "r3", "K1", "K2"))


#First exploring no MPA, just varied fishing effort
my_colors <- c("#000000", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")


#this is figure 3
p1 <- ggplot(outcome_combined_long %>% filter(area_mpa == 0), aes(x = fishing_p1, y = rel_biomass, color = model_version)) +
  geom_line(lwd = 0.75) +
  scale_color_manual(values = my_colors, name = "Model version") + 
  labs(x = "Fishing effort", y = "Relative biomass") +
  theme_minimal() +
  theme(text = element_text(size = 20), legend.position = "bottom")


p2<-ggplot(outcome_harvest_long %>% filter(area_mpa == 0), aes(x=fishing_p1, y = rel_harvest, color=model_version)) +
  geom_line(lwd=0.75) +
  scale_color_manual(values = my_colors, name = "Model version") + 
  labs(x = "Fishing effort", y = "Relative harvest") +
  theme_minimal() +
  theme(text = element_text(size=20), legend.position = "bottom")

ggarrange(p1, p2, nrow=1, ncol=2, common.legend = TRUE) # 8 x 4


#now with MPAs

##Fishing effort is X axis, facetted by MPA level
outcome_combined_long <- outcome_combined_long %>% 
  mutate(area_mpa = round(area_mpa, 2))

outcome_harvest_long <- outcome_harvest_long %>% 
  mutate(area_mpa = round(area_mpa, 2))


# Define custom labels for the facets
custom_labels <- c(`0.05` = "5% spatial closure",
                   `0.1` = "10% spatial closure",
                   `0.3` = "30% spatial closure",
                   `0.5` = "50% spatial closure")


#this is figure 4
ggplot(outcome_combined_long %>%
         filter(area_mpa %in% c(0.05, 0.1, 0.3, 0.5)),
       aes(x = fishing_p1, y = rel_biomass, col = model_version)) +
  geom_line(lwd=0.75) +
  facet_wrap(~area_mpa, labeller = labeller(area_mpa = custom_labels)) +  # Use custom labels
  scale_color_manual(values = my_colors, name = "Model version") + 
  labs(x = "Fishing effort", y = "Relative biomass") +
  theme_minimal() +
  theme(text = element_text(size=20), legend.position = "bottom") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), legend.position = "bottom") #save 8x8



#this is figure 5
ggplot(outcome_harvest_long %>%
         filter(area_mpa %in% c(0.05, 0.1, 0.3, 0.5)),
       aes(x = fishing_p1, y = rel_harvest, col = model_version)) +
  geom_line(lwd=0.75) +
  facet_wrap(~area_mpa, labeller = labeller(area_mpa = custom_labels)) +  # Use custom labels
  scale_color_manual(values = my_colors, name = "Model version") + 
  labs(x = "Fishing effort", y = "Relative harvest") +
  theme_minimal() +
  theme(text = element_text(size=20), legend.position = "bottom") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), legend.position = "bottom") #save 8x8






















