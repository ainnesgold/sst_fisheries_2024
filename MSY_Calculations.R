##MSY calculations

library(tidyverse)
library(ggpubr)

projections <- read.csv("anomaly_df.csv")
SST_dev <- projections[['anomaly']]

patch_area <- c(1)
number_patches <- length(patch_area)
timesteps <- 70 #length(SST_dev)
r <- 0.3
K <- 101.3
catchability <- 1

population <- array(NA, dim = c(timesteps, number_patches))
harvest <- array(NA, dim = c(timesteps, number_patches))

patch_area_sequences <- seq(1, 1, by = 1)
fishing_effort_sequences <- seq(0, 1, by = 0.01)

parameter_grid <- expand.grid(patch_area = patch_area_sequences,
                              fishing_effort = fishing_effort_sequences)

############################### BASELINE MODEL ###############################
# Saving outputs
outcome_population <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_population_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))
outcome_harvest_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))



for (iter in 1:nrow(parameter_grid)) {
  population[1] <- c(10)  # Initialize population for each iteration
  
  for (t in 2:timesteps) {
    harvest[t] <- population[t-1] * (1 - exp(-parameter_grid[['fishing_effort']][[iter]] * catchability))
    population[t] <- population[t-1] + r * population[t-1] * (1 - population[t-1] / K) - harvest[t]
  }
  
  #save the whole trajectory
  outcome_population[iter, ] <- population
  outcome_harvest[iter, ] <- harvest
  
  #save the mean of the last 20 years
  outcome_population_end[iter] <- mean(population[(t-19):t])
  outcome_harvest_end[iter] <- mean(harvest[(t-19):t])
}

##TRAJECTORY PLOTS
# Transpose outcome_population for correct plotting
outcome_population <- t(outcome_population)
outcome_harvest <- t(outcome_harvest)

# Plotting with modified line characteristics and legend position
matplot(outcome_population, type = "l", col = 1:nrow(parameter_grid), lty = 1, lwd = 2, xlab = "Time", ylab = "Population")
legend("right", legend = parameter_grid$fishing_effort, col = 1:nrow(parameter_grid), lty = rep(1, nrow(parameter_grid)), ncol=2, lwd = 2, cex = 0.8)


matplot(outcome_harvest, type = "l", col = 1:nrow(parameter_grid), lty = 1, lwd = 2, xlab = "Time", ylab = "Harvest")
legend("right", legend = parameter_grid$fishing_effort, col = 1:nrow(parameter_grid), lty = rep(1, nrow(parameter_grid)), ncol=2, lwd = 2, cex = 0.8)


##PLOTTING EQUILIBRIUM VALUES VS FISHING EFFORT
outcome <- cbind(parameter_grid, outcome_population_end, outcome_harvest_end)

p1<-ggplot(outcome, aes(x=fishing_effort, y=outcome_population_end)) +
  geom_line(lwd=1) +
  labs(x="Fishing effort", y = bquote("Biomass"~(g/m^2))) +
  theme_minimal() +
  theme(text = element_text(size=20),
        legend.position = "bottom")


p2<-ggplot(outcome, aes(x=fishing_effort, y=outcome_harvest_end)) +
  geom_line(lwd=1) +
  labs(x="Fishing effort", y = bquote("Harvest"~(g/m^2))) +
  theme_minimal() +
  theme(text = element_text(size=20),
        legend.position = "bottom")  



figure<-ggarrange(p1+rremove("xlab"), p2+rremove("xlab"),
                  nrow=2, ncol=1, common.legend = TRUE, legend = "top")

annotate_figure(figure, bottom = text_grob("Fishing effort",
                                           size = 20))

#calculate max harvest for optimal nonspatial mgt scenario
# Step 1: Identify the maximum value of outcome$outcome_harvest_end
max_harvest <- max(outcome$outcome_harvest_end)
max_harvest

# Step 2: Find the index (or indices) of the maximum value
max_harvest_index <- which(outcome$outcome_harvest_end == max_harvest)

# Step 3: Extract the corresponding fishing effort value(s)
fishing_effort_at_max_harvest <- outcome$fishing_effort[max_harvest_index]
fishing_effort_at_max_harvest

# Step 4: Find the corresponding population (outcome_population_end) at these fishing effort value(s)
population_at_max_harvest <- outcome$outcome_population_end[max_harvest_index]

# Display the result
population_at_max_harvest



############################### R1 MODEL ###############################
# Saving outputs
outcome_population <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_population_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))
outcome_harvest_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))

population <- array(NA, dim = c(timesteps, number_patches))
harvest <- array(NA, dim = c(timesteps, number_patches))

r_temp_1 <- array(NA, dim = c(timesteps, 1))


for (iter in 1:nrow(parameter_grid)) {
  population[1] <- c(10)  # Initialize population for each iteration
  
  for (t in 2:timesteps) {
    r_temp_1[t] <- 0.3 + 0*SST_dev[t] + -0.0037*SST_dev[t]^2
    harvest[t] <- population[t-1] * (1 - exp(-parameter_grid[['fishing_effort']][[iter]] * catchability))
    population[t] <- population[t] <- population[t-1] + r_temp_1[t] * population[t-1] * 
      (1 - population[t-1] / K) - harvest[t]
  }
  
  #save the whole trajectory
  outcome_population[iter, ] <- population
  outcome_harvest[iter, ] <- harvest
  
  #save the mean of the last 20 years
  outcome_population_end[iter] <- mean(population[(t-19):t], na.rm = TRUE)
  outcome_harvest_end[iter] <- mean(harvest[(t-19):t], na.rm = TRUE)
}

outcome_population <- t(outcome_population)
outcome_harvest <- t(outcome_harvest)
outcome <- cbind(parameter_grid, outcome_population_end, outcome_harvest_end)

#calculate max harvest for optimal nonspatial mgt scenario
# Step 1: Identify the maximum value of outcome$outcome_harvest_end
max_harvest_r1 <- max(outcome$outcome_harvest_end)
max_harvest_r1

# Step 2: Find the index (or indices) of the maximum value
max_harvest_index <- which(outcome$outcome_harvest_end == max_harvest_r1)

# Step 3: Extract the corresponding fishing effort value(s)
fishing_effort_at_max_harvest_r1 <- outcome$fishing_effort[max_harvest_index]
fishing_effort_at_max_harvest_r1

# Step 4: Find the corresponding population (outcome_population_end) at these fishing effort value(s)
population_at_max_harvest_r1 <- outcome$outcome_population_end[max_harvest_index]
population_at_max_harvest_r1




############################### R2 MODEL ###############################
# Saving outputs
outcome_population <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_population_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))
outcome_harvest_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))

r_temp_2 <- array(NA, dim = c(timesteps, 1))
population <- array(NA, dim = c(timesteps, number_patches))
harvest <- array(NA, dim = c(timesteps, number_patches))


for (iter in 1:nrow(parameter_grid)) {
  population[1] <- c(10)  # Initialize population for each iteration
  
  for (t in 2:timesteps) {
    r_temp_2[t] <- 0.3 + 0*(SST_dev[t] - 1) + -0.0037*(SST_dev[t] - 1 )^2
    harvest[t] <- population[t-1] * (1 - exp(-parameter_grid[['fishing_effort']][[iter]] * catchability))
    population[t] <- population[t] <- population[t-1] + r_temp_2[t] * population[t-1] * 
      (1 - population[t-1] / K) - harvest[t]
  }
  
  #save the whole trajectory
  outcome_population[iter, ] <- population
  outcome_harvest[iter, ] <- harvest
  
  #save the mean of the last 20 years
  outcome_population_end[iter] <- mean(population[(t-19):t], na.rm = TRUE)
  outcome_harvest_end[iter] <- mean(harvest[(t-19):t], na.rm = TRUE)
}

outcome_population <- t(outcome_population)
outcome_harvest <- t(outcome_harvest)
outcome <- cbind(parameter_grid, outcome_population_end, outcome_harvest_end)

#calculate max harvest for optimal nonspatial mgt scenario
# Step 1: Identify the maximum value of outcome$outcome_harvest_end
max_harvest_r2 <- max(outcome$outcome_harvest_end)
max_harvest_r2

# Step 2: Find the index (or indices) of the maximum value
max_harvest_index <- which(outcome$outcome_harvest_end == max_harvest_r2)

# Step 3: Extract the corresponding fishing effort value(s)
fishing_effort_at_max_harvest_r2 <- outcome$fishing_effort[max_harvest_index]
fishing_effort_at_max_harvest_r2

# Step 4: Find the corresponding population (outcome_population_end) at these fishing effort value(s)
population_at_max_harvest_r2 <- outcome$outcome_population_end[max_harvest_index]

# Display the result
population_at_max_harvest_r2



############################### R3 MODEL ###############################
# Saving outputs
outcome_population <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_population_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))
outcome_harvest_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))

r_temp_3 <- array(NA, dim = c(timesteps, 1))
population <- array(NA, dim = c(timesteps, number_patches))
harvest <- array(NA, dim = c(timesteps, number_patches))


for (iter in 1:nrow(parameter_grid)) {
  population[1] <- c(10)  # Initialize population for each iteration
  
  for (t in 2:timesteps) {
    r_temp_3[t] <- 0.3 + 0*(SST_dev[t] + 1) + -0.0037*(SST_dev[t] + 1 )^2
    harvest[t] <- population[t-1] * (1 - exp(-parameter_grid[['fishing_effort']][[iter]] * catchability))
    population[t] <- population[t] <- population[t-1] + r_temp_3[t] * population[t-1] * 
      (1 - population[t-1] / K) - harvest[t]
  }
  
  #save the whole trajectory
  outcome_population[iter, ] <- population
  outcome_harvest[iter, ] <- harvest
  
  #save the mean of the last 20 years
  outcome_population_end[iter] <- mean(population[(t-19):t], na.rm = TRUE)
  outcome_harvest_end[iter] <- mean(harvest[(t-19):t], na.rm = TRUE)
}

outcome_population <- t(outcome_population)
outcome_harvest <- t(outcome_harvest)
outcome <- cbind(parameter_grid, outcome_population_end, outcome_harvest_end)

#calculate max harvest for optimal nonspatial mgt scenario
# Step 1: Identify the maximum value of outcome$outcome_harvest_end
max_harvest_r3 <- max(outcome$outcome_harvest_end)
max_harvest_r3

# Step 2: Find the index (or indices) of the maximum value
max_harvest_index <- which(outcome$outcome_harvest_end == max_harvest_r3)

# Step 3: Extract the corresponding fishing effort value(s)
fishing_effort_at_max_harvest_r3 <- outcome$fishing_effort[max_harvest_index]
fishing_effort_at_max_harvest_r3

# Step 4: Find the corresponding population (outcome_population_end) at these fishing effort value(s)
population_at_max_harvest_r3 <- outcome$outcome_population_end[max_harvest_index]

# Display the result
population_at_max_harvest_r3


############################### K1 MODEL ###############################
# Saving outputs
outcome_population <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_population_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))
outcome_harvest_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))

K_temp_1 <- array(NA, dim = c(timesteps, 1))
population <- array(NA, dim = c(timesteps, number_patches))
harvest <- array(NA, dim = c(timesteps, number_patches))


for (iter in 1:nrow(parameter_grid)) {
  population[1] <- c(10)  # Initialize population for each iteration
  
  for (t in 2:timesteps) {
    K_temp_1[t] <- -4.95243768 * SST_dev[t-1] + 101.3
    if (K_temp_1[t] < 10) {
      K_temp_1[t] <- 10
    } else if (K_temp_1[t] > 101.3) {
      K_temp_1[t] <- 101.3
    }
    
    harvest[t] <- population[t-1] * (1 - exp(-parameter_grid[['fishing_effort']][[iter]] * catchability))
    population[t] <- population[t] <- population[t-1] + r * population[t-1] * 
      (1 - population[t-1] / K_temp_1[t]) - harvest[t]
  }
  
  #save the whole trajectory
  outcome_population[iter, ] <- population
  outcome_harvest[iter, ] <- harvest
  
  #save the mean of the last 20 years
  outcome_population_end[iter] <- mean(population[(t-19):t], na.rm = TRUE)
  outcome_harvest_end[iter] <- mean(harvest[(t-19):t], na.rm = TRUE)
}

outcome_population <- t(outcome_population)
outcome_harvest <- t(outcome_harvest)
outcome <- cbind(parameter_grid, outcome_population_end, outcome_harvest_end)

#calculate max harvest for optimal nonspatial mgt scenario
# Step 1: Identify the maximum value of outcome$outcome_harvest_end
max_harvest_K1 <- max(outcome$outcome_harvest_end)
max_harvest_K1

# Step 2: Find the index (or indices) of the maximum value
max_harvest_index <- which(outcome$outcome_harvest_end == max_harvest_K1)

# Step 3: Extract the corresponding fishing effort value(s)
fishing_effort_at_max_harvest_K1 <- outcome$fishing_effort[max_harvest_index]
fishing_effort_at_max_harvest_K1

# Step 4: Find the corresponding population (outcome_population_end) at these fishing effort value(s)
population_at_max_harvest_K1 <- outcome$outcome_population_end[max_harvest_index]

# Display the result
population_at_max_harvest_K1




############################### K2 MODEL ###############################
# Saving outputs
outcome_population <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = timesteps, nrow = nrow(parameter_grid))
outcome_population_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))
outcome_harvest_end <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))

K_temp_2 <- array(NA, dim = c(timesteps, 1))
population <- array(NA, dim = c(timesteps, number_patches))
harvest <- array(NA, dim = c(timesteps, number_patches))


for (iter in 1:nrow(parameter_grid)) {
  population[1] <- c(10)  # Initialize population for each iteration
  
  for (t in 2:timesteps) {
    K_temp_2[t] <- 101.3 + 0*SST_dev[t-1] + -0.7*SST_dev[t-1]^2
    if (K_temp_2[t] < 10) {
      K_temp_2[t] <- 10 }
    else if (K_temp_2[t] > 101.3) {
      K_temp_2[t] <- 101.3
    }
    
    harvest[t] <- population[t-1] * (1 - exp(-parameter_grid[['fishing_effort']][[iter]] * catchability))
    population[t] <- population[t] <- population[t-1] + r * population[t-1] * 
      (1 - population[t-1] / K_temp_2[t]) - harvest[t]
  }
  
  #save the whole trajectory
  outcome_population[iter, ] <- population
  outcome_harvest[iter, ] <- harvest
  
  #save the mean of the last 20 years
  outcome_population_end[iter] <- mean(population[(t-19):t], na.rm = TRUE)
  outcome_harvest_end[iter] <- mean(harvest[(t-19):t], na.rm = TRUE)
}

outcome_population <- t(outcome_population)
outcome_harvest <- t(outcome_harvest)
outcome <- cbind(parameter_grid, outcome_population_end, outcome_harvest_end)

#calculate max harvest for optimal nonspatial mgt scenario
# Step 1: Identify the maximum value of outcome$outcome_harvest_end
max_harvest_K2 <- max(outcome$outcome_harvest_end)
max_harvest_K2

# Step 2: Find the index (or indices) of the maximum value
max_harvest_index <- which(outcome$outcome_harvest_end == max_harvest_K2)

# Step 3: Extract the corresponding fishing effort value(s)
fishing_effort_at_max_harvest_K2 <- outcome$fishing_effort[max_harvest_index]
fishing_effort_at_max_harvest_K2

# Step 4: Find the corresponding population (outcome_population_end) at these fishing effort value(s)
population_at_max_harvest_K2 <- outcome$outcome_population_end[max_harvest_index]

# Display the result
population_at_max_harvest_K2




#Make a dataframe of all the model versions optimal nonspatial management harvest and biomass values

optimal_df <- data.frame(
  model_version = c("baseline", "r1", "r2", "r3", "K1", "K2"),
  fishing_effort = c(fishing_effort_at_max_harvest, fishing_effort_at_max_harvest_r1, fishing_effort_at_max_harvest_r2, fishing_effort_at_max_harvest_r3,
                     fishing_effort_at_max_harvest_K1, fishing_effort_at_max_harvest_K2),
  max_harvest = c(max_harvest, max_harvest_r1, max_harvest_r2, max_harvest_r3, max_harvest_K1, max_harvest_K2),
  population_at_max_harvest = c(population_at_max_harvest, population_at_max_harvest_r1, population_at_max_harvest_r2, population_at_max_harvest_r3, population_at_max_harvest_K1, population_at_max_harvest_K2)
)


