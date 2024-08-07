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
