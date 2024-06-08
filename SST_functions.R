#SST functions
#This code creates the Figure 1 SST functions

#growth rate plot
#0.3 commonly used as growth rate
#r = 0 when 9 degree difference - nonlethal critical thermal maxima was 36 degrees
#try this function with different optimal temps (that will be calculated in the anomaly part of sst_cmip6.R)
quad_r <- function (x, a=0.3, b=0, c=-0.0037) {
  a + b*x + c*x^2
}

quad_r2 <- function (x, a=0.3, b=0, c=-0.0037) {
  a + b*(x-1) + c*(x-1)^2
}

quad_r3 <- function (x, a=0.3, b=0, c=-0.0037) {
  a + b*(x+1) + c*(x+1)^2
}

par(mar = c(5, 5, 4, 2) + 0.1)  # Adjust the values as needed

curve(quad_r, from = -10, to = 10, xlab = "Anomaly (°C)", ylab = "Intrinsic growth rate", 
      col = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, lwd = 2)

curve(quad_r, from = -10, to = 10, xlab="Anomaly (°C)", ylab = "", col="#E41A1C", lty = "dashed", lwd = 4)
curve(quad_r2, from = -10, to = 10, add = TRUE, col="#377EB8", lty = "dotted", lwd = 4)
curve(quad_r3, from = -10, to = 10, add = TRUE, col="#4DAF4A", lty = "dotdash", lwd = 4)
title(ylab= "Intrinsic growth rate", line=2, cex.lab=1.2)
legend("topright", legend=c("r1", "r2", "r3"), col=c("#E41A1C", "#377EB8", "#4DAF4A"), lty=c("dashed", "dotted", "dotdash"), cex=1, text.width = 2, lwd = 4)


scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), name = "Model version") +
  



##Piecewise linear function for K
# Define the linear equation
linear_equation <- function(x) {
  -4.95243768 * x +101.3
}

# Create a function to restrict y between 10 and 101.3
#10 g/m2 is maintained so matter how warm it gets (Darling et al. 2017)
#max carrying capacity set to 101.3 g/m2 (MacNeil et al. 2015)
restricted_y <- function(x) {
  y <- linear_equation(x)
  if (y < 10) {
    y <- 10
  } else if (y > 101.3) {
    y <- 101.3
  }
  return(y)
}

# Test the restricted_y function with a value of x
x_value <- 1
result <- restricted_y(x_value)
cat("For x =", x_value, ", y =", result, "\n")

##quadratic - same points used to establish linear
#max is "current" temp and max K
# 3 degree deviation (2040 temp) associated with 6% decline in max
# Define the quad_K function
quad_K <- function(x, a = 101.3, b = 0, c = -0.7) {
  y <- a + b * x + c * x^2
  y <- ifelse(y < 10, 10, y)
  return(y)
}



# Generate x values
x_values <- seq(-20, 20, by = 1)  # Adjust the range and step size as needed

# Apply the restricted_y function to each x value
y_values <- sapply(x_values, restricted_y)

# Plot the function
par(mar = c(5, 5, 4, 2) + 0.1)

plot(x_values, y_values, type = "l", col = "#984EA3", lty="longdash", lwd = 4, cex.lab=1.2,
     xlab = "Anomaly (°C)", ylab = bquote('Carrying capacity'~(g/m^2)))

# Add the quad_K plot on top of the restricted_y plot
# Plot the quad_K function
curve(quad_K, from = -20, to = 20, add=TRUE, col = "#FF7F00", lty = "twodash", lwd = 4)

# Add a legend
legend("topright", legend = c("K1", "K2"), col = c("#984EA3", "#FF7F00"), lty = c("longdash", "twodash"), lwd = 4)


  
