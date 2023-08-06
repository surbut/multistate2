## coefficient analysis

## smoothed coefficients
c=coefinput(start = "Ht",stop = "Cad",ages = c(40:80),modelfit = abinom)
m=melt(c,id.vars="age")
ggplotly(ggplot(m,aes(age,exp(value),col=variable,fill=variable))+stat_smooth()+
  geom_hline(yintercept = 1))

coefplotsmoothashr=
function(ages,start,stop,modelfit){
  require(ggplot2)
  agenames=as.character(c(ages))
  s=sapply(agenames,function(x){modelfit$model_list[[x]][[stop]][[start]][,"Estimate"]})
  s=data.frame(s)
  e=sapply(agenames,function(x){modelfit$model_list[[x]][[stop]][[start]][,"Std. Error"]})
  colnames(s)=colnames(e)=agenames
  s=t(s)
  e=t(e)
  melt=melt(s)
  g=ggplot(melt,aes(Var1,value,col=Var2))+geom_point()+stat_smooth() ## does not account for standard error
  
  m=ggplot_build(g)$data[[1]]
  return(list("mat"=m,"plot"=g))}


#RW1 MODEL:
http://faculty.washington.edu/jonno/AV-Smoothing.pdf
φt −φt−1 ∼N(0,σ2), smooth towards the previous value.
RW2 MODEL:
  (φt − φt−1) − (φt−1 − φt−2) ∼ N(0, σ2), smooth towards the previous slope.

## loess inverse variance

coefplotsmoothashr=
  function(ages,start,stop,modelfit){
    require(ggplot2)
    agenames=as.character(c(ages))
    s=sapply(agenames,function(x){modelfit$model_list[[x]][[stop]][[start]][,"Estimate"]})
    s=data.frame(s)
    e=sapply(agenames,function(x){modelfit$model_list[[x]][[stop]][[start]][,"Std. Error"]})
    colnames(s)=colnames(e)=agenames
    s=t(s)
    coeflong=melt(s)
    e=t(e)
    errorlong=melt(e)
    data <- data.frame(coeflong,weights=1/errorlong$value)
   
    
    g=ggplot(data[data$Var2%in%"cad.prs"], aes(x = Var1, y = value, weight = weights,col=Var2,group=Var2)) +
      geom_point() + # Scatter plot of the data
      geom_smooth(method = "loess") + # Smoother
      labs(x = "Time Periods", y = "Coefficients", title = "Smoothed Coefficients Over Time")
    
    gnoweight=ggplot(data, aes(x = Var1, y = value,col=Var2,group=Var2)) +
    geom_point() + # Scatter plot of the data
    geom_smooth(method = "loess") + # Smoother
    labs(x = "Time Periods", y = "Coefficients", title = "Smoothed Coefficients Over Time")
    
    m=ggplot_build(g)$data[[1]]
    return(list("mat"=m,"plot"=g,"noweights"=gnoweight))}



###
# 
# # Sample data
# x <- ages
# y <- coefficients
# weights <- rep(1, length(x))
# 
# # Design matrix
# X <- cbind(1, x)
# 
# # Weighted design matrix and response
# WX <- sqrt(weights) * X
# Wy <- sqrt(weights) * y
# 
# # Solve the normal equations
# beta <- solve(t(WX) %*% WX) %*% t(WX) %*% Wy


calculate_weights <- function(target_age, ages, standard_errors, span) {
  w_error <- 1 / standard_errors
  relative_distance <- abs((ages - target_age) / (span / 2))
  w_distance <- ifelse(relative_distance < 1, (1 - relative_distance^3)^3, 0)
  w_error * w_distance
}

# Example span (e.g., 5 years)
span <- 5

# Example ages and standard errors
ages <- 40:80
standard_errors <- errorlong[,3]

# Calculate weights for target age 20
target_age <- 20
weights <- calculate_weights(target_age, ages, standard_errors, span)

# Use these weights in weighted least squares regression as previously described

rover_loess <- function(ages, coefficients, standard_errors, window_width = 5, degree = 1) {
  n <- length(ages)
  smoothed_coefficients <- numeric(n)
  
  for (i in 1:n) {
    # Find the observations within the window
    in_window <- abs(ages - ages[i]) <= window_width
    
    # Determine distance weights using the tricube weight function
    max_distance <- max(abs(ages[in_window] - ages[i]))
    distance_weights <- (1 - (abs(ages - ages[i]) / max_distance) ^ 3) ^ 3
    distance_weights[!in_window] <- 0
    
    # Multiply by the inverse of the standard errors
    weights <- distance_weights * (1 / standard_errors)
    
    # Create the design matrix
    X <- matrix(1, n, degree + 1)
    for (d in 1:degree) {
      X[, d + 1] <- ages ^ d
    }
    
    # Weighted design matrix and response
    WX <- sqrt(weights) * X
    Wy <- sqrt(weights) * coefficients
    
    # Solve the normal equations
    beta <- solve(t(WX) %*% WX) %*% t(WX) %*% Wy
    
    # Make the prediction
    smoothed_coefficients[i] <- sum(beta * c(1, ages[i] ^ (1:degree)))
  }
  
  return(smoothed_coefficients)
}

# Example usage:
set.seed(123)
ages <- coeflong$Var1[errorlong$Var2%in%"statin_now"]
coefficients <- coeflong$value[errorlong$Var2%in%"statin_now"]
standard_errors <- errorlong$value[errorlong$Var2%in%"statin_now"]


smoothed_coefficients <- custom_loess(ages, coefficients, standard_errors,degree = 2, window_width = 10)

# Plot original and smoothed coefficients
plot(ages, coefficients, main = "Custom LOESS Fit")
lines(ages, smoothed_coefficients, col = "red")
#The degree parameter controls the degree of the polynomial to be fit. Setting degree = 1 fits a linear polynomial, while degree = 2 fits a quadratic polynomial.

The weights combine distance weighting with inverse standard error weighting, so coefficients that are both closer to the target age and more precise will have more influence on the smoothed values.



library(ggplot2)

# Create a data frame to hold the original data and smoothed coefficients
data <- data.frame(
  Age = ages,
  Coefficients = coefficients,
  Custom_LOESS = custom_loess(ages, coefficients, standard_errors, window_width = 10)
)

# Add the standard LOESS fit
data$Standard_LOESS <- predict(loess(Coefficients ~ Age, data, weights = 1 / standard_errors, span = 0.2))
data$Standard_LOESS_noweight =predict(loess(Coefficients ~ Age, data, span = 0.2))
# Create the ggplot
p <- ggplot(data, aes(x = Age, y = Coefficients)) +
  geom_point(aes(color = "Original Data")) +
  geom_line(aes(y = Custom_LOESS, color = "Custom LOESS")) +
  geom_smooth(aes(y = Coefficients, color = "Standard LOESS"), method = "loess", se = FALSE) +
  geom_smooth(aes(y = Coefficients, color = "Standard LOESS_noweight"), method = "loess", se = FALSE) +
  
  labs(title = "Comparison of LOESS Fits",
       y = "Coefficients",
       color = "Legend") +
  theme_minimal()

print(p)


library(ggplot2)
library(tidyr)

# Create a function to apply the custom LOESS to each set of coefficients
apply_custom_loess <- function(coefficients, standard_errors, window_width = 5) {
  ages <- 1:nrow(coefficients)
  smoothed_coefficients <- matrix(0, nrow(coefficients), ncol(coefficients))
  
  for (i in 1:ncol(coefficients)) {
    smoothed_coefficients[, i] <- custom_loess(ages, coefficients[, i], standard_errors[, i], window_width)
  }
  
  return(smoothed_coefficients)
}

# Example coefficients and standard errors (41 ages x 5 sets)
set.seed(123)
coefficients <- s
standard_errors <- e

# Apply custom LOESS smoothing
smoothed_coefficients <- apply_custom_loess(coefficients, standard_errors)

# Create a data frame for ggplot
data <- data.frame(
  Age = rep(1:41, 7),
  Coefficient = as.vector(coefficients),
  Custom_LOESS = as.vector(smoothed_coefficients),
  Coefficient_Set = rep(1:7, each = 41)
)

# Convert to a long format for ggplot
data_long <- pivot_longer(data, c(Coefficient, Custom_LOESS), names_to = "Method", values_to = "Value")

# Create the ggplot
p <- ggplot(data_long, aes(x = Age, y = Value, color = Method)) +
  geom_line() +
  facet_wrap(~Coefficient_Set, ncol = 1, scales = "free_y") +
  labs(title = "Custom LOESS Fit for Different Coefficient Sets",
       y = "Coefficients",
       color = "Method") +
  theme_minimal()

print(p)

