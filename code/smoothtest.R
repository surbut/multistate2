custom_loess <- function(ages, coefficients, standard_errors, window_width, degree = 1) {
  n <- length(ages)
  smoothed_coefficients <- numeric(n)
  inverse_variance <- 1 / (standard_errors^2)
  
  for (i in 1:n) {
    # Determine the neighboring ages for smoothing
    distances <- abs(ages - ages[i])
    
    # Adjust the window width at the extremes
    if (ages[i] <= min(ages) + 5 || ages[i] >= max(ages) - 5) {
      adjusted_window_width <- window_width * 2
    } else {
      adjusted_window_width <- window_width
    }
    
    # Determine the neighbors within the window
    neighbors <- which(distances <= adjusted_window_width)
    
    # Tricube weight function within the window
    tricube_weights <- (1 - (distances[neighbors] / adjusted_window_width) ^ 3) ^ 3
    tricube_weights[distances[neighbors] > adjusted_window_width] <- 0
    
    # Final weights
    weights <- tricube_weights * inverse_variance[neighbors]
    
    # Create the design matrix
    X <- matrix(1, length(neighbors), degree + 1)
    for (d in 1:degree) {
      X[, d + 1] <- ages[neighbors] ^ d
    }
    
    # Weighted design matrix and response
    WX <- sqrt(weights) * X
    Wy <- sqrt(weights) * coefficients[neighbors]
    
    # Solve the normal equations
    beta <- solve(t(WX) %*% WX) %*% t(WX) %*% Wy
    
    # Make the prediction
    smoothed_coefficients[i] <- sum(beta * c(1, ages[i] ^ (1:degree)))
  }
  
  return(smoothed_coefficients)
}

# Function to apply custom and standard LOESS to each set of coefficients


# Apply smoothing (specify the degree of the polynomial if other than 1)
#smoothing_results <- apply_smoothing(coefficients, standard_errors, window_width = 5, degree = 2)

# Function to apply custom and standard LOESS to each set of coefficients
apply_smoothing <-
  function(ages,
           coefficients,
           standard_errors,
           window_width = window_width,
           degree = degree,
           span = span) {
    smoothed_coefficients <-
      matrix(0, nrow(coefficients), ncol(coefficients))
    standard_loess_weighted <-
      matrix(0, nrow(coefficients), ncol(coefficients))
    standard_loess_unweighted <-
      matrix(0, nrow(coefficients), ncol(coefficients))
    
    for (i in 1:ncol(coefficients)) {
      smoothed_coefficients[, i] <-
        custom_loess(ages,
                     coefficients[, i],
                     standard_errors[, i],
                     window_width,
                     degree = degree)
      standard_loess_weighted[, i] <-
        predict(loess(
          coefficients[, i] ~ ages,
          weights = 1 / standard_errors[, i] ^ 2,
          span = span
        ))
      standard_loess_unweighted[, i] <-
        predict(loess(coefficients[, i] ~ ages, span = span))
    }
  
    return(
      list(
        custom_loess = smoothed_coefficients,
        standard_loess_weighted = standard_loess_weighted,
        standard_loess_unweighted = standard_loess_unweighted
      )
    )
  }



plot_coef = function(ages, coefficients, smoothing_results) {
  # Create a data frame for ggplot
  data <- data.frame(
    Age = rep(ages, ncol(coefficients)),
    Coefficient = as.vector(coefficients),
    Custom_LOESS = as.vector(smoothing_results$custom_loess),
    Standard_LOESS_Weighted = as.vector(smoothing_results$standard_loess_weighted),
    Standard_LOESS_Unweighted = as.vector(smoothing_results$standard_loess_unweighted),
    Coefficient_Set = rep(colnames(coefficients), each = length(ages))
    
  )
  
  # Convert to a long format for ggplot
  data_long <-
    pivot_longer(
      data,
      c(
        Coefficient,
        Custom_LOESS,
        Standard_LOESS_Weighted,
        Standard_LOESS_Unweighted
      ),
      names_to = "Method",
      values_to = "Value"
    )
  
  p <- ggplot(data_long, aes(x = Age, y = Value)) +
    geom_point(data = filter(data_long, Method == "Coefficient"), aes(color = Method)) +
    geom_line(data = filter(data_long, Method != "Coefficient"), aes(color = Method)) +
    #geom_smooth(data = filter(data_long, Method == "Coefficient"), aes(color = "red"))+
    facet_wrap( ~ Coefficient_Set, nrow=2, scales = "free_y") +
    labs(title = "Comparison of LOESS Fits for Different Coefficient Sets",
         y = "Coefficients",
         color = "Method") +
    theme_minimal()
}



coefplotsmooth2 = function(ages,
                           start,
                           stop,
                           modelfit,
                           window_width = 10,
                           span = 0.2,
                           degree = 2) {
  require(ggplot2)
  agenames = as.character(c(ages))
  s = sapply(agenames, function(x) {
    modelfit$model_list[[x]][[stop]][[start]][, "Estimate"]
  })
  s = t(s)
  
  e = sapply(agenames, function(x) {
    modelfit$model_list[[x]][[stop]][[start]][, "Std. Error"]
  })
  e = t(e)
  rownames(s) = rownames(e) = agenames
  
  
  smoothing_results <-
    apply_smoothing(
      ages = ages,
      coefficients = s,
      standard_errors = e,
      window_width = window_width,
      degree = degree,
      span = span
    )
  g = plot_coef(
    ages = ages,
    coefficients = s,
    smoothing_results =     smoothing_results
  ) + ggtitle(paste0(start, "_to_", stop, "_Coefficients"))
  w = smoothing_results$standard_loess_weighted
  m = smoothing_results$custom_loess
  colnames(m) = colnames(s)
  rownames(m) = ages
  return(list(
    "custom_smooth" = m,
    "weighted" = w,
    "plot" = g,
    "errors" = e,
    "unsmoothed_coefficients" = s
  ))
}



function(ages,start,stop,modelfit,window_width = 10,span = 0.2,degree=2){
  require(ggplot2)
  agenames=as.character(c(ages))
  s=sapply(agenames,function(x){modelfit$model_list[[x]][[stop]][[start]][,"Estimate"]})
  s=t(s)
  
  e=sapply(agenames,function(x){modelfit$model_list[[x]][[stop]][[start]][,"Std. Error"]})
  e=t(e)
  rownames(s)=rownames(e)=agenames
  
  
  smoothing_results <- apply_smoothing(ages = ages,coefficients = s,standard_errors = e,window_width = window_width,degree = degree,span = span)
  g=plot_coef(ages = ages,coefficients = s,smoothing_results = smoothing_results)+ggtitle(paste0(start,"_to_",stop,"_Coefficients"))
  w=smoothing_results$standard_loess_weighted
  m=smoothing_results$custom_loess
  colnames(m)=colnames(s)
  rownames(m)=ages
  return(list("mat"=m,"weighted"=w,"plot"=g,"errors"=e,"unsmoothed_coefficients"=s))}
