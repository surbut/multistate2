df_frame <- data.table(...)
# Function to create the filtering expression
create_filter_expr <- function(states, age) {
  s <- sapply(states, function(state) {paste0(age, " < ", state, "_0_censor_age")})
  s2 <- paste0(s, collapse = " & ")
  return(s2)
}
# Function to apply the filtering expression to the data.table
get_at_risk <- function(df, filter_expr) {
  filter_expr_eval <- eval(parse(text = paste0("df[", filter_expr, "]")))
  return(filter_expr_eval)
}
# Function to get at-risk individuals for each age in the age range
get_at_risk_by_age <- function(df, age_range, states) {
  at_risk_list <- list()
  for (age in age_range) {
    filter_expr <- create_filter_expr(states, age)
    atrisk <- get_at_risk(df, filter_expr)
    at_risk_list[[as.character(age)]] <- atrisk
  }
  return(at_risk_list)
}
# Example of usage
nstates <- c("Cad", "Ht", "HyperLip", "Dm")
age_range <- 1:60 # You can set the age range as needed
at_risk_by_age_list <- get_at_risk_by_age(df_frame, age_range, nstates)
# Accessing the at-risk data.table for a specific age
age_30_at_risk <- at_risk_by_age_list[["30"]]