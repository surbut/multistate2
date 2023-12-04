# how to return the array that has the score for every person


pm=readRDS("../output/predictedrsiskboot.rds")
s=statusarray(df_frame = data.table(test),ages = ages,nstates = nstates)
s2=s[c(1:40),,c(1,2,3,4,7,8,9,10)]

get_person_cross_threshold <- function(pm, s2, person_idx) {
  # Convert the risk array for the given person to a data frame
  risk_df <- as.data.frame(as.table(pm[, person_idx, , ]))
  colnames(risk_df) <- c("Bootstrap", "Age", "State", "Risk")
  
  # Convert the state indicator array for the given person to a data frame
  state_df <- as.data.frame(as.table(s2[, person_idx, ]))
  colnames(state_df) <- c("Age", "State", "Indicator")
  
  # Extract the state info for the given person
  current_state_df <- state_df %>%
    filter(Indicator == 1) %>%
    select(-Indicator)
  
  # Merge the data to get risks based on the state the person is in
  merged_df <- risk_df %>%
    inner_join(current_state_df, by = c("Age", "State"))
  
  # Compute mean and SE for the person across ages
  summary_df <- merged_df %>%
    group_by(Age) %>%
    summarise(
      mean_risk = mean(Risk),
      se_risk = sd(Risk) / sqrt(n()),
      #lower_95 = mean_risk - 1.96 * se_risk,  # For 95% CI
      #upper_95 = mean_risk + 1.96 * se_risk,  # For 95% CI
      lower_95 = quantile(Risk,0.025),  # For 95% CI
      upper_95 = quantile(Risk,0.975),  # For 95% CI
      num_exceeded=sum(Risk>thresh),
      
      .groups = "drop"
    )
  
  return(list(mean_risk=summary_df$mean_risk,se=summary_df$se_risk,num_exceeded=summary_df$num_exceeded))}

states=sapply(c(1:nrow(test)),function(x){get_person_cross_threshold(pm,s2,as.character(test$identifier[x]))$mean_risk})

saveRDS(states,"~/multistate2/output/state_occupancy_risk.rds")
