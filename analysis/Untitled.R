file_paths <- list.files(path = "~/Desktop/output_gp/", pattern = "\\.tab\\.tsv\\.gz$", full.names = TRUE)
merged_has <- data.table(sample_id = integer(0))

for (file_path in file_paths) {
  file_name <- gsub(".tab.tsv.gz$", "", basename(file_path))
  file_data <-fread(file_path,select = c("sample_id", "has_disease"))
  fd = file_data[, c("sample_id", "has_disease")]
  names(fd)[2] = paste0(file_name, ":has_disease")
   
  # Merge data based on patient_id
  if (ncol(merged_has) == 1) {
    merged_has <- fd
  } else {
    merged_has <-
      merge(merged_has,
            fd,
            by =  "sample_id",
            all = TRUE)
  }
}

merged_has[is.na(merged_has)]=0


colnames(merged_has)[-1]=stringr::str_split_fixed(names(merged_has)[-1],pattern = ":",n = 2)[,1]
saveRDS(merged_has,file = "~/Library/CloudStorage/Dropbox-Personal/merged_has.rds")

merged_censor <- data.table(sample_id = integer(0))

for (file_path in file_paths) {
  file_name <- gsub(".tab.tsv.gz$", "", basename(file_path))
  file_data <-fread(file_path,select = c("sample_id","censor_age"))
  fd = file_data[, c("sample_id", "censor_age")]
  
  names(fd)[2] = paste0(file_name, ":censor_age")
  # Merge data based on patient_id
  if (ncol(merged_censor) == 1) {
    merged_censor <- fd
  } else {
    merged_censor <-
      merge(merged_censor,
            fd,
            by =  "sample_id",
            all = TRUE)
  }
}

saveRDS(merged_censor,file="~/Desktop/output_gp/merged_censor.rds")

#####

merged_has=readRDS("~/Desktop/output_gp/merged_has.rds")
merged_censor=readRDS("~/Desktop/output_gp/merged_censor.rds")

colnames(merged_censor)[-1]=stringr::str_split_fixed(names(merged_censor)[-1],pattern = ":",n = 2)[,1]

disease_df <- as.data.frame(merged_has)
age_df <- as.data.frame(merged_censor)


a=apply(disease_df,2,function(x){sum(is.na(x))})

# Calculate the co-occurrence matrix
disease_matrix=disease_df[,a==0]
time_matrix=age_df[,a==0]

disease_mx=as.matrix(disease_matrix[,-1])

time_matrix=time_matrix[,-1]

co_occurrence_matrix <- t(disease_mx)%*%(disease_mx)



# First, you'd calculate the average age of occurrence for each disease.
age_matrix=age_df[,-184]
avg_age <- colMeans(age_matrix, na.rm = TRUE)  # Assuming NA for no occurrence

# Then, you could use these average ages in correlation calculations, comparing across diseases.
age_correlation_matrix <- cor(age_matrix, use = "pairwise.complete.obs")  # Handles missing values
# Set the diagonal to NA because it's just the count of each disease, not a co-occurrence
diag(co_occurrence_matrix) <- NA

diag(co_occurrence_matrix) <- NA

# You might want to visualize the matrix for easier interpretation
library(ggplot2)
library(reshape2)

co_occurrence_melted <- melt(co_occurrence_matrix)
ggplot(data = co_occurrence_melted, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(survival)

# Assuming 'age_matrix' contains ages and 'disease_matrix' contains occurrence information (1 or 0)
# We would create a survival object for each disease like so:

surv_objects <- list()

for (i in 1:ncol(disease_matrix)) {
  # The survival object considers both time to event and censoring information
  surv_objects[[i]] <- Surv(time = age_matrix[, i], event = disease_matrix[, i])
}


install.packages("survival")
library(survival)

# Initialize a matrix to store the concordance results.
concordance_results <- matrix(NA, ncol = length(surv_objects), nrow = length(surv_objects))

# Adjusting the loop that fills in the concordance_results matrix
for (i in 1:10) {
  for (j in i:10) {  # still exploiting the symmetry of the matrix
    surv_object_i <- Surv(time = age_matrix[, i], event = disease_matrix[, i])
    surv_object_j <- Surv(time = age_matrix[, j], event = disease_matrix[, j])
    
    # Assuming surv_object_i and surv_object_j are correctly defined as before,
    # now focus on the result extraction:
    concordance_result <- concordancefit(surv_object_i, surv_object_j)
    
    # The 'concordance_result' might be a list or otherwise complex structure. 
    # We need to make sure we're extracting the actual statistic correctly.
    # This is hypothetical and depends on the actual structure of 'concordance_result':
    concordance_stat <- concordance_result$concordance  # Or another appropriate key
    
    # Now, the crucial part: ensuring we're inserting this scalar value correctly:
    if (!is.null(concordance_stat) && length(concordance_stat) == 1) {
      concordance_results[i, j] <- concordance_stat
      concordance_results[j, i] <- concordance_stat
    } else {
      # Handle the case where the statistic wasn't extracted as expected
      warning(paste("Unexpected result for pair", i, j, "- check the 'concordance_result' structure."))
    }
  }
}


# You might want to visualize this as a heatmap too
library(ggplot2)
library(reshape2)

co_occurrence_melted <- melt(co_occurrence_matrix)
ggplot(data = co_occurrence_melted, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Co-occurrence count")

rownames(disease_df)=disease_df$sample_id
disease_df=disease_df[,-1]
rownames(age_df)=age_df[,"sample_id"]
age_df=age_df[,-1]
disease_df$person_id=rownames(disease_df)
age_df$person_id=rownames(age_df)
disease_long <- disease_df %>% gather(key = "disease", value = "occurrence", -person_id)
age_long <- age_df %>% gather(key = "disease", value = "age", -person_id)
# Merge the datasets
combined_data <- merge(disease_long, age_long, by = c("person_id", "disease"))
# Calculate correlation between age and occurrence for each disease
correlations <- by(data = combined_data, INDICES = combined_data$disease, FUN = function(df) {
cor(df$age, df$occurrence, use = "pairwise.complete.obs")  # Adjust for missing data
})

correlation_vector <- sapply(correlations, c)

# Now, create a data frame with the disease names and their corresponding correlation
correlation_df <- data.frame(
  disease = names(correlation_vector),
  correlation = as.numeric(correlation_vector)
)

sorted_correlation_df <- correlation_df[order(-correlation_df$correlation), ]

library(ggplot2)

ggplot(sorted_correlation_df, aes(x = reorder(disease, correlation), y = correlation)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flips the axes for a horizontal plot
  theme_minimal() +
  labs(title = "Correlation between Age and Disease Occurrence",
       x = "Disease",
       y = "Correlation")
