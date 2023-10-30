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

colnames(merged_censor)[-1]=colnames(merged_has)[-1]=stringr::str_split_fixed(names(merged_censor)[-1],pattern = ":",n = 2)[,1]

disease_df <- as.data.frame(merged_has)
age_df <- as.data.frame(merged_censor)


a=apply(disease_df,2,function(x){sum(is.na(x))})

# Calculate the co-occurrence matrix
disease_matrix=disease_df[,a==0]
time_matrix=age_df[,a==0]

disease_mx=as.matrix(disease_matrix[,-1])
censor_mx= (disease_mx + 1) %% 2
age_mx=time_matrix[,-1]

co_occurrence_matrix <- t(disease_mx)%*%(disease_mx)

###
# Convert the age matrix into a bin matrix
age_bins <- matrix(NA, nrow = nrow(age_mx), ncol = ncol(age_mx))
for (i in 1:ncol(age_mx)) {
  age_bins[, i] <- cut(age_mx[, i], breaks = seq(0, 100, by = 10), labels = FALSE)
}


library(data.table)

co_occurrence_by_age <- list()

# Get unique age bins; assuming age_bins is a matrix with the same dimensions as disease_mx
unique_age_bins <- unique(as.vector(age_bins))

# Loop through each unique age bin
for(bin in unique_age_bins) {
  # Create a mask for rows corresponding to the current age bin and uncensored data
  mask <- (age_bins == bin)

  # Subset the disease matrix based on the mask
  # Create a matrix of the same dimensions filled with NA
  disease_subset <- matrix(0, nrow=nrow(disease_mx), ncol=ncol(disease_mx))

  # Fill in the values where mask is TRUE
  disease_subset[mask] <- disease_mx[mask]

  # Compute the co-occurrence matrix using fast matrix multiplication
  co_occurrence_matrix <- t(disease_subset) %*% disease_subset

  # Add co-occurrence matrix to list
  co_occurrence_by_age[[as.character(bin)]] <- co_occurrence_matrix
}

  saveRDS(co_occurrence_by_age,file = "~/Library/CloudStorage/Dropbox-Personal//cooccurencebyage.rds")
