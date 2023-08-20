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