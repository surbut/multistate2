#||||||||||||||||||||||||||||||||
#
# automate calculation of 30yFRS
#
#||||||||||||||||||||||||||||||||

# please note
# for a single individual, we can calculate either the
# original version
# recalibrated version
# refitted version

 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


fun.frs_30y <- function(dat = dat, # data set containing the variables
                         id = "zz_nr", # ID variable
                         sex = "mysex", # must be a factor with two levels, named "male" and "female"
                         age = "myage", 
                         sysbp = "mysysbp",
                         treat = "mytreat", # blood pressure treatment: must be logical! 
                         smoking = "mysmoking", # must be a logical variable!!
                         diabetes = "mydiabetes", # must be a logical variable!!
                         totalchol = "mytotalchol",
                         hdl = "myhdl"
){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  library(survival)
 source("create_values_30y.R")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #............................
  # input checks
  #............................
  # data set must be there
  if(missing(dat)){
    stop("argument 'dat' is missing. please specify data set to use.")
  }
  
  if(!exists(deparse(substitute(dat)))){
    stop(paste("argument 'dat': data set", deparse(substitute(dat)), "does not exist. please specify valid data set to use."))
  }
  
  if(!is.data.frame(dat)){
    stop(paste("argument dat: data set", deparse(substitute(dat)), "is not a valid R data frame. please specify a valid R data frame to use."))
  }
  
  myVars <- c(sex, age, sysbp, treat, smoking, diabetes, hdl, totalchol)
  
  
  if(any(!myVars %in% colnames(dat))){
    stop(paste0(" variables '", paste0(myVars[!myVars %in% colnames(dat)], collapse = "' and '"), "' not available in ", deparse(substitute(dat)), ". please specify valid column names of 'dat'."))
  }else{
    myVars <- colnames(dat)[match(myVars, colnames(dat))]
  }
  
  
  # there must be no missing values in any of the covariates
  if(sum(is.na(dat[ , myVars])) != 0){
    stop("No missing values allowed in any of the covariates!")
  }
  
  
    # prep : logaritmize
  dat$ln_age <- log(dat[, age])
  dat$ln_sysbp <- log(dat[, sysbp])
  dat$ln_totalchol <- log(dat[, totalchol])
  dat$ln_hdl <- log(dat[, hdl])
  
  dat$tmp_sex <- ifelse(dat[,sex] =="male", 1, 0)
  
  # predicotrs
  the_preds <- dat[, c("tmp_sex", "ln_age", "ln_sysbp", treat, smoking, diabetes,"ln_totalchol", "ln_hdl")]
  
  #::::::::::::::::::::::::::::
  # original ----
  #::::::::::::::::::::::::::::
  message("calculation original version....")
  
  
  # xbeta
  linpred <- t(apply(the_preds, 1, function(x){x*coeffs_cvd_orig}))
  xbeta <- rowSums(linpred)
  
  # dxbeta
  linpredd <- t(apply(the_preds, 1, function(x){x*coeffs_compete_orig}))
  dxbeta <- rowSums(linpredd)
  
  # diff and exp
  W5 <- exp(xbeta-mxbeta_orig)
  X5 <- exp(dxbeta-mdxbeta_orig)
  
  # iterate over individuals
  frs_orig <-vapply(seq_along(W5), function(theid){
    columnF <- col2_orig^W5[theid]
    columnK <- col3_orig^X5[theid]
    
    head_of_columnM <- W5[theid]*(-log(col2_orig[1]))
    columnM <- columnF*columnK*columnG_orig*W5[theid]
    therisk <-sum(c(head_of_columnM, head(columnM, -1)))
    therisk
  }, FUN.VALUE = double(1))
  
  rm(linpred, linpredd, xbeta, dxbeta, W5, X5)
  
 
  #::::::::::::::::::::::::::::
  # recalibrated ----
  #::::::::::::::::::::::::::::
  message("calculation recalibrated version....")
  
  
  # xbeta
  linpred <- t(apply(the_preds, 1, function(x){x*coeffs_cvd_orig}))
  xbeta <- rowSums(linpred)
  
  # dxbeta
  linpredd <- t(apply(the_preds, 1, function(x){x*coeffs_compete_orig}))
  dxbeta <- rowSums(linpredd)
 
  # diff and exp
  W5 <- exp(xbeta-mxbeta_recali)
  X5 <- exp(dxbeta-mdxbeta_recali)
  
  # iterate over individuals
  frs_recali <-vapply(seq_along(W5), function(theid){
    columnF <- col2^W5[theid]
    columnK <- col3^X5[theid]
    
    head_of_columnM <- W5[theid]*(-log(col2[1]))
    columnM <- columnF*columnK*columnG*W5[theid]
    therisk <-sum(c(head_of_columnM, head(columnM, -1)))
    therisk
  }, FUN.VALUE = double(1))
  
  rm(linpred, linpredd, xbeta, dxbeta, W5, X5)
  
  #::::::::::::::::::::::::::::
  # refitted ----
  #::::::::::::::::::::::::::::
  message("calculation refitted version....")
  
  # xbeta
  linpred <- t(apply(the_preds, 1, function(x){x*coeffs_cvd_refit}))
  xbeta <- rowSums(linpred)
  
  # dxbeta
  linpredd <- t(apply(the_preds, 1, function(x){x*coeffs_compete_refit}))
  dxbeta <- rowSums(linpredd)
  
  # diff and exp
  W5 <- exp(xbeta-mxbeta_refit)
  X5 <- exp(dxbeta-mdxbeta_refit)
  
  # iterate over individuals
  frs_refit <- vapply(seq_along(W5), function(theid){
    columnF <- col2^W5[theid]
    columnK <- col3^X5[theid]
    
    head_of_columnM <- W5[theid]*(-log(col2[1]))
    columnM <- columnF*columnK*columnG*W5[theid]
    therisk <-sum(c(head_of_columnM, head(columnM, -1)))
    therisk
  }, FUN.VALUE = double(1))
  
  
  as.data.frame(cbind(id =dat[, id], frs_orig =frs_orig*100, frs_recali=frs_recali*100, frs_refit=frs_refit*100))
  
  # end function
}
