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
                         id = "id", # ID variable
                         sex = "mysex", # must be a factor with two levels, named "male" and "female"
                         age = "myage", 
                         sbp = "mysbp",
                         treat = "mytreat", # blood pressure treatment: must be logical! 
                         smoking = "mysmoking", # must be a logical variable!!
                         diabetes = "mydiabetes", # must be a logical variable!!
                         totalchol = "mytotalchol",
                         hdl = "myhdl"){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  library(survival)
 #source("~/multistate2/code/frs30_URBUT/create_values_30y.R")
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
  
  myVars <- c(sex, age, sbp, treat, smoking, diabetes, hdl, totalchol)
  
  
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
  dat$ln_sysbp <- log(dat[,sbp])
  dat$ln_totalchol <- log(dat[, totalchol])
  dat$ln_hdl <- log(dat[, hdl])
  
  
  dat$tmp_sex <- ifelse(dat[,sex] =="male", 1, 0)
  #print(head(dat))
  # predicotrs
  the_preds <- dat[, c("tmp_sex", "ln_age", "ln_sysbp", "mytreat", "mysmoking", "mydiabetes","ln_totalchol", "ln_hdl")]
  
  #::::::::::::::::::::::::::::
  # original ----
  #::::::::::::::::::::::::::::
  #message("calculation original version....")
  
  
  # CVD
  coeffs_cvd_orig <- c("mysexmale" = 0.34362,
                       "ln_age" = 2.63588,
                       "ln_sysbp" = 1.8803,
                       "mytreatTRUE" = 0.5232,
                       "mysmokingTRUE" = 0.59397, 
                       "mydiabetesTRUE" = 0.68602,
                       "ln_totalchol" = 1.12673,
                       "ln_hdl" = -0.90941)
  # competing events
  coeffs_compete_orig <- c("mysexmale" = 0.48123,
                           "ln_age" = 3.39222,
                           "ln_sysbp" = 1.39862,
                           "mytreatTRUE" = 0.19035,
                           "mysmokingTRUE" = 0.99858, 
                           "mydiabetesTRUE" = 0.49756,
                           "ln_totalchol" = -0.00439,
                           "ln_hdl" = 0.16081)
  
  
  
  # xbeta
  linpred <- t(apply(the_preds, 1, function(x){x*coeffs_cvd_orig}))
  xbeta <- rowSums(linpred)
  
  # dxbeta
  linpredd <- t(apply(the_preds, 1, function(x){x*coeffs_compete_orig}))
  dxbeta <- rowSums(linpredd)
  
  
  col2=scan(file = "~/multistate2/code/frs30_URBUT/col2_recali.csv")
  col3=scan(file = "~/multistate2/code/frs30_URBUT/col3_recali.csv")
  
  # log (i-1)-log(i)
  auxi<- c(col2[-1],NA)
  columnG <- log(col2)-log(auxi)
  
  columnG_orig <- scan(file = "~/multistate2/code/frs30_URBUT/columnG_orig.csv")
  #write.csv(columnG, file = "columnG_recali.csv", row.names = FALSE)
  
  
  # these files have been manually created by copying from the excelsheet
  col2_orig <- scan(file = "~/multistate2/code/frs30_URBUT/col2_orig.csv")
  col3_orig <- scan(file = "~/multistate2/code/frs30_URBUT/col3_orig.csv")

  
  #::::::::::::::::::::::::::::
  #  mxbeta and mdxbeta ----
  #::::::::::::::::::::::::::::
  # taken form excel sheet
  mxbeta_orig  <- 21.29326612
  mdxbeta_orig  <- 20.12840698
  
  
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
  

  as.data.frame(cbind(id =dat[,id], frs_orig =frs_orig*100))
  

  
  
  #::::::::::::::::::::::::::::
  # recalibrated ----
  #::::::::::::::::::::::::::::
  #message("calculation recalibrated version....")
  
  
  # xbeta
  linpred <- t(apply(the_preds, 1, function(x){x*coeffs_cvd_orig}))
  xbeta <- rowSums(linpred)
  
  # dxbeta
  linpredd <- t(apply(the_preds, 1, function(x){x*coeffs_compete_orig}))
  dxbeta <- rowSums(linpredd)
  
  mxbeta_recali <- 22.31
 
  mdxbeta_recali <- 21.13

  
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
  
  as.data.frame(cbind(id =dat[, id], frs_orig =frs_orig*100, frs_recali=frs_recali*100))
  
  #
}

