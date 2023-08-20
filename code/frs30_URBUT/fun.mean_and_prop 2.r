#|||||||||||||||||||||||||||||||||||||||||
#
# some magic function
#
# for adaptation of the FRS 1998, calulate the mean values and proprtions of the rsik factors in the dataset at hand
#
#|||||||||||||||||||||||||||||||||||||||||



fun.mean_and_prop <- function(dat = dat,
                              myVars = c("myage", "bp_cat", "mysmoking", "mydiabetes", "totalchol_cat", "ldl_cat", "hdl_cat")
                              ){
  
  
  # metric variables, 
  metricVars <- myVars[vapply(myVars, function(x){is.numeric(dat[ , x])}, FUN.VALUE = logical(1))]
  
  # binary variables: factor variables with two levels or logical
  binVars <- myVars[vapply(myVars, function(x){(is.factor(dat[ , x]) & nlevels(dat[ , x]) == 2L) | is.logical(dat[, x]) }, FUN.VALUE = logical(1))]
  
  # categorical variables factor variables with more levels
  catVars <- myVars[vapply(myVars, function(x){is.factor(dat[ , x]) & nlevels(dat[ , x]) > 2L}, FUN.VALUE = logical(1))]
  
  
  
  #-------------------------
  # metric vars
  #-------------------------
  
  metric_list <- vector(length = length(metricVars), mode = "list")
  names(metric_list) <- metricVars
  
  if(length(metricVars) != 0){
  for(imetricVar in metricVars){
    metric_list[[imetricVar]] <- mean(dat[, imetricVar]) # no na.rm because there shouldn't be any missing values!!
    }
  # end if metric != 0 
  }
  
  #-------------------------
  # binary vars
  #-------------------------
  
  bin_list <- vector(length = length(binVars), mode = "list")
  names(bin_list) <- binVars
  
  if(length(binVars) != 0){
     
    for(ibinVar in binVars){
       if(is.logical(dat[, ibinVar])){
          tmp <- factor(dat[, ibinVar], levels= c(FALSE, TRUE), labels = c("FALSE", "TRUE")) 
          bin_list[[ibinVar]] <- table(tmp)/sum(table(tmp))
       }else{
      bin_list[[ibinVar]] <- table(dat[, ibinVar])/sum(table(dat[, ibinVar]))
       }
    }
  # end if length binVars
  }
  
  #-------------------------
  # catvars
  #-------------------------

  cat_list <- vector(length = length(catVars), mode = "list")
  names(cat_list) <- catVars
    
  if(length(catVars) != 0){
      for(icatVar in catVars){
      cat_list[[icatVar]] <- table(dat[, icatVar])/sum(table(dat[, icatVar]))
    }
  # end if length catVars
  }
  
  return(c(metric_list, bin_list, cat_list))
  
 # end function 
}