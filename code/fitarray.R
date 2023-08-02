fitfunc = function(df_frame, ages, nstates, mode,covariates) {
  nar = array(
    data = NA,
    dim = c(length(ages), length(nstates), length(nstates)),
    dimnames = list(ages, nstates, nstates)
  )
  event = array(
    data = NA,
    dim = c(length(ages), length(nstates), length(nstates)),
    dimnames = list(ages, nstates, nstates)
  )
  
  mean = array(
    data = NA,
    dim = c(length(ages), length(nstates), length(nstates)),
    dimnames = list(ages, nstates, nstates)
  )
  endlist = startlist = vector("list", length = length(nstates))
  names(endlist) = names(startlist) = nstates
  
  agelist = vector("list", length = length(ages))
  timeatrisk = agelist
  
  names(agelist) = names(timeatrisk) = ages
  for (i in 1:length(agelist)) {
    agelist[[i]] = startlist
    for (j in 1:length(startlist)) {
      agelist[[i]][[j]] = endlist
    }
  }
  
  for (i in 1:length(timeatrisk)) {
    timeatrisk[[i]] = startlist
  }
  
  
  
  
  
  
  df_frame$cad.prs = scale(df_frame$cad.prs)
  
  
  # from Health (1) to Health directly alone
  
  for (i in 1:length(ages)) {
    age = ages[i]
    agename = as.character(age)
    nx = age + 1
    atrisk = df_frame[age < Cad_0_censor_age &
                        age < Ht_0_censor_age &
                        age < HyperLip_0_censor_age &
                        age < Dm_0_censor_age &
                        age < Death_Censor_Age, ]
    atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                 atrisk$statin_age <= nx, 1, 0)
    atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                  atrisk$htn_age <= nx, 1, 0)
    
    NAR = dim(atrisk)[1]

    #atrisk=atrisk%>%mutate(yearsinstate=age) all same years in state
    
    if (nrow(atrisk) < 10) {
      agelist[[agename]][["Cad"]][["Health"]] = NA
      
      agelist[[agename]][["Health"]][["Health"]] = NA
      
      agelist[[agename]][["death"]][["Health"]] = NA
      
      agelist[[agename]][["Ht"]][["Health"]] = NA
      
      
      agelist[[agename]][["HyperLip"]][["Health"]] = NA
      
      agelist[[agename]][["Dm"]][["Health"]] = NA
    } else{
      
      my_list <- covariates
      
      # Remove the word "banana"
      new_list <- gsub("\\+yearsinstate", "", my_list)
      
      ## from Health (1) to Health 
      censored = dim(atrisk[which(
        Cad_0_censor_age > nx &
          nx < Ht_0_censor_age &
          nx < HyperLip_0_censor_age &
          nx < Dm_0_censor_age & 
          nx < Death_Censor_Age
      ), ])[1]
      
      nar[i, "Health", 1] = NAR
      event[i, "Health", 1] = censored
      mean[i, "Health", 1] = censored / NAR
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Cad_0_censor_age > nx &
            nx < Ht_0_censor_age &
            nx < HyperLip_0_censor_age &
            nx < Dm_0_censor_age & nx < Death_Censor_Age,
          1,
          0
        ) ~", new_list)),
        data = atrisk
      )
      agelist[[agename]][["Health"]][["Health"]] =summary(fit2)$coefficients
      
      rm(censored)
      
      ## from Health (1) to CAD directly alone
      
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 #&
                                    #nx < Death_Censor_Age
                                    ), ])[1]
      nar[i, "Cad", 1] = NAR
      event[i, "Cad", 1] = censored
      mean[i, "Cad", 1] = censored / NAR
      
      rm(censored)
      fit2 = glm(
        family = mode,as.formula(paste0("
        ifelse(
          Cad_0_censor_age <= nx &
            Cad_0_Any == 2 & nx < Death_Censor_Age,
          1,
          0
        ) ~", new_list)),
        ## here age represents time in state
        data = atrisk
      )
      agelist[[agename]][["Cad"]][["Health"]] =summary(fit2)$coefficients
      
      ## from Health (1) to Death directly alone (but could die with CAD, make CAD and death marginal, 
      ## i.e., don't make conditions that other conditions don't develop that year)
      
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      nar[i, "death", 1] = NAR
      event[i, "death", 1] = censored
      mean[i, "death", 1] = censored / NAR
      rm(censored)
      fit2 = glm(as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0)~",new_list)),
                 family=mode,data=atrisk)
      agelist[[agename]][["Death"]][["Health"]]=summary(fit2)$coefficients
      #rm(atrisk)
    
    
    
    ### from health to HT

      censored = dim(atrisk[which(
        Ht_0_censor_age <= nx & Ht_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < HyperLip_0_censor_age &
          nx < Dm_0_censor_age & nx < Death_Censor_Age
      ), ])[1]
      nar[i, "Ht", 1] = NAR
      event[i, "Ht", 1] = censored
      mean[i, "Ht", 1] = censored / NAR
      rm(censored)
      
      ### from health to HyperLip
      
      
      censored = dim(atrisk[which(
        HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < Ht_0_censor_age & nx < Dm_0_censor_age &
          nx < Death_Censor_Age
      ), ])[1]
      
      nar[i, "HyperLip", 1] = NAR
      event[i, "HyperLip", 1] = censored
      mean[i, "HyperLip", 1] = censored / NAR
      rm(censored)
      
      
      ## from health to Dm
      
      censored = dim(atrisk[which(
        Dm_0_censor_age <= nx & Dm_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < Ht_0_censor_age & nx < HyperLip_0_censor_age &
          nx < Death_Censor_Age
      ), ])[1]
      
      nar[i, "Dm", 1] = NAR
      event[i, "Dm", 1] = censored
      mean[i, "Dm", 1] = censored / NAR
      rm(censored)
      rm(atrisk)
    }
    ## From one risk states to 2,3,4,CAD or death
    
    # From HT
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Ht_0_censor_age &
                        Ht_0_Any == 2 &
                        age < HyperLip_0_censor_age &
                        age < Dm_0_censor_age, ]
    NAR = dim(atrisk)[1]
    if (nrow(atrisk) < 10) {
      agelist[[agename]][["Cad"]][["Ht"]] = NA
      
      agelist[[agename]][["Ht"]][["Ht"]] = NA
      
      agelist[[agename]][["death"]][["Ht"]] = NA
      
      agelist[[agename]][["Ht&HyperLip"]][["Ht"]] = NA
      
      agelist[[agename]][["Ht&Dm"]][["Ht"]] = NA
      
      agelist[[agename]][["Ht&HyperLip&Dm"]][["Ht"]] = NA
      
      
      
      #se_high[i,"death",2]=NA
    } else{
      atrisk = atrisk %>% mutate(yearsinstate = (age - Ht_0_censor_age))
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)
      
      timeatrisk[[agename]][["Ht"]] = atrisk$yearsinstate
      
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    nx < Death_Censor_Age), ])[1]
      nar[i, "Cad", 2] = NAR
      
      event[i, "Cad", 2] = censored
      mean[i, "Cad", 2] = censored / NAR
      rm(censored)
      
      
      fit2 = glm(
        family = mode,
        ## in order to parse the formulat
        as.formula(paste0("ifelse(
          Cad_0_censor_age <= nx &
            Cad_0_Any == 2 & nx < Death_Censor_Age,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["Cad"]][["Ht"]] = summary(fit2)$coefficients
      
      ## from HT (2) to Ht directly alone
      censored = dim(atrisk[which(
        Ht_0_censor_age <= nx & Ht_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < HyperLip_0_censor_age &
          nx < Dm_0_censor_age & nx < Death_Censor_Age
      ), ])[1]
      nar[i, "Ht", 2] = NAR
      event[i, "Ht", 2] = censored
      mean[i, "Ht", 2] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Ht_0_censor_age <= nx & Ht_0_Any == 2
          &
            Cad_0_censor_age > nx &
            nx < HyperLip_0_censor_age &
            nx < Dm_0_censor_age & nx < Death_Censor_Age,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      
      agelist[[agename]][["Ht"]][["Ht"]] =summary(fit2)$coefficients
      #
      
      ## from HT (2) to Death directly alone (but could die with CAD)
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any  == 2 &
                                    Ht_0_censor_age < nx &
                                    Ht_0_Any == 2), ])[1]
      nar[i, "death", 2] = NAR
      event[i, "death", 2] = censored
      mean[i, "death", 2] = censored / NAR
      rm(censored)
      
      
      fit2 = glm(as.formula(paste0("ifelse(Death_Censor_Age <= nx &Death_censor_Any == 2, 1, 0) ~",covariates)),
        family=mode,data=atrisk)
      agelist[[agename]][["death"]][["Ht"]] =summary(fit2)$coefficients
      
      ## from HT (2) to Ht and Hyperlip directly 
      censored = dim(atrisk[which(HyperLip_0_censor_age <= nx &
                                    HyperLip_0_Any  == 2 &
                                   Cad_0_censor_age > nx &
                                    nx < Dm_0_censor_age), ])[1]
      nar[i, "Ht&HyperLip", 2] = NAR
      event[i, "Ht&HyperLip", 2] = censored
      mean[i, "Ht&HyperLip", 2] = censored / NAR
      
      ## from HT (2) to Ht and Dm directly 
      censored = dim(atrisk[which(Dm_0_censor_age <= nx &
                                    Dm_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < HyperLip_0_censor_age), ])[1]
      nar[i, "Ht&Dm", 2] = NAR
      event[i, "Ht&Dm", 2] = censored
      mean[i, "Ht&Dm", 2] = censored / NAR
      
      ## from HT (2) to Ht and Dm and Hyperip directly 
      censored = dim(atrisk[which(Dm_0_censor_age <= nx &
                                    Dm_0_Any  == 2 &
                                    HyperLip_0_censor_age <= nx &
                                    HyperLip_0_Any  == 2 &
                                    Cad_0_censor_age > nx), ])[1]
      nar[i, "Ht&HyperLip&Dm", 2] = NAR
      event[i, "Ht&HyperLip&Dm", 2] = censored
      mean[i, "Ht&HyperLip&Dm", 2] = censored / NAR
      
      
      rm(atrisk)
      
    }
    
    
    
    ##  from Hyperlip
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2 &
                        age < Ht_0_censor_age &
                        age < Dm_0_censor_age, ]
    NAR = dim(atrisk)[1]
    
    
    if (nrow(atrisk) < 10) {
      agelist[[agename]][["Cad"]][["HyperLip"]] = NA
      
      agelist[[agename]][["HyperLip"]][["HyperLip"]] = NA
      
      agelist[[agename]][["death"]][["HyperLip"]] = NA
      
      agelist[[agename]][["HyperLip&Dm"]][["HyperLip"]] = NA
      
      agelist[[agename]][["Ht&HyperLip"]][["HyperLip"]] = NA
      
      agelist[[agename]][["Ht&HyperLip&Dm"]][["HyperLip"]] = NA
      
    } else{
      ## from Hyperlip (3) to CAD marginal (i.e., could develop other states in the interim year)
      atrisk = atrisk %>% mutate(yearsinstate = age - HyperLip_0_censor_age)
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)
      
      timeatrisk[[agename]][["HyperLip"]] = atrisk$yearsinstate
      
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2
                                  & nx < Death_Censor_Age), ])[1]
      
      
      nar[i, "Cad", 3] = NAR
      event[i, "Cad", 3] = censored
      mean[i, "Cad", 3] = censored / NAR
      
      rm(censored)
      
      fit2 = glm(
        as.formula(paste0("ifelse(Cad_0_censor_age <= nx &Cad_0_Any == 2, 1, 0) ~",
          covariates)),
        data = atrisk,family=mode)
      
      agelist[[agename]][["Cad"]][["HyperLip"]] =summary(fit2)$coefficients
      
      
      ## from Hyperlip (2) to Hyperlip directly alone
      
      
      censored = dim(atrisk[which(
        HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < Ht_0_censor_age & nx < Dm_0_censor_age &
          nx < Death_Censor_Age
      ), ])[1]
      
      nar[i, "HyperLip", 3] = NAR
      event[i, "HyperLip", 3] = censored
      mean[i, "HyperLip", 3] = censored / NAR
      rm(censored)
      
      
      fit2 = glm(
        family = mode,as.formula(paste0("
        ifelse(
          HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
          & Cad_0_censor_age > nx &
            nx < Ht_0_censor_age &
            nx < Dm_0_censor_age &
            nx < Death_Censor_Age,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      
      agelist[[agename]][["HyperLip"]][["HyperLip"]] =summary(fit2)$coefficients
      
      ## from HyperLip (2) to Death directly alone (but could die with CAD)
      
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      nar[i, "death", 3] = NAR
      event[i, "death", 3] = censored
      mean[i, "death", 3] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0)~",covariates)),
        data = atrisk
      )
      
      agelist[[agename]][["death"]][["HyperLip"]] =summary(fit2)$coefficients
      
      
      
      ## from HyperLip (2) to Ht and Hyperlip directly 
      censored = dim(atrisk[which(Ht_0_censor_age <= nx &
                                    Ht_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < Dm_0_censor_age), ])[1]
      nar[i, "Ht&HyperLip", "HyperLip"] = NAR
      event[i, "Ht&HyperLip", "HyperLip"] = censored
      mean[i, "Ht&HyperLip", "HyperLip"] = censored / NAR
      
      ## from Hyperlip (2) to Hyperlip and Dm directly 
      censored = dim(atrisk[which(Dm_0_censor_age <= nx &
                                    Dm_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < Ht_0_censor_age), ])[1]
      nar[i, "HyperLip&Dm", 3] = NAR
      event[i, "HyperLip&Dm", 3] = censored
      mean[i, "HyperLip&Dm", 3] = censored / NAR
      
      ## from HT (2) to Ht and Dm and Hyperip directly 
      censored = dim(atrisk[which(Dm_0_censor_age <= nx &
                                    Dm_0_Any  == 2 &
                                    Ht_0_censor_age <= nx &
                                    Ht_0_Any  == 2 &
                                    Cad_0_censor_age > nx), ])[1]
      nar[i, "Ht&HyperLip&Dm", 3] = NAR
      event[i, "Ht&HyperLip&Dm", 3] = censored
      mean[i, "Ht&HyperLip&Dm", 3] = censored / NAR
      
      
      
      
      rm(atrisk)
      
    }
    
    
    
    ## from Dm2
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age &
                        Dm_0_Any == 2 &
                        age < Ht_0_censor_age &
                        age < HyperLip_0_censor_age, ]
    NAR = dim(atrisk)[1]
    
    
    if (nrow(atrisk) < 10) {
      agelist[[agename]][["Cad"]][["Dm"]] = NA
      
      agelist[[agename]][["Dm"]][["Dm"]] = NA
      
      agelist[[agename]][["death"]][["Dm"]] = NA
      
      agelist[[agename]][["Ht&Dm"]][["Dm"]] = NA
      
      agelist[[agename]][["HyperLip&Dm"]][["Dm"]] = NA
      
      agelist[[agename]][["Ht&HyperLip&Dm"]][["Dm"]] = NA
      
      
    } else{
      # from Dm (4) to CAD directly alone
      
      atrisk = atrisk %>% mutate(yearsinstate = (age - Dm_0_censor_age))
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)
      
      timeatrisk[[agename]][["Dm"]] = atrisk$yearsinstate
      
      censored = dim(atrisk[which(
        Cad_0_censor_age <= nx & Cad_0_Any == 2 &
          Dm_0_censor_age <= nx & Dm_0_Any == 2
        #&nx<HyperLip_0_censor_age&nx<Ht_0_censor_age
        & nx < Death_Censor_Age
      ), ])[1]
      
      nar[i, "Cad", 4] = NAR
      event[i, "Cad", 4] = censored
      mean[i, "Cad", 4] = censored / NAR
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Cad_0_censor_age <= nx & Cad_0_Any == 2 &
            Dm_0_censor_age <= nx & Dm_0_Any == 2
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      
      agelist[[agename]][["Cad"]][["Dm"]] =summary(fit2)$coefficients
      
      
      
      rm(censored)
      
      ## from DM (4) to Dm directly alone
      censored = dim(atrisk[which(
        Dm_0_censor_age <= nx &
          Dm_0_Any == 2 &
          Cad_0_censor_age > nx &
          nx < HyperLip_0_censor_age &
          nx < Ht_0_censor_age & nx < Death_Censor_Age
      ), ])[1]
      
      nar[i, "Dm", 4] = NAR
      event[i, "Dm", 4] = censored
      mean[i, "Dm", 4] = censored / NAR
      agelist[[agename]][["Dm"]][["Dm"]] =summary(fit2)$coefficients
      
      #rm(fit)
      rm(censored)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Dm_0_censor_age <= nx &
            Dm_0_Any == 2 &
            Cad_0_censor_age > nx &
            nx < HyperLip_0_censor_age & nx < Ht_0_censor_age &
            nx < Death_Censor_Age
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      
      agelist[[agename]][["Dm"]][["Dm"]] =summary(fit2)$coefficients
      
      ## from Dm to Death directly alone (but could die with CAD)
      
      censored = dim(atrisk[which(
        Death_Censor_Age <= nx &
          Death_censor_Any == 2), ])[1]
      
      nar[i, "death", 4] = NAR
      event[i, "death", 4] = censored
      mean[i, "death", 4] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["death"]][["Dm"]] =summary(fit2)$coefficients
      
      
      
      
      ## from Dm (4) to Ht and Dm directly 
      censored = dim(atrisk[which(Ht_0_censor_age <= nx &
                                    Ht_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < Dm_0_censor_age), ])[1]
      nar[i, "Ht&Dm", "Dm"] = NAR
      event[i, "Ht&Dm", "Dm"] = censored
      mean[i, "Ht&Dm", "Dm"] = censored / NAR
      
      ## from Dm (4) to Hyperlip and Dm directly 
      censored = dim(atrisk[which(HyperLip_0_censor_age <= nx &
                                    HyperLip_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < Ht_0_censor_age), ])[1]
      nar[i, "HyperLip&Dm",4] = NAR
      event[i, "HyperLip&Dm", 4] = censored
      mean[i, "HyperLip&Dm", 4] = censored / NAR
      
      ## from Dm (2) to Ht and Dm and Hyperip directly 
      censored = dim(atrisk[which(HyperLip_0_censor_age <= nx &
                                    HyperLip_0_Any  == 2 &
                                    Ht_0_censor_age <= nx &
                                    Ht_0_Any  == 2 &
                                    Cad_0_censor_age > nx), ])[1]
      nar[i, "Ht&HyperLip&Dm", 4] = NAR
      event[i, "Ht&HyperLip&Dm", 4] = censored
      mean[i, "Ht&HyperLip&Dm", 4] = censored / NAR
      
      
      
      
      rm(atrisk)
    }
    
    ## from CAD to death (to do, should we add RF plus CAD as a starting RF? (no, because that would add the complexity of having to add RF + CAD as an ending state and right now we consider CAD as an 'absorbing' state i.e., if you are diagnosed with an additional RF in the interin  year between diagnosis, not considered, ending with CAD is the same as ending with CAD + ...))
    ## from CAD marginal to Death directly
    
    atrisk = df_frame[age > Cad_0_censor_age &
                        Cad_0_Any == 2 & Death_Censor_Age > age, ]
    NAR = dim(atrisk)[1]
    
    if (nrow(atrisk) < 10) {
      agelist[[agename]][["Cad"]][["Cad"]] = NA
      
      
      
      agelist[[agename]][["death"]][["Cad"]] = NA
    } else{
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      atrisk = atrisk %>% mutate(yearsinstate = age - Cad_0_censor_age)
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)
      
      timeatrisk[[agename]][["Cad"]] = atrisk$yearsinstate
      
      #fit=glm(family=mode,ifelse( Death_Censor_Age<=nx&Death_censor_Any==2
      # ,1,0)~cad.prs.lev,data = atrisk)
      
      nar[i, "death", 5] = NAR
      event[i, "death", 5] = censored
      mean[i, "death", 5] = censored / NAR
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~",covariates)),
        data = atrisk
      )
      
      
      censored = dim(atrisk[which(Death_Censor_Age > nx &
                                    Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2), ])[1]
      
      nar[i, "Cad", 5] = NAR
      event[i, "Cad", 5] = censored
      mean[i, "Cad", 5] = censored / NAR
      agelist[[agename]][["Cad"]][["Cad"]] =summary(fit2)$coefficients
      
      
      rm(censored)
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Death_Censor_Age > nx &
            Cad_0_censor_age <= nx & Cad_0_Any == 2,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["Cad"]][["Cad"]] =summary(fit2)$coefficients
      
      
      
      
      
      
      
      ## from HyperLip (2) to Ht and Hyperlip directly 
      censored = dim(atrisk[which(Ht_0_censor_age <= nx &
                                    Ht_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < Dm_0_censor_age), ])[1]
      nar[i, "Ht&HyperLip", "HyperLip"] = NAR
      event[i, "Ht&HyperLip", "HyperLip"] = censored
      mean[i, "Ht&HyperLip", "HyperLip"] = censored / NAR
      
      ## from Hyperlip (2) to Hyperlip and Dm directly 
      censored = dim(atrisk[which(Dm_0_censor_age <= nx &
                                    Dm_0_Any  == 2 &
                                    Cad_0_censor_age > nx &
                                    nx < Ht_0_censor_age), ])[1]
      nar[i, "HyperLip&Dm", 3] = NAR
      event[i, "HyperLip&Dm", 3] = censored
      mean[i, "HyperLip&Dm", 3] = censored / NAR
      
      ## from HT (2) to Ht and Dm and Hyperip directly 
      censored = dim(atrisk[which(Dm_0_censor_age <= nx &
                                    Dm_0_Any  == 2 &
                                    Ht_0_censor_age <= nx &
                                    Ht_0_Any  == 2 &
                                    Cad_0_censor_age > nx), ])[1]
      nar[i, "Ht&HyperLip&Dm", 2] = NAR
      event[i, "Ht&HyperLip&Dm", 2] = censored
      mean[i, "Ht&HyperLip&Dm", 2] = censored / NAR
      
      rm(atrisk)
    }
    
    ### 2 risk states
    ## from Ht and DM to CAD marginal
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age & Dm_0_Any == 2 &
                        age > Ht_0_censor_age & Ht_0_Any == 2
                      & age < HyperLip_0_censor_age, ]
    NAR = dim(atrisk)[1]
    
    if (nrow(atrisk) < 10)
    {
      agelist[[agename]][["Cad"]][["Ht&Dm"]] = NA
      
      agelist[[agename]][["Ht&Dm"]][["Ht&Dm"]] = NA
      agelist[[agename]][["death"]][["Ht&Dm"]] = NA
      
      
    } else{
      atrisk$yearsinstate = apply(atrisk[, c("Ht_0_censor_age", "Dm_0_censor_age")], 1, function(x) {
        age - max(x[1], x[2])
      })
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)
      
      
      timeatrisk[[agename]][["Ht&Dm"]] = atrisk$yearsinstate
      
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    nx < Death_Censor_Age), ])[1]
      nar[i, "Cad", 9] = NAR
      event[i, "Cad", 9] = censored
      mean[i, "Cad", 9] = censored / NAR
      rm(censored)
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Death_Censor_Age > nx & Cad_0_censor_age <= nx & Cad_0_Any == 2
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["Cad"]][["Ht&Dm"]] =summary(fit2)$coefficients
      ## from Ht and Dm to Death directly alone (but could die with CAD)
      
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      #fit=glm(family=mode,ifelse( Death_Censor_Age<=nx&Death_censor_Any==2
      #                ,1,0)~cad.prs.lev,data = atrisk)
      nar[i, "death", 9] = NAR
      event[i, "death", 9] = censored
      mean[i, "death", 9] = censored / NAR
      
      
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["death"]][["Ht&Dm"]] =summary(fit2)$coefficients
      
      
      ## from Ht and Dm to Ht and Dm directly
      
      censored = dim(atrisk[which(
        nx < Cad_0_censor_age &
          nx >= Dm_0_censor_age & Dm_0_Any == 2 &
          nx >= Ht_0_censor_age & Ht_0_Any == 2
        & nx < HyperLip_0_censor_age
      ), ])[1]
      
      nar[i, "Ht&Dm", 9] = NAR
      event[i, "Ht&Dm", 9] = censored
      mean[i, "Ht&Dm", 9] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          nx < Cad_0_censor_age &
            nx >= Dm_0_censor_age & Dm_0_Any == 2 &
            nx >= Ht_0_censor_age & Ht_0_Any == 2
          & nx < HyperLip_0_censor_age
          ,
          1,
          0
        ) ~", covariates)),
        data = atrisk
      )
      
      agelist[[agename]][["Ht&Dm"]][["Ht&Dm"]] =summary(fit2)$coefficients
      rm(atrisk)
      
    }
    
    
    
    ## from HyperLip and DM to CAD marginal
    ## from HyperLip and Dm to Death directly alone (but could die with CAD) ###### BAD
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age &
                        Dm_0_Any == 2 &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2 &
                        age < Ht_0_censor_age, ]
    NAR = dim(atrisk)[1]
    
    if (nrow(atrisk) < 10)
    {
      agelist[[agename]][["Cad"]][["HyperLip&Dm"]] = NA
      
      agelist[[agename]][["HyperLip&Dm"]][["HyperLip&Dm"]] = NA
      agelist[[agename]][["death"]][["HyperLip&Dm"]] = NA
      
      
    } else{
      atrisk$yearsinstate = apply(atrisk[, c("HyperLip_0_censor_age", "Dm_0_censor_age")], 1, function(x) {
        age - max(x[1], x[2])
      })
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)
      
      timeatrisk[[agename]][["HyperLip&Dm"]] = atrisk$yearsinstate
      
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      nar[i, "death", 8] = NAR
      event[i, "death", 8] = censored
      mean[i, "death", 8] = censored / NAR
      
      rm(censored)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["death"]][["HyperLip&Dm"]] =summary(fit2)$coefficients
      
      
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    Death_Censor_Age > nx), ])[1]
      nar[i, "Cad", 8] = NAR
      event[i, "Cad", 8] = censored
      mean[i, "Cad", 8] = censored / NAR
      
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Cad_0_censor_age <= nx & Cad_0_Any == 2 & Death_Censor_Age > nx
          ,
          1,
          0
        ) ~",covariates)),data = atrisk
      )
      agelist[[agename]][["Cad"]][["HyperLip&Dm"]] =summary(fit2)$coefficients
      
      
      censored = dim(atrisk[which(
        nx < Cad_0_censor_age &
          nx > Dm_0_censor_age &
          Dm_0_Any == 2 &
          nx > HyperLip_0_censor_age &
          HyperLip_0_Any == 2 & nx < Ht_0_censor_age
      ), ])[1]
      nar[i, "HyperLip&Dm", 8] = NAR
      event[i, "HyperLip&Dm", 8] = censored
      mean[i, "HyperLip&Dm", 8] = censored / NAR
      
      #rm(fit)
      
      rm(censored)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          nx < Cad_0_censor_age &
            nx > Dm_0_censor_age &
            Dm_0_Any == 2 &
            nx > HyperLip_0_censor_age &
            HyperLip_0_Any == 2 & nx < Ht_0_censor_age
          ,
          1,
          0
        ) ~", covariates)),
        data = atrisk
      )
      agelist[[agename]][["HyperLip&Dm"]][["HyperLip&Dm"]] =summary(fit2)$coefficients
      rm(atrisk)
      
    }
    
    
    
    ## from HyperLip and Ht to CAD marginal
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Ht_0_censor_age &
                        Ht_0_Any == 2 &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2
                      & age < Dm_0_censor_age, ]
    NAR = dim(atrisk)[1]
    
    
    if (nrow(atrisk) < 10)
    {
      agelist[[agename]][["Cad"]][["Ht&HyperLip"]] = NA
      
      agelist[[agename]][["Ht&HyperLip"]][["Ht&HyperLip"]] = NA
      agelist[[agename]][["death"]][["Ht&HyperLip"]] = NA
      
    } else{
      atrisk$yearsinstate = apply(atrisk[, c("HyperLip_0_censor_age", "Ht_0_censor_age")], 1, function(x) {
        age - max(x[1], x[2])
      })
      timeatrisk[[agename]][["Ht&HyperLip"]] = atrisk$yearsinstate
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)
      
      
      censored = dim(atrisk[which(
        Cad_0_censor_age <= nx &
          Cad_0_Any == 2 &
          Ht_0_censor_age <= nx &
          Ht_0_Any == 2 &
          HyperLip_0_censor_age <= nx &
          HyperLip_0_Any == 2 & nx < Death_Censor_Age
      ), ])[1]
      
      
      
      
      nar[i, "Cad", 7] = NAR
      event[i, "Cad", 7] = censored
      mean[i, "Cad", 7] = censored / NAR
      
      
      rm(censored)
      ##rm(fit)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Cad_0_censor_age <= nx &
            Cad_0_Any == 2 &
            Ht_0_censor_age <= nx &
            Ht_0_Any == 2 &
            HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2 &
            nx < Death_Censor_Age,1,0) ~",covariates)),
        data = atrisk
      )
      
      
      agelist[[agename]][["Cad"]][["Ht&HyperLip"]] =summary(fit2)$coefficients
      
      
      ## from  HyperLip and Ht to HyperLip and Ht
      
      censored = dim(atrisk[which(
        Ht_0_censor_age <= nx &
          Ht_0_Any == 2 &
          HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < Dm_0_censor_age & nx < Death_Censor_Age
      ), ])[1]
      
      
      
      nar[i, "Ht&HyperLip", 7] = NAR
      event[i, "Ht&HyperLip", 7] = censored
      mean[i, "Ht&HyperLip", 7] = censored / NAR
      
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Ht_0_censor_age <= nx &
            Ht_0_Any == 2 &
            HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
          &
            Cad_0_censor_age > nx &
            nx < Dm_0_censor_age & nx < Death_Censor_Age
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      
      
      
      agelist[[agename]][["Ht&HyperLip"]][["Ht&HyperLip"]] =summary(fit2)$coefficients
      ## from HyperLip and HT to Death directly alone (but could die with CAD)
      
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      
      nar[i, "death", 7] = NAR
      event[i, "death", 7] = censored
      mean[i, "death", 7] = censored / NAR
      
      
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["death"]][["Ht&HyperLip"]] =summary(fit2)$coefficients
      rm(atrisk)
      
    }
    
    
    ## three risk states
    ## from HyperLip and Ht $ Dmto CAD marginal
    
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Ht_0_censor_age &
                        Ht_0_Any == 2 &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2
                      & age > Dm_0_censor_age & Dm_0_Any == 2, ]
    NAR = dim(atrisk)[1]
    
    
    if (nrow(atrisk) < 10)
    {
      agelist[[agename]][["Cad"]][["Ht&HyperLip&Dm"]] = NA
      
      
      agelist[[agename]][["Ht&HyperLip&Dm"]][["Ht&HyperLip&Dm"]] = NA
      
      agelist[[agename]][["death"]][["death"]] = NA
      
      
      
    } else{
      atrisk$yearsinstate = apply(atrisk[, c("HyperLip_0_censor_age",
                                             "Ht_0_censor_age",
                                             "Dm_0_censor_age")], 1, function(x) {
                                               age - max(x[3], max(x[1], x[2]))
                                             })
      atrisk$statin_now = ifelse(atrisk$statin == 1 &
                                   atrisk$statin_age <= nx, 1, 0)
      atrisk$antihtn_now = ifelse(atrisk$antihtn == 1 &
                                    atrisk$htn_age <= nx, 1, 0)
      
      timeatrisk[[agename]][["Ht&HyperLip&Dm"]] = atrisk$yearsinstate
      
      censored = nrow(atrisk[which(Cad_0_censor_age <= nx &
                                     Cad_0_Any == 2 &
                                     nx < Death_Censor_Age), ])
      
      
      nar[i, "Cad", 10] = NAR
      event[i, "Cad", 10] = censored
      mean[i, "Cad", 10] = censored / NAR
      #rm(fit)
      rm(censored)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Cad_0_censor_age <= nx & Cad_0_Any == 2 & nx < Death_Censor_Age
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      
      agelist[[agename]][["Cad"]][["Ht&HyperLip&Dm"]] =summary(fit2)$coefficients
      
      
      # STAY
      censored = dim(atrisk[which(
        Ht_0_censor_age <= nx &
          Ht_0_Any == 2 &
          HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
        &
          Cad_0_censor_age > nx &
          Dm_0_censor_age <= nx &
          Dm_0_Any == 2 & nx < Death_Censor_Age
      ), ])[1]
      
      
      nar[i, "Ht&HyperLip&Dm", 10] = NAR
      event[i, "Ht&HyperLip&Dm", 10] = censored
      mean[i, "Ht&HyperLip&Dm", 10] = censored / NAR
      
      
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(
          Ht_0_censor_age <= nx &
            Ht_0_Any == 2 &
            HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
          &
            Cad_0_censor_age > nx &
            Dm_0_censor_age <= nx &
            Dm_0_Any == 2 & nx < Death_Censor_Age
          ,
          1,
          0
        ) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["Ht&HyperLip&Dm"]][["Ht&HyperLip&Dm"]] =summary(fit2)$coefficients
      
      
      ## from all 3 to death (but could die with CAD)
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      nar[i, "death", 10] = NAR
      event[i, "death", 10] = censored
      mean[i, "death", 10] = censored / NAR
      
      
      rm(censored)
      fit2 = glm(
        family = mode,
        as.formula(paste0("ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0) ~",covariates)),
        data = atrisk
      )
      agelist[[agename]][["death"]][["Ht&HyperLip&Dm"]] =summary(fit2)$coefficients
      
      
    }
   # print(i)
    
  }
  
  
  mylist = list(
    "events" = event,
    "rates" = mean,
    "AR" = nar,
    "model_list" = agelist,
    "yearsinstate" = timeatrisk
    
  )
  return(mylist)
}




plotfunc_prs_sex=function(array,title)
{
  m=melt(array)
  
  # Modify the column names to remove backticks
  colnames(m) <- c("Ages", "PRS Percentile", "Sex", "Predicted per-year Risk")
  
  m$Sex=factor(m$Sex,levels = c("female","male"),labels = c("Female","Male"))
  m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=c(0.2,0.5,0.80),labels = c("Low","Intermediate","High"))
  #m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=levels(as.factor(m$`PRS Percentile`)))
  
  library(ggplot2)
  
  # Create a mapping for line types based on sex
  line_type_mapping <- c("Male" = "solid", "Female" = "dashed")
  
  # Plot the data
  ggplot(m, aes(Ages, `Predicted per-year Risk`, fill = factor(`PRS Percentile`),color = factor(`PRS Percentile`), linetype = factor(Sex))) +
    #geom_point()+
  invisible(stat_smooth(aes(fill=`PRS Percentile`,linetype=Sex,color=`PRS Percentile`))) +
    
    scale_color_manual(values = c("Low" = "#5CB85CFF", "Intermediate" = "#EEA236FF", "High" ="#D43F3AFF"))+
    scale_fill_manual(values = c("Low" = "#5CB85CFF", "Intermediate" = "#EEA236FF", "High" ="#D43F3AFF")) +
    scale_linetype_manual(values = line_type_mapping) +guides(fill = "none")+
    labs(color = "PRS Percentile", linetype = "Sex")+ylab(paste0("Predicted per-year Risk of Transitiont from ",title))+
    theme_minimal()
  
}


## nosmoother

plotfunc_prs_sex_nosmooth=function(array,title)
{
  m=melt(array)
  
  # Modify the column names to remove backticks
  colnames(m) <- c("Ages", "PRS Percentile", "Sex", "Predicted per-year Risk")
  
  m$Sex=factor(m$Sex,levels = c("female","male"),labels = c("Female","Male"))
  m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=c(0.2,0.5,0.80),labels = c("Low","Intermediate","High"))
  #m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=levels(as.factor(m$`PRS Percentile`)))
  
  library(ggplot2)
  
  # Create a mapping for line types based on sex
  line_type_mapping <- c("Male" = "solid", "Female" = "dashed")
  
  # Plot the data
  ggplot(m, aes(Ages, `Predicted per-year Risk`, fill = factor(`PRS Percentile`),color = factor(`PRS Percentile`), linetype = factor(Sex))) +
    geom_point()+
    #invisible(stat_smooth(aes(fill=`PRS Percentile`,linetype=Sex,color=`PRS Percentile`))) +
    
    scale_color_manual(values = c("Low" = "#5CB85CFF", "Intermediate" = "#EEA236FF", "High" ="#D43F3AFF"))+
    scale_fill_manual(values = c("Low" = "#5CB85CFF", "Intermediate" = "#EEA236FF", "High" ="#D43F3AFF")) +
    scale_linetype_manual(values = line_type_mapping) +guides(fill = "none")+
    labs(color = "PRS Percentile", linetype = "Sex")+ylab(paste0("Predicted per-year Risk of Transitiont from ",title," to CAD"))+
    theme_minimal()
  
}




plotfunc_prs_nosex=function(array,title)
{
  m=melt(array)
  
  # Modify the column names to remove backticks
  colnames(m) <- c("Ages", "PRS Percentile", "Predicted per-year Risk")
  
  #m$Sex=factor(m$Sex,levels = c("female","male"),labels = c("Female","Male"))
  m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=c(0.2,0.5,0.80),labels = c("Low","Intermediate","High"))
  #m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=levels(as.factor(m$`PRS Percentile`)))
  
  library(ggplot2)
  
  # Create a mapping for line types based on sex
  #line_type_mapping <- c("Male" = "solid", "Female" = "dashed")
  
  # Plot the data
  ggplot(m, aes(Ages, `Predicted per-year Risk`, fill = factor(`PRS Percentile`),color = factor(`PRS Percentile`))) +
    invisible(stat_smooth(aes(fill=`PRS Percentile`,color=`PRS Percentile`)))+
    
    scale_color_manual(values = c("Low" = "#5CB85CFF", "Intermediate" = "#EEA236FF", "High" ="#D43F3AFF"))+
    scale_fill_manual(values = c("Low" = "#5CB85CFF", "Intermediate" = "#EEA236FF", "High" ="#D43F3AFF")) +
   guides(fill = "none")+
    labs(color = "PRS Percentile")+ylab(paste0("Predicted per-year Risk of Transitiont from ",title," to CAD"))+
    theme_minimal()
  
}



absrisk=function(age,start,stop,cad.prs,sex,modelfit){
  agename=as.character(age)
  mod=modelfit$model_list[[agename]][[stop]][[start]]
  if (start=="Health") {
    new.frame=as.matrix(data.frame(1,cad.prs,sex))
  } else {
    new.frame=as.matrix(data.frame(1,cad.prs,sex))
  }
  ## x%*%B to get log(p/1-p)
  pover1minp=exp(new.frame%*%mod[,1])
  return(pover1minp/(1+pover1minp))
}

absrisk_nosex=function(age,start,stop,cad.prs,modelfit){
  agename=as.character(age)
  mod=modelfit$model_list[[agename]][[stop]][[start]]
  if (start=="Health") {
    new.frame=as.matrix(data.frame(1,cad.prs))
  } else {
    new.frame=as.matrix(data.frame(1,cad.prs))
  }
  ## x%*%B to get log(p/1-p)
  pover1minp=exp(new.frame%*%mod[,1])
  return(pover1minp/(1+pover1minp))
}

stateriskfunc=function(ages,prs_quants,start,stop,modelfit){
  sexes=c(0,1)
  
  riskmat=array(data = NA,dim = c(length(ages),length(prs_quants),length(sexes)),dimnames = list(ages,pnorm(prs_quants),c("female","male")))
  for(i in 1:length(ages)){
    age=ages[i]
    for(j in 1:length(prs_quants)){
      prs=prs_quants[j]
      for(r in 1:length(sexes)){
        sex=sexes[r]
        risk=absrisk(age = age,start = start,stop = stop,cad.prs = prs,sex = sex,modelfit = modelfit)
        riskmat[i,j,r]=risk
      }
    }
  }
  return(riskmat)}



stateriskfunc_nosex=function(ages,prs_quants,start,stop,modelfit){
 
  
  riskmat=array(data = NA,dim = c(length(ages),length(prs_quants)),dimnames = list(ages,pnorm(prs_quants)))
  for(i in 1:length(ages)){
    age=ages[i]
    for(j in 1:length(prs_quants)){
      prs=prs_quants[j]
      risk=absrisk_nosex(age = age,start = start,stop = stop,cad.prs = prs,modelfit = modelfit)
        riskmat[i,j]=risk
      }
    
  }
  return(riskmat)}



plotfunc_prs_sex_compare=function(train,title,test)
{
  m=melt(train)
  
  # Modify the column names to remove backticks
  colnames(m) <- c("Ages", "PRS Percentile", "Sex", "Predicted per-year Risk")
  
  m$Sex=factor(m$Sex,levels = c("female","male"),labels = c("Female","Male"))
  m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=c(0.2,0.5,0.80),labels = c("Low","Intermediate","High"))
  
  
  m2=melt(test)
  # Modify the column names to remove backticks
  colnames(m2) <- c("Ages", "PRS Percentile", "Sex", "Predicted per-year Risk")
  
  m2$Sex=factor(m2$Sex,levels = c("female","male"),labels = c("Female","Male"))
  m2$`PRS Percentile`=factor(m2$`PRS Percentile`,levels=c(0.2,0.5,0.80),labels = c("Low","Intermediate","High"))
  
  
  library(ggplot2)
  
  # Create a mapping for line types based on sex
  line_type_mapping <- c("Male" = "solid", "Female" = "dashed")
  
  # Plot the data
  ggplot(m, aes(Ages, `Predicted per-year Risk`, fill = factor(`PRS Percentile`),color = factor(`PRS Percentile`), linetype = factor(Sex))) +
    stat_smooth(aes(fill=`PRS Percentile`,linetype=Sex,color=`PRS Percentile`)) +
    geom_point(data = m2,mapping = aes(color=`PRS Percentile`,shape=Sex))+
    scale_color_manual(values = c("Low" = "#5CB85CFF", "Intermediate" = "#EEA236FF", "High" ="#D43F3AFF"))+
    scale_fill_manual(values = c("Low" = "#5CB85CFF", "Intermediate" = "#EEA236FF", "High" ="#D43F3AFF")) +
    scale_linetype_manual(values = line_type_mapping) +guides(fill = "none")+
    labs(color = "PRS Percentile", linetype = "Sex")+ylab(paste0("Predicted per-year Risk of Transitiont from ",title," to CAD"))+
    theme_minimal()
  
}



plotfunc_prs_sex_deciles=function(array,title)
{
  m=melt(array)
  
  # Modify the column names to remove backticks
  colnames(m) <- c("Ages", "PRS Percentile", "Sex", "Predicted per-year Risk")
  
  m$Sex=factor(m$Sex,levels = c("female","male"),labels = c("Female","Male"))
  #m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=c(0.2,0.5,0.80),labels = c("Low","Intermediate","High"))
  m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=levels(as.factor(m$`PRS Percentile`)))
  
  library(ggplot2)
  
  # Create a mapping for line types based on sex
  line_type_mapping <- c("Male" = "solid", "Female" = "dashed")
  
  # Plot the data
  ggplot(m, aes(Ages, `Predicted per-year Risk`, fill = factor(`PRS Percentile`),color = factor(`PRS Percentile`), linetype = factor(Sex))) +
    stat_smooth(aes(fill=`PRS Percentile`,linetype=Sex,color=`PRS Percentile`)) +
   scale_color_aaas()+
    scale_fill_aaas()+ scale_linetype_manual(values = line_type_mapping) +guides(fill = "none")+
    labs(color = "PRS Percentile", linetype = "Sex")+ylab(paste0("Predicted per-year Risk of Transitiont from ",title," to CAD"))+
    theme_minimal()
  
}

## for diabetest fit only prs



plotfunc_prs_sex_dm=function(array,title)
{
  m=melt(array)
  
  # Modify the column names to remove backticks
  colnames(m) <- c("Ages", "PRS Percentile", "Sex", "Predicted per-year Risk")
  
  #m$Sex=factor(m$Sex,levels = c("female","male"),labels = c("Female","Male"))
  #m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=c(0.2,0.5,0.80),labels = c("Low","Intermediate","High"))
  m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=levels(as.factor(m$`PRS Percentile`)))
  
  library(ggplot2)
  
  # Create a mapping for line types based on sex
  line_type_mapping <- c("Male" = "solid", "Female" = "dashed")
  
  # Plot the data
  ggplot(m, aes(Ages, `Predicted per-year Risk`, fill = factor(`PRS Percentile`),color = factor(`PRS Percentile`))) +
    stat_smooth(aes(fill=`PRS Percentile`,color=`PRS Percentile`)) +
    scale_color_aaas()+
    scale_fill_aaas()+guides(fill = "none")+
    labs(color = "PRS Percentile", linetype = "Sex")+ylab(paste0("Predicted per-year Risk of Transitiont from ",title," to CAD"))+
    theme_minimal()
  
}




plotfunc_prs_sex_deciles=function(array,title)
{
  m=melt(array)
  
  # Modify the column names to remove backticks
  colnames(m) <- c("Ages", "PRS Percentile", "Sex", "Predicted per-year Risk")
  
  m$Sex=factor(m$Sex,levels = c("female","male"),labels = c("Female","Male"))
  #m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=c(0.2,0.5,0.80),labels = c("Low","Intermediate","High"))
  m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=levels(as.factor(m$`PRS Percentile`)))
  
  library(ggplot2)
  
  # Create a mapping for line types based on sex
  line_type_mapping <- c("Male" = "solid", "Female" = "dashed")
  
  # Plot the data
  ggplot(m, aes(Ages, `Predicted per-year Risk`, fill = factor(`PRS Percentile`),color = factor(`PRS Percentile`), linetype = factor(Sex))) +
    stat_smooth(aes(fill=`PRS Percentile`,linetype=Sex,color=`PRS Percentile`)) +
    scale_color_aaas()+
    scale_fill_aaas()+ scale_linetype_manual(values = line_type_mapping) +guides(fill = "none")+
    labs(color = "PRS Percentile", linetype = "Sex")+ylab(paste0("Predicted per-year Risk of Transitiont from ",title," to CAD"))+
    theme_minimal()
  
}

## for diabetest fit only prs



plotfunc_prs_sex_dm=function(array,title)
{
  m=melt(array)
  
  # Modify the column names to remove backticks
  colnames(m) <- c("Ages", "PRS Percentile", "Predicted per-year Risk")
  
  #m$Sex=factor(m$Sex,levels = c("female","male"),labels = c("Female","Male"))
  #m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=c(0.2,0.5,0.80),labels = c("Low","Intermediate","High"))
  m$`PRS Percentile`=factor(m$`PRS Percentile`,levels=levels(as.factor(m$`PRS Percentile`)))
  
  library(ggplot2)
  
  # Create a mapping for line types based on sex
  line_type_mapping <- c("Male" = "solid", "Female" = "dashed")
  
  # Plot the data
  ggplot(m, aes(Ages, `Predicted per-year Risk`, fill = factor(`PRS Percentile`),color = factor(`PRS Percentile`))) +
    stat_smooth(aes(fill=`PRS Percentile`,color=`PRS Percentile`)) +
    scale_color_aaas()+
    scale_fill_aaas()+guides(fill = "none")+
    labs(color = "PRS Percentile", linetype = "Sex")+ylab(paste0("Predicted per-year Risk of Transitiont from ",title," to CAD"))+
    theme_minimal()
  
}

# prs_quants=qnorm(c(seq(0.1,0.9,by=0.1)))
# prs_quants=qnorm(c(0.25,0.5,0.75))
# 
# s=stateriskfunc_nosex(ages = ages,start = "Dm",stop = "Cad",prs_quants = prs_quants,modelfit = aprs)
# plotfunc_prs_sex_dm(s,title="DM")












absrisk_smoking=function(age,start,stop,cad.prs,smoke,sex,modelfit){
  agename=as.character(age)
  mod=modelfit$model_list[[agename]][[stop]][[start]]
  if (start=="Health") {
    new.frame=as.matrix(data.frame(1,cad.prs,sex,smoke))
  } else {
    new.frame=as.matrix(data.frame(1,cad.prs,sex,smoke))
  }
  ## x%*%B to get log(p/1-p)
  pover1minp=exp(new.frame%*%mod[,1])
  return(pover1minp/(1+pover1minp))
}

stateriskfunc_smoking=function(ages,prs_quants,start,stop,modelfit){
  sexes=c(0,1)
  smoking=c(0,1)
  riskmat=array(data = NA,dim = c(length(ages),length(prs_quants),length(sexes),length(smoking)),dimnames = list(ages,pnorm(prs_quants),c("female","male"),c("none","smoke")))
  for(i in 1:length(ages)){
    age=ages[i]
    for(j in 1:length(prs_quants)){
      prs=prs_quants[j]
      for(r in 1:length(sexes)){
        sex=sexes[r]
        for(s in 1:length(smoking)){
          smoker=smoking[s]
        risk=absrisk_smoking(age = age,start = start,stop = stop,cad.prs = prs,sex = sex,modelfit = modelfit,smoke = smoker)
        riskmat[i,j,r,s]=risk
      }
    }
    }
  }
  return(riskmat)}





## bpmed at baseline




absrisk_smoking_bp=function(age,start,stop,cad.prs,smoke,sex,anti_htn,modelfit){
  agename=as.character(age)
  mod=modelfit$model_list[[agename]][[stop]][[start]]
  
  if (nrow(mod)<5) {
    new.frame=as.matrix(data.frame(1,cad.prs,sex,smoke))
  } else {
    new.frame=as.matrix(data.frame(1,cad.prs,sex,smoke,anti_htn))
  }
  ## x%*%B to get log(p/1-p)
  pover1minp=exp(new.frame%*%mod[,1])
  return(pover1minp/(1+pover1minp))
}

stateriskfunc_smoking_bp=function(ages,prs_quants,start,stop,modelfit){
  sexes=c(0,1)
  smoking=c(0,1)
  bpmeds=c(0,1)
  riskmat=array(data = NA,dim = c(length(ages),length(prs_quants),length(sexes),length(smoking),length(bpmeds)),dimnames = list(ages,pnorm(prs_quants),c("female","male"),c("none","smoke"),c("none","medstart")))
  for(i in 1:length(ages)){
    age=ages[i]
    for(j in 1:length(prs_quants)){
      prs=prs_quants[j]
      for(r in 1:length(sexes)){
        sex=sexes[r]
        for(s in 1:length(smoking)){
          smoker=smoking[s]
          for(h in 1:length(bpmeds)){
            ah=bpmeds[h]
          risk=absrisk_smoking_bp(age = age,start = start,stop = stop,cad.prs = prs,sex = sex,modelfit = modelfit,smoke = smoker,anti_htn=ah)
          riskmat[i,j,r,s,h]=risk}
        }
      }
    }
  }
  return(riskmat)}

# a=plotfunc_prs_sex(s2[,,,1,1],title = "Ht and Non-Smoker and No Med")
# a
# 
# b=plotfunc_prs_sex(s2[,,,1,2],title = "Ht and Non-Smoker and Start Med")
# b
# 
# 
# b=plotfunc_prs_sex(s2[,,,2,2],title = "Ht and Smoker and Start Med")



multipleprsfunc=function(s,prsprobs){
  m=melt(s)
  m$sex=factor(m$Var3,levels = c("female","male"),labels = c("Female","Male"))
  m$prs=factor(m$Var2,levels=prsprobs,labels=prsprobs)
  line_type_mapping <- c("Male" = "solid", "Female" = "dashed")
  g=ggplot(m,aes(Var1,value,color=prs,linetype=sex))+geom_smooth(aes(colour=prs,linetype=sex,fill=prs))+
  scale_linetype_manual(values = line_type_mapping) +
    guides(fill = "none")
  g
}



return_coefficients_states=function(fixed,start,stop){
  
  matrices_list <- list()
  
  # Iterate over each element in the middle dimension of fixed$model_list
  for (i in seq_along(fixed$model_list)) {
    model <- fixed$model_list[[i]]
    
    # Check if the model has the desired structure ($Cad$Health)
    if (is.list(model) && stop %in% names(model) && is.list(model[[stop]]) && start %in% names(model[[stop]])) {
      matrix <- model[[stop]][[start]]
      
      # Check if the extracted element is a matrix
      if (is.matrix(matrix)) {
        matrices_list[[i]] <- matrix
      }
    }
  }
  
  # Remove any NULL elements from the list
  matrices_list <- matrices_list[!sapply(matrices_list, is.null)]
  return(matrices_list)}




compute_empiricalrisk=function(age,age2,cat,df_frame){
  
  atrisk = df_frame[age < Cad_0_censor_age &
                      age < Ht_0_censor_age &
                      age < HyperLip_0_censor_age &
                      age < Dm_0_censor_age, ]
  
  atriskhigh=atrisk[atrisk$int==cat&atrisk$smoke==0,]
  
  rate=mean(atriskhigh$Cad_0_Any==2&atriskhigh$Cad_0_censor_age>age&atriskhigh$Cad_0_censor_age<age2)
  return(rate)
  
}


## pce func


## grab predicited ascvd 

compute_pce_predictedrisk=function(age,df_frame,cat){
  
  df_frame$phenos.enrollment=as.numeric(df_frame$phenos.enrollment)
  atrisk = df_frame[age < Cad_0_censor_age &
                      age < Ht_0_censor_age &
                      age < HyperLip_0_censor_age &
                      age < Dm_0_censor_age  ,]
  
  lb=age-1
  ub=age+1
  atrisk=atrisk[phenos.enrollment>lb&phenos.enrollment<ub,]
  
  atriskint=atrisk[atrisk$int==cat,]
  #rate=mean(na.omit(atriskint$as2))
  rate=mean(na.omit(atriskint$ascvd_10y_accaha))
  return(rate)
}

## matrix and projection functions

matriskfun=function(smoothedplot,ages,quantiles){
  
  ## extract fits
  ggp_data <- ggplot_build(smoothedplot)$data[[1]] 
  cats=levels(as.factor(ggp_data$group))
  print(all.equal(length(cats),length(quantiles)*2))
  yearlyriskmat=matrix(NA,nrow=length(ages),ncol=length(cats))
  yearlyriskredmat=matrix(NA,nrow=length(ages),ncol=length(cats))
  conmat=matrix(NA,nrow=length(ages),ncol=length(cats))
  conredmat=matrix(NA,nrow=length(ages),ncol=length(cats))
  gdt=data.table(ggp_data)
  for(i in 1:length(ages)){
    age=ages[i]
    for(c in 1:length(cats))
    { 
      category=cats[c]
      fit=gdt[round(x,0)==age&group==category,]
      
      yearlyriskmat[i,c]=mean(fit$y)
      yearlyriskredmat[i,c]=0.8*mean(fit$y)
      conmat[i,c]=1-mean(fit$y)
      conredmat[i,c]=1-0.8*mean(fit$y)
    }
  }
  rownames(yearlyriskmat)=rownames(conmat)=rownames(yearlyriskredmat)=rownames(conredmat)=ages
  return(list("yearlyrisk"=yearlyriskmat,"yearlynotrisk"=conmat,"yearlyreducedmat"=yearlyriskredmat,"yearlyreducednotmad"=conredmat))
}

projection_withmat=function(survivalmat,agestart,agestop){
  start.ind=which(rownames(survivalmat)==agestart)
  stop.ind=which(rownames(survivalmat)==agestop)
  rf=data.frame(survivalmat)[start.ind:stop.ind,]
  apply(rf,2,function(x){1-prod(x)})}



compute_projection_level=function(agestart,agestop,smoothedplot,category,ages){
  m=matriskfun(smoothedplot = smoothedplot,quantiles = prs_quants,ages = ages)
  p=projection_withmat(m$yearlynotrisk,agestart = agestart,agestop = agestop)
  return(p[category])
}


compute_pce_predictedrisk_newint=
  function(age,df_frame,cat){
    
    df_frame$phenos.enrollment=as.numeric(df_frame$phenos.enrollment)
    atrisk = df_frame[age < Cad_0_censor_age &
                        age < Ht_0_censor_age &
                        age < HyperLip_0_censor_age &
                        age < Dm_0_censor_age  ,]
    
    lb=age-1
    ub=age+1
    atrisk=atrisk[phenos.enrollment>lb&phenos.enrollment<ub,]
    
    atriskint=atrisk[atrisk$newint==cat,]
    #rate=mean(na.omit(atriskint$as2))
    rate=mean(na.omit(atrisk$ascvd_10y_accaha))
    return(rate)
  }

projection_with_plot=function(plot,ages,quantiles,agestart,agestop){
  m=matriskfun(smoothedplot = plot,ages=ages,quantiles = quantiles)
  p=projection_withmat(m$yearlynotrisk,agestart = agestart,agestop = agestop)
  p=data.frame(p)
  g=as.factor(c("Female","Male"))
  pround=as.factor(round(pnorm(quantiles),1))
  rownames(p)=levels(interaction(g,pround))
  return(p)
}


xyear_projection=function(smoothedplot,category,agestart,agestop)
{
  suppressMessages(expr, classes = "message")
  riskmat=NULL
  for(i in 1:length(seq(agestart:agestop))){
    #print(i)
    age=seq(agestart,agestop)[i]
    #print(age)
    riskmat[i]=extract_smoothed_risk(smoothedplot,category = category,age = age)
  }
  rf=data.frame(riskmat)
  
  rf$survival=1-rf$riskmat
  #suppressMessages()
  return(1-prod(rf$survival))
}

extract_smoothed_risk=function(plot,category,age){
  ggp_data <- ggplot_build(plot)$data[[1]] 
  gdt=data.table(ggp_data)
  fit=gdt[round(x,0)==age&group==category,]
  return(mean(fit$y))}
