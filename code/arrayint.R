array_int = function(df_frame, ages, nstates, mode) {
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
  coef_cont = array(
    data = NA,
    dim = c(length(ages), length(nstates), length(nstates)),
    dimnames = list(ages, nstates, nstates)
  )
  se_cont = array(
    data = NA,
    dim = c(length(ages), length(nstates), length(nstates)),
    dimnames = list(ages, nstates, nstates)
  )
  
  df_frame$cad.prs=scale(df_frame$cad.prs)
  
  
  # from Health (1) to Health directly alone
  
  for (i in 1:length(ages)) {
    age = ages[i]
    nx = age + 1
    atrisk = df_frame[age < Cad_0_censor_age &
                        age < Ht_0_censor_age &
                        age < HyperLip_0_censor_age &
                        age < Dm_0_censor_age, ]
    NAR = dim(atrisk)[1]
    
    if (nrow(atrisk) < 1) {
      coef_cont[i, "Cad", 1] = NA
      se_cont[i, "Cad", 1] = NA
      
      
      coef_cont[i, "Health", 1] = NA
      se_cont[i, "Health", 1] = NA
      
      
      coef_cont[i, "death", 1] = NA
      se_cont[i, "death", 1] = NA
      
    } else{
      censored = dim(atrisk[which(
        Cad_0_censor_age > nx &
          nx < Ht_0_censor_age &
          nx < HyperLip_0_censor_age &
          nx < Dm_0_censor_age & nx < Death_Censor_Age
      ), ])[1]
      
      nar[i, "Health", 1] = NAR
      event[i, "Health", 1] = censored
      mean[i, "Health", 1] = censored / NAR
      fit2 = glm(
        family = mode,
        ifelse(
          Cad_0_censor_age > nx &
            nx < Ht_0_censor_age &
            nx < HyperLip_0_censor_age &
            nx < Dm_0_censor_age & nx < Death_Censor_Age,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "Health", 1] = as.numeric(coef(fit2))[2]
      se_cont[i, "Health", 1] = sqrt(diag(vcov(fit2)))[1]
      
      rm(censored)
      ## from Health (1) to CAD directly alone
      
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    nx < Death_Censor_Age), ])[1]
      nar[i, "Cad", 1] = NAR
      event[i, "Cad", 1] = censored
      mean[i, "Cad", 1] = censored / NAR
      
      rm(censored)
      fit2 = glm(
        family = mode,
        ifelse(
          Cad_0_censor_age <= nx &
            Cad_0_Any == 2 & nx < Death_Censor_Age,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "Cad", 1] = as.numeric(coef(fit2))[2]
      se_cont[i, "Cad", 1] = sqrt(diag(vcov(fit2)))[1]
      
      ## from Health (1) to Death directly alone (but could die with CAD, make CAD and death marginal, i.e., don't make conditions that other conditions don't develop that year)
      
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      nar[i, "death", 1] = NAR
      event[i, "death", 1] = censored
      mean[i, "death", 1] = censored / NAR
      rm(censored)
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", 1] = coef(fit2)[1]
      se_cont[i, "death", 1] = sqrt(diag(vcov(fit2)))[1]
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
    if (nrow(atrisk) < 1) {
      coef_cont[i, "Cad", 2] = NA
      #coef_high[i,"Cad",2]=NA
      se_cont[i, "Cad", 2] = NA
      #se_high[i,"Cad",2]=NA
      
      coef_cont[i, "Ht", 2] = NA
      #coef_high[i,"Ht",2]=NA
      se_cont[i, "Ht", 2] = NA
      #se_high[i,"Ht",2]=NA
      
      coef_cont[i, "death", 2] = NA
      #coef_high[i,"death",2]=NA
      se_cont[i, "death", 2] = NA
      #se_high[i,"death",2]=NA
    } else{
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    nx < Death_Censor_Age), ])[1]
      nar[i, "Cad", 2] = NAR
      event[i, "Cad", 2] = censored
      mean[i, "Cad", 2] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        ifelse(
          Cad_0_censor_age <= nx &
            Cad_0_Any == 2 & nx < Death_Censor_Age,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "Cad", 2] = coef(fit2)[1]
      se_cont[i, "Cad", 2] = sqrt(diag(vcov(fit2)))[1]
      
      
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
        ifelse(
          Ht_0_censor_age <= nx & Ht_0_Any == 2
          &
            Cad_0_censor_age > nx &
            nx < HyperLip_0_censor_age &
            nx < Dm_0_censor_age & nx < Death_Censor_Age,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "Ht", 2] = coef(fit2)[1]
      se_cont[i, "Ht", 2] = sqrt(diag(vcov(fit2)))[1]
      
      
      ## from HT (2) to Death directly alone (but could die with CAD)
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any  == 2 &
                                    Ht_0_censor_age < nx &
                                    Ht_0_Any == 2), ])[1]
      nar[i, "death", 2] = NAR
      event[i, "death", 2] = censored
      mean[i, "death", 2] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", 2] = coef(fit2)[1]
      se_cont[i, "death", 2] = sqrt(diag(vcov(fit2)))[1]
      
      rm(atrisk)
      
    }
    
    
    
    ##  from Hyperlip
    
    
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > HyperLip_0_censor_age &
                        HyperLip_0_Any == 2 &
                        age < Ht_0_censor_age &
                        age < Dm_0_censor_age, ]
    NAR = dim(atrisk)[1]
    
    
    if (nrow(atrisk) < 1) {
      coef_cont[i, "Cad", 3] = NA
      #coef_high[i,"Cad",3]=NA
      se_cont[i, "Cad", 3] = NA
      #se_high[i,"Cad",3]=NA
      
      coef_cont[i, "HyperLip", 3] = NA
      #coef_high[i,"HyperLip",3]=NA
      se_cont[i, "HyperLip", 3] = NA
      #se_high[i,"HyperLip",3]=NA
      
      coef_cont[i, "death", 3] = NA
      #coef_high[i,"death",3]=NA
      se_cont[i, "death", 3] = NA
      #se_high[i,"death",3]=NA
      
    } else{
      ## from Hyperlip (3) to CAD marginal (i.e., could develop other states in the interim year)
      
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2
                                  & nx < Death_Censor_Age), ])[1]
      
      
      nar[i, "Cad", 3] = NAR
      event[i, "Cad", 3] = censored
      mean[i, "Cad", 3] = censored / NAR
      
      rm(censored)
      
      fit2 = glm(
        family = mode,
        ifelse(Cad_0_censor_age <= nx &
                 Cad_0_Any == 2, 1, 0) ~ cad.prs,
        data = atrisk
      )
      
      coef_cont[i, "Cad", 3] = coef(fit2)[1]
      
      se_cont[i, "Cad", 3] = sqrt(diag(vcov(fit2)))[1]
      
      
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
        family = mode,
        ifelse(
          HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
          &
            Cad_0_censor_age > nx &
            nx < Ht_0_censor_age &
            nx < Dm_0_censor_age &
            nx < Death_Censor_Age,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "HyperLip", 3] = coef(fit2)[1]
      se_cont[i, "HyperLip", 3] = sqrt(diag(vcov(fit2)))[1]
      
      ## from HyperLip (2) to Death directly alone (but could die with CAD)
      
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      nar[i, "death", 3] = NAR
      event[i, "death", 3] = censored
      mean[i, "death", 3] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~ cad.prs,
        data = atrisk
      )
      
      coef_cont[i, "death", 3] = coef(fit2)[1]
      se_cont[i, "death", 3] = sqrt(diag(vcov(fit2)))[1]
      rm(atrisk)
      
    }
    
    
    
    ## from Dm2
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age &
                        Dm_0_Any == 2 &
                        age < Ht_0_censor_age &
                        age < HyperLip_0_censor_age, ]
    NAR = dim(atrisk)[1]
    
    
    
    if (nrow(atrisk) < 1) {
      coef_cont[i, "Cad", 4] = NA
      #coef_high[i,"Cad",4]=NA
      se_cont[i, "Cad", 4] = NA
      #se_high[i,"Cad",4]=NA
      
      coef_cont[i, "Dm", 4] = NA
      #coef_high[i,"Dm",4]=NA
      se_cont[i, "Dm", 4] = NA
      #se_high[i,"Dm",4]=NA
      
      coef_cont[i, "death", 4] = NA
      #coef_high[i,"death",4]=NA
      se_cont[i, "death", 4] = NA
      #se_high[i,"death",4]=NA
      
      
    } else{
      # from Dm (4) to CAD directly alone
      
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
        ifelse(
          Cad_0_censor_age <= nx & Cad_0_Any == 2 &
            Dm_0_censor_age <= nx & Dm_0_Any == 2
          ,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      
      coef_cont[i, "Cad", 4] = coef(fit2)[1]
      se_cont[i, "Cad", 4] = sqrt(diag(vcov(fit2)))[1]
      
      
      #rm(fit)
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
      #coef_mid[i,"Dm",4]=coef(fit)[2]
      #coef_high[i,"Dm",4]=coef(fit)[3]
      #se_mid[i,"Dm",4]=sqrt(diag(vcov(fit)))[2]
      #se_high[i,"Dm",4]=sqrt(diag(vcov(fit)))[3]
      
      #rm(fit)
      rm(censored)
      
      fit2 = glm(
        family = mode,
        ifelse(
          Dm_0_censor_age <= nx &
            Dm_0_Any == 2 &
            Cad_0_censor_age > nx &
            nx < HyperLip_0_censor_age & nx < Ht_0_censor_age &
            nx < Death_Censor_Age
          ,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      
      coef_cont[i, "Dm", 4] = coef(fit2)[1]
      se_cont[i, "Dm", 4] = sqrt(diag(vcov(fit2)))[1]
      
      ## from Dm to Death directly alone (but could die with CAD)
      
      censored = dim(atrisk[which(
        Death_Censor_Age <= nx &
          Death_censor_Any == 2#&Dm_0_censor_age<=nx&Dm_0_Any==2
        &
          nx < Ht_0_censor_age & nx < HyperLip_0_censor_age
      ), ])[1]
      
      nar[i, "death", 4] = NAR
      event[i, "death", 4] = censored
      mean[i, "death", 4] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", 4] = coef(fit2)[1]
      se_cont[i, "death", 4] = sqrt(diag(vcov(fit2)))[1]
      
      rm(atrisk)
    }
    
    ## from CAD to death (to do, should we add RF plus CAD as a starting RF? (no, because that would add the complexity of having to add RF + CAD as an ending state and right now we consider CAD as an 'absorbing' state i.e., if you are diagnosed with an additional RF in the interin  year between diagnosis, not considered, ending with CAD is the same as ending with CAD + ...))
    ## from CAD marginal to Death directly
    
    atrisk = df_frame[age > Cad_0_censor_age &
                        Cad_0_Any == 2 & Death_Censor_Age > age, ]
    NAR = dim(atrisk)[1]
    
    
    if (nrow(atrisk) < 1) {
      coef_cont[i, "Cad", 5] = NA
      #coef_high[i,"Cad",5]=NA
      se_cont[i, "Cad", 5] = NA
      #se_high[i,"Cad",5]=NA
      
      
      coef_cont[i, "death", 5] = NA
      #coef_high[i,"death",5]=NA
      se_cont[i, "death", 5] = NA
      #se_high[i,"death",5]=NA
    } else{
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      #fit=glm(family=mode,ifelse( Death_Censor_Age<=nx&Death_censor_Any==2
      # ,1,0)~cad.prs.lev,data = atrisk)
      
      nar[i, "death", 5] = NAR
      event[i, "death", 5] = censored
      mean[i, "death", 5] = censored / NAR
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", 5] = coef(fit2)[1]
      se_cont[i, "death", 5] = sqrt(diag(vcov(fit2)))[1]
      
      censored = dim(atrisk[which(Death_Censor_Age > nx &
                                    Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2), ])[1]
      
      nar[i, "Cad", 5] = NAR
      event[i, "Cad", 5] = censored
      mean[i, "Cad", 5] = censored / NAR
      
      rm(censored)
      fit2 = glm(
        family = mode,
        ifelse(
          Death_Censor_Age > nx &
            Cad_0_censor_age <= nx & Cad_0_Any == 2,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "Cad", 5] = coef(fit2)[1]
      se_cont[i, "Cad", 5] = sqrt(diag(vcov(fit2)))[1]
      rm(atrisk)
    }
    
    ### 2 risk states
    ## from Ht and DM to CAD marginal
    
    atrisk = df_frame[age < Cad_0_censor_age &
                        age > Dm_0_censor_age & Dm_0_Any == 2 &
                        age > Ht_0_censor_age & Ht_0_Any == 2
                      & age < HyperLip_0_censor_age, ]
    NAR = dim(atrisk)[1]
    
    
    if (nrow(atrisk) < 1)
    {
      coef_cont[i, "Cad", 9] = NA
      #coef_high[i,"Cad",9]=NA
      se_cont[i, "Cad", 9] = NA
      #se_high[i,"Cad",9]=NA
      
      coef_cont[i, "Ht&Dm", 9] = NA
      #coef_high[i,"Ht&Dm",9]=NA
      se_cont[i, "Ht&Dm", 9] = NA
      #se_high[i,"Ht&Dm",9]=NA
      
      coef_cont[i, "death", 9] = NA
      #coef_high[i,"death",9]=NA
      se_cont[i, "death", 9] = NA
      #se_high[i,"death",9]=NA
      
      
    } else{
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    nx < Death_Censor_Age), ])[1]
      nar[i, "Cad", 9] = NAR
      event[i, "Cad", 9] = censored
      mean[i, "Cad", 9] = censored / NAR
      rm(censored)
      fit2 = glm(
        family = mode,
        ifelse(
          Death_Censor_Age > nx & Cad_0_censor_age <= nx & Cad_0_Any == 2
          ,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "Cad", 9] = coef(fit2)[1]
      se_cont[i, "Cad", 9] = sqrt(diag(vcov(fit2)))[1]
      
      ## from Ht and Dm to Death directly alone (but could die with CAD)
      
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      #fit=glm(family=mode,ifelse( Death_Censor_Age<=nx&Death_censor_Any==2
      #                ,1,0)~cad.prs.lev,data = atrisk)
      nar[i, "death", 9] = NAR
      event[i, "death", 9] = censored
      mean[i, "death", 9] = censored / NAR
      
      #coef_mid[i,"death",9]=coef(fit)[2]
      #coef_high[i,"death",9]=coef(fit)[3]
      #se_mid[i,"death",9]=sqrt(diag(vcov(fit)))[2]
      #se_high[i,"death",9]=sqrt(diag(vcov(fit)))[3]
      
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", 9] = coef(fit2)[1]
      se_cont[i, "death", 9] = sqrt(diag(vcov(fit2)))[1]
      ## from Ht and Dm to Ht and Dm directly
      
      censored = dim(atrisk[which(
        nx < Cad_0_censor_age &
          nx >= Dm_0_censor_age & Dm_0_Any == 2 &
          nx >= Ht_0_censor_age & Ht_0_Any == 2
        & nx < HyperLip_0_censor_age
      ), ])[1]
      #fit=glm(family=mode,ifelse(nx<Cad_0_censor_age&nx>=Dm_0_censor_age&Dm_0_Any==2&nx>=Ht_0_censor_age&Ht_0_Any==2
      #           &nx<HyperLip_0_censor_age
      #               ,1,0)~cad.prs.lev,data = atrisk)
      nar[i, "Ht&Dm", 9] = NAR
      event[i, "Ht&Dm", 9] = censored
      mean[i, "Ht&Dm", 9] = censored / NAR
      #coef_mid[i,"Ht&Dm",9]=coef(fit)[2]
      #coef_high[i,"Ht&Dm",9]=coef(fit)[3]
      #se_mid[i,"Ht&Dm",9]=sqrt(diag(vcov(fit)))[2]
      #se_high[i,"Ht&Dm",9]=sqrt(diag(vcov(fit)))[3]
      
      
      #rm(fit)
      rm(censored)
      
      fit2 = glm(
        family = mode,
        ifelse(
          nx < Cad_0_censor_age &
            nx >= Dm_0_censor_age & Dm_0_Any == 2 &
            nx >= Ht_0_censor_age & Ht_0_Any == 2
          & nx < HyperLip_0_censor_age
          ,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      
      coef_cont[i, "Ht&Dm", 9] = coef(fit2)[1]
      se_cont[i, "Ht&Dm", 9] = sqrt(diag(vcov(fit2)))[1]
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
    
    if (nrow(atrisk) < 1)
    {
      coef_cont[i, "Cad", 8] = NA
      #coef_high[i,"Cad",8]=NA
      se_cont[i, "Cad", 8] = NA
      #se_high[i,"Cad",8]=NA
      
      coef_cont[i, "HyperLip&Dm", 8] = NA
      #coef_high[i,"HyperLip&Dm" ,8]=NA
      se_cont[i, "HyperLip&Dm" , 8] = NA
      #se_high[i,"HyperLip&Dm" ,8]=NA
      
      coef_cont[i, "death", 8] = NA
      #coef_high[i,"death",8]=NA
      se_cont[i, "death", 8] = NA
      #se_high[i,"death",8]=NA
      
      
    } else{
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      nar[i, "death", 8] = NAR
      event[i, "death", 8] = censored
      mean[i, "death", 8] = censored / NAR
      
      rm(censored)
      
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", 8] = coef(fit2)[1]
      se_cont[i, "death", 8] = sqrt(diag(vcov(fit2)))[1]
      
      
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
        ifelse(
          Cad_0_censor_age <= nx & Cad_0_Any == 2 & Death_Censor_Age > nx
          ,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "Cad", 8] = coef(fit2)[1]
      se_cont[i, "Cad", 8] = sqrt(diag(vcov(fit2)))[1]
      
      
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
        ifelse(
          nx < Cad_0_censor_age &
            nx > Dm_0_censor_age &
            Dm_0_Any == 2 &
            nx > HyperLip_0_censor_age &
            HyperLip_0_Any == 2 & nx < Ht_0_censor_age
          ,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "HyperLip&Dm", 8] = coef(fit2)[1]
      se_cont[i, "HyperLip&Dm", 8] = sqrt(diag(vcov(fit2)))[1]
      
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
    
    
    if (nrow(atrisk) < 1)
    {
      coef_cont[i, "Cad", 7] = NA
      #coef_high[i,"Cad",7]=NA
      se_cont[i, "Cad", 7] = NA
      #se_high[i,"Cad",7]=NA
      
      coef_cont[i, "Ht&HyperLip", 7] = NA
      #coef_high[i,"Ht&HyperLip",7]=NA
      se_cont[i, "Ht&HyperLip", 7] = NA
      #se_high[i,"Ht&HyperLip",7]=NA
      
      coef_cont[i, "death", 7] = NA
      #coef_high[i,"death",7]=NA
      se_cont[i, "death", 7] = NA
      #se_high[i,"death",7]=NA
      
      
    } else{
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
        ifelse(
          Cad_0_censor_age <= nx &
            Cad_0_Any == 2 &
            Ht_0_censor_age <= nx &
            Ht_0_Any == 2 &
            HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2 &
            nx < Death_Censor_Age
          ,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      
      
      coef_cont[i, "Cad", 7] = coef(fit2)[1]
      se_cont[i, "Cad", 7] = sqrt(diag(vcov(fit2)))[1]
      
      
      ## from  HyperLip and Ht to HyperLip and Ht
      
      censored = dim(atrisk[which(
        Ht_0_censor_age <= nx &
          Ht_0_Any == 2 &
          HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
        &
          Cad_0_censor_age > nx &
          nx < Dm_0_censor_age & nx < Death_Censor_Age
      ), ])[1]
      
      #fit=glm(family=mode,ifelse(Ht_0_censor_age<=nx&Ht_0_Any==2&HyperLip_0_censor_age<=nx&HyperLip_0_Any==2
      #       &Cad_0_censor_age>nx&nx<Dm_0_censor_age&nx< Death_Censor_Age
      #       ,1,0)~cad.prs.lev,data = atrisk)
      
      
      ##coef_mid[i,"Ht&HyperLip",7]=coef(fit)[2]
      #coef_high[i,"Ht&HyperLip",7]=coef(fit)[3]
      #se_mid[i,"Ht&HyperLip",7]=sqrt(diag(vcov(fit)))[2]
      #se_high[i,"Ht&HyperLip",7]=sqrt(diag(vcov(fit)))[3]
      
      nar[i, "Ht&HyperLip", 7] = NAR
      event[i, "Ht&HyperLip", 7] = censored
      mean[i, "Ht&HyperLip", 7] = censored / NAR
      
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        ifelse(
          Ht_0_censor_age <= nx &
            Ht_0_Any == 2 &
            HyperLip_0_censor_age <= nx & HyperLip_0_Any == 2
          &
            Cad_0_censor_age > nx &
            nx < Dm_0_censor_age & nx < Death_Censor_Age
          ,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      
      
      coef_cont[i, "Ht&HyperLip", 7] = coef(fit2)[1]
      se_cont[i, "Ht&HyperLip", 7] = sqrt(diag(vcov(fit2)))[1]
      ## from HyperLip and HT to Death directly alone (but could die with CAD)
      
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      #fit=glm(family=mode,ifelse( Death_Censor_Age<=nx&Death_censor_Any==2
      ##            ,1,0)~cad.prs.lev,data = atrisk)
      
      
      ##coef_mid[i,"death",7]=coef(fit)[2]
      #coef_high[i,"death",7]=coef(fit)[3]
      #se_mid[i,"death",7]=sqrt(diag(vcov(fit)))[2]
      #se_high[i,"death",7]=sqrt(diag(vcov(fit)))[3]
      
      nar[i, "death", 7] = NAR
      event[i, "death", 7] = censored
      mean[i, "death", 7] = censored / NAR
      
      
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", 7] = coef(fit2)[1]
      se_cont[i, "death", 7] = sqrt(diag(vcov(fit2)))[1]
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
    
    
    if (nrow(atrisk) < 1)
    {
      coef_cont[i, "Cad", 10] = NA
      #coef_high[i,"Cad",10]=NA
      se_cont[i, "Cad", 10] = NA
      #se_high[i,"Cad",10]=NA
      
      coef_cont[i, "Ht&HyperLip&Dm", 10] = NA
      #coef_high[i,"Ht&HyperLip&Dm",10]=NA
      se_cont[i, "Ht&HyperLip&Dm", 10] = NA
      #se_high[i,"Ht&HyperLip&Dm",10]=NA
      
      coef_cont[i, "death", 10] = NA
      #coef_high[i,"death",10]=NA
      se_cont[i, "death", 10] = NA
      #se_high[i,"death",10]=NA
      
      
    } else{
      # CAD marg
      #fit=glm(family=mode,ifelse(Cad_0_censor_age<=nx&Cad_0_Any==2&nx< Death_Censor_Age
      #                 ,1,0)~cad.prs.lev,data = atrisk)
      censored = nrow(atrisk[which(Cad_0_censor_age <= nx &
                                     Cad_0_Any == 2 &
                                     nx < Death_Censor_Age), ])
      
      #coef_mid[i,"Cad",10]=coef(fit)[2]
      #coef_high[i,"Cad",10]=coef(fit)[3]
      #se_mid[i,"Cad",10]=sqrt(diag(vcov(fit)))[2]
      #se_high[i,"Cad",10]=sqrt(diag(vcov(fit)))[3]
      
      
      nar[i, "Cad", 10] = NAR
      event[i, "Cad", 10] = censored
      mean[i, "Cad", 10] = censored / NAR
      #rm(fit)
      rm(censored)
      
      fit2 = glm(
        family = mode,
        ifelse(
          Cad_0_censor_age <= nx & Cad_0_Any == 2 & nx < Death_Censor_Age
          ,
          1,
          0
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "Cad", 10] = coef(fit2)[1]
      se_cont[i, "Cad", 10] = sqrt(diag(vcov(fit2)))[1]
      
      
      
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
        ifelse(
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
        ) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "Ht&HyperLip&Dm", 10] = coef(fit2)[1]
      se_cont[i, "Ht&HyperLip&Dm", 10] = sqrt(diag(vcov(fit2)))[1]
      
      
      ## from all 3 to death (but could die with CAD)
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      nar[i, "death", 10] = NAR
      event[i, "death", 10] = censored
      mean[i, "death", 10] = censored / NAR
      
      
      rm(censored)
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", 10] = coef(fit2)[1]
      se_cont[i, "death", 10] = sqrt(diag(vcov(fit2)))[1]
      
    }
    #print(i)
    
  }
  
  
  mylist = list(
    "events" = event,
    "rates" = mean,
    "AR" = nar,
    "coef_cont" = coef_cont,
    "se_cont" = se_cont
  )
  return(mylist)
}
