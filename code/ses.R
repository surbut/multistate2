
nstates=c("Health","oneRF","twoRF","threeRF","Cad","death")

arrayrf = function(df_frame, ages, nstates, mode) {
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
      se_cont[i, "Health", 1] = sqrt(diag(vcov(fit2)))[2]
      
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
      se_cont[i, "Cad", 1] = sqrt(diag(vcov(fit2)))[2]
      
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
      coef_cont[i, "death", 1] = coef(fit2)[2]
      se_cont[i, "death", 1] = sqrt(diag(vcov(fit2)))[2]
      rm(atrisk)
    }
    
    
    ## From one risk states to 2,3,4,CAD or death
    
    
   htage=df_frame[which(age < df_frame$Cad_0_censor_age &
                        age > df_frame$Ht_0_censor_age &
                        df_frame$Ht_0_Any == 2 &
                        age < df_frame$HyperLip_0_censor_age &
                        age < df_frame$Dm_0_censor_age),]
    
    hlage=df_frame[which(age < df_frame$Cad_0_censor_age &
                        age > df_frame$HyperLip_0_censor_age &
                        df_frame$HyperLip_0_Any == 2 &
                        age < df_frame$Ht_0_censor_age &
                        age < df_frame$Dm_0_censor_age),]
    
    dmage=df_frame[which(age < df_frame$Cad_0_censor_age &
                        age > df_frame$Dm_0_censor_age &
                        df_frame$Dm_0_Any == 2 &
                        age < df_frame$Ht_0_censor_age &
                        age < df_frame$HyperLip_0_censor_age),]


    u=union(hlage$identifier,union(htage$identifier,dmage$identifier))
    df_frame=data.table(df_frame)
    atrisk=df_frame[identifier%in%u,]
    NAR = dim(atrisk)[1]
    
    if (nrow(atrisk) < 1) {
      coef_cont[i, "Cad", "oneRF"] = NA
      #coef_high[i,"Cad",2]=NA
      se_cont[i, "Cad", "oneRF"] = NA
      #se_high[i,"Cad",2]=NA
      
      coef_cont[i, "oneRF", "oneRF"] = NA
      #coef_high[i,"Ht",2]=NA
      se_cont[i, "oneRF", "oneRF"] = NA
      #se_high[i,"Ht",2]=NA
      
      coef_cont[i, "death", "oneRF"] = NA
      #coef_high[i,"death",2]=NA
      se_cont[i, "death", "oneRF"] = NA
      #se_high[i,"death",2]=NA
    } else{
      
      # transition to CAD
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    nx < Death_Censor_Age), ])[1]
      nar[i, "Cad", "oneRF"] = NAR
      event[i, "Cad", "oneRF"] = censored
      mean[i, "Cad", "oneRF"] = censored / NAR
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
      coef_cont[i, "Cad", "oneRF"] = coef(fit2)[2]
      se_cont[i, "Cad", "oneRF"] = sqrt(diag(vcov(fit2)))[2]
      
      
      ## from single (2) to single directly alone
      
      
      htage_n=atrisk[(nx < atrisk$Cad_0_censor_age &
                              atrisk$Ht_0_censor_age <= nx &
                              atrisk$Ht_0_Any == 2 &
                              nx < atrisk$HyperLip_0_censor_age &
                              nx < atrisk$Dm_0_censor_age),]
      
      hlage_n=atrisk[(nx < atrisk$Cad_0_censor_age &
                              atrisk$HyperLip_0_censor_age <= nx &
                              atrisk$HyperLip_0_Any == 2 &
                              nx < atrisk$Ht_0_censor_age &
                              nx < atrisk$Dm_0_censor_age),]
      
      dmage_n=atrisk[(nx < atrisk$Cad_0_censor_age &
                              atrisk$Dm_0_censor_age <= nx &
                              atrisk$Dm_0_Any == 2 &
                              nx < atrisk$Ht_0_censor_age &
                              nx < atrisk$HyperLip_0_censor_age),]
      
      u=union(hlage_n$identifier,union(htage_n$identifier,dmage_n$identifier))
      censormat=df_frame[identifier%in%u,]
      atrisk$singlerisk=rep(0,nrow(atrisk))
      atrisk[identifier%in%u,"singlerisk"]=1
      #a=apply(atrisk[identifier%in%u,c("Dm_0_censor_age","Ht_0_censor_age","HyperLip_0_censor_age")],1,function(x){sum(x<=nx)})
      censored = dim(censormat)[1]
      
      nar[i, "oneRF", "oneRF"] = NAR
      event[i, "oneRF", "oneRF"] = censored
      mean[i, "oneRF", "oneRF"] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        singlerisk
         ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "oneRF", 2] = coef(fit2)[2]
      se_cont[i, "oneRF", 2] = sqrt(diag(vcov(fit2)))[2]
      
      
      ## from HT (2) to Death directly alone (but could die with CAD)
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any  == 2 &
                                    Ht_0_censor_age < nx &
                                    Ht_0_Any == 2), ])[1]
      nar[i, "death", "oneRF"] = NAR
      event[i, "death", "oneRF"] = censored
      mean[i, "death", "oneRF"] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", "oneRF"] = coef(fit2)[2]
      se_cont[i, "death", "oneRF"] = sqrt(diag(vcov(fit2)))[2]
      
      rm(atrisk)
      
    }
    
    
## from CAD to death (to do, should we add RF plus CAD as a starting RF? (no, because that would add the complexity of having to add RF + CAD as an ending state and right now we consider CAD as an 'absorbing' state i.e., if you are diagnosed with an additional RF in the interin  year between diagnosis, not considered, ending with CAD is the same as ending with CAD + ...))
    ## from CAD marginal to Death directly
    
    atrisk = df_frame[age > Cad_0_censor_age &
                        Cad_0_Any == 2 & Death_Censor_Age > age, ]
    NAR = dim(atrisk)[1]
    
    
    if (nrow(atrisk) < 1) {
      coef_cont[i, "Cad", "Cad"] = NA
      #coef_high[i,"Cad",5]=NA
      se_cont[i, "Cad", "Cad"] = NA
      #se_high[i,"Cad",5]=NA
      
      
      coef_cont[i, "death", "Cad"] = NA
      #coef_high[i,"death",5]=NA
      se_cont[i, "death", "Cad"] = NA
      #se_high[i,"death",5]=NA
    } else{
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      #fit=glm(family=mode,ifelse( Death_Censor_Age<=nx&Death_censor_Any==2
      # ,1,0)~cad.prs.lev,data = atrisk)
      
      nar[i, "death", "Cad"] = NAR
      event[i, "death", "Cad"] = censored
      mean[i, "death", "Cad"] = censored / NAR
      rm(censored)
      #rm(fit)
      
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2
               , 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", "Cad"] = coef(fit2)[2]
      se_cont[i, "death", "Cad"] = sqrt(diag(vcov(fit2)))[2]
      
      censored = dim(atrisk[which(Death_Censor_Age > nx &
                                    Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2), ])[1]
      
      nar[i, "Cad", "Cad"] = NAR
      event[i, "Cad", "Cad"] = censored
      mean[i, "Cad", "Cad"] = censored / NAR
      
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
      coef_cont[i, "Cad", "Cad"] = coef(fit2)[2]
      se_cont[i, "Cad", "Cad"] = sqrt(diag(vcov(fit2)))[2]
      rm(atrisk)
    }
    
    ### 2 risk states
    ## from Ht and DM to CAD marginal
    
    
    atrisk_dmh = df_frame[which(age < Cad_0_censor_age &
                        age > Dm_0_censor_age & Dm_0_Any == 2 &
                        age > Ht_0_censor_age & Ht_0_Any == 2
                      & age < HyperLip_0_censor_age), ]
    
    atrisk_dmhl = df_frame[which(age < Cad_0_censor_age &
                                  age > Dm_0_censor_age & Dm_0_Any == 2 &
                                  age > HyperLip_0_censor_age & HyperLip_0_Any == 2
                                & age < Ht_0_censor_age), ]
    atrisk_hthl = df_frame[which(age < Cad_0_censor_age &
                                   age > Ht_0_censor_age & Ht_0_Any == 2 &
                                   age > HyperLip_0_censor_age & HyperLip_0_Any == 2
                                 & age < Dm_0_censor_age), ]
    
    
    
    
    u=union(atrisk_dmh$identifier,union(atrisk_dmhl$identifier,atrisk_hthl$identifier))
    atrisk=df_frame[identifier%in%u,]
    NAR = dim(atrisk)[1]
    
    
    
    if (nrow(atrisk) < 1)
    {
      coef_cont[i, "Cad", "twoRF",] = NA
      #coef_high[i,"Cad",9]=NA
      se_cont[i, "Cad", "twoRF",] = NA
      #se_high[i,"Cad",9]=NA
      
      coef_cont[i, "twoRF",, "twoRF"] = NA
      #coef_high[i,"Ht&Dm",9]=NA
      se_cont[i, "twoRF", "twoRF"] = NA
      #se_high[i,"Ht&Dm",9]=NA
      
      coef_cont[i, "death", "twoRF"] = NA
      #coef_high[i,"death",9]=NA
      se_cont[i, "death", "twoRF"] = NA
      #se_high[i,"death",9]=NA
      
      
    } else{
      censored = dim(atrisk[which(Cad_0_censor_age <= nx &
                                    Cad_0_Any == 2 &
                                    nx < Death_Censor_Age), ])[1]
      nar[i, "Cad", "twoRF"] = NAR
      event[i, "Cad", "twoRF"] = censored
      mean[i, "Cad", "twoRF"] = censored / NAR
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
      coef_cont[i, "Cad", "twoRF"] = coef(fit2)[2]
      se_cont[i, "Cad", "twoRF"] = sqrt(diag(vcov(fit2)))[2]
      
      ## from Ht and Dm to Death directly alone (but could die with CAD)
      
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      #fit=glm(family=mode,ifelse( Death_Censor_Age<=nx&Death_censor_Any==2
      #                ,1,0)~cad.prs.lev,data = atrisk)
      nar[i, "death", "twoRF"] = NAR
      event[i, "death", "twoRF"] = censored
      mean[i, "death", "twoRF"] = censored / NAR
      
      #coef_mid[i,"death","twoRF"]=coef(fit)[2]
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
      coef_cont[i, "death", "twoRF"] = coef(fit2)[2]
      se_cont[i, "death", "twoRF"] = sqrt(diag(vcov(fit2)))[2]
      
      
      ## from two risk to two risk
      
      ## from single (2) to single directly alone
      
      
      hthlage_n=atrisk[(nx < atrisk$Cad_0_censor_age &
                        atrisk$Ht_0_censor_age <= nx &
                        atrisk$Ht_0_Any == 2 &
                        atrisk$HyperLip_0_censor_age <= nx &
                        atrisk$HyperLip_0_censor_age ==2 &
                        nx < atrisk$Dm_0_censor_age),]
      
      hldmage_n=atrisk[(nx < atrisk$Cad_0_censor_age &
                        atrisk$HyperLip_0_censor_age <= nx &
                        atrisk$HyperLip_0_Any == 2 &
                        atrisk$Dm_0_censor_age <= nx &
                          atrisk$Dm_0_Any==2 &
                        nx < atrisk$Ht_0_censor_age),]
      
      dmhtage_n=atrisk[(nx < atrisk$Cad_0_censor_age &
                        atrisk$Dm_0_censor_age <= nx &
                        atrisk$Dm_0_Any == 2 &
                        atrisk$Ht_0_censor_age <= nx &
                        atrisk$Ht_0_Any == 2 &
                        nx < atrisk$HyperLip_0_censor_age),]
      
      u=union(hthlage_n$identifier,union(hldmage_n$identifier,dmhtage_n$identifier))
      censormat=df_frame[identifier%in%u,]
      atrisk$doublerisk=rep(0,nrow(atrisk))
      atrisk[identifier%in%u,"doublerisk"]=1
      censored = dim(censormat)[1]
      
      nar[i, "twoRF", "twoRF"] = NAR
      event[i, "twoRF","twoRF"] = censored
      mean[i, "twoRF", "twoRF"] = censored / NAR
      rm(censored)
      
      fit2 = glm(
        family = mode,
        doublerisk
        ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "twoRF", "twoRF"] = coef(fit2)[2]
      se_cont[i, "twoRF", "twoRF"] = sqrt(diag(vcov(fit2)))[2]
      rm(atrisk)}
      
      
      
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
      coef_cont[i, "Cad", "threeRF"] = NA
      #coef_high[i,"Cad",10]=NA
      se_cont[i, "Cad", "threeRF"] = NA
      #se_high[i,"Cad",10]=NA
      
      coef_cont[i, "threeRF", "threeRF"] = NA
      #coef_high[i,"Ht&HyperLip&Dm",10]=NA
      se_cont[i, "threeRF", "threeRF"] = NA
      #se_high[i,"Ht&HyperLip&Dm",10]=NA
      
      coef_cont[i, "death", "threeRF"] = NA
      #coef_high[i,"death",10]=NA
      se_cont[i, "death", "threeRF"] = NA
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
      
      
      nar[i, "Cad", "threeRF"] = NAR
      event[i, "Cad", "threeRF"] = censored
      mean[i, "Cad", "threeRF"] = censored / NAR
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
      coef_cont[i, "Cad", "threeRF"] = coef(fit2)[2]
      se_cont[i, "Cad", "threeRF"] = sqrt(diag(vcov(fit2)))[2]
      
      
      
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
      
      
      nar[i, "threeRF","threeRF"] = NAR
      event[i, "threeRF","threeRF"] = censored
      mean[i, "threeRF","threeRF"] = censored / NAR
      
      
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
      coef_cont[i, "threeRF","threeRF"] = coef(fit2)[2]
      se_cont[i, "threeRF","threeRF"] = sqrt(diag(vcov(fit2)))[2]
      
      
      ## from all 3 to death (but could die with CAD)
      censored = dim(atrisk[which(Death_Censor_Age <= nx &
                                    Death_censor_Any == 2), ])[1]
      
      nar[i, "death", "threeRF"] = NAR
      event[i, "death", "threeRF"] = censored
      mean[i, "death", "threeRF"] = censored / NAR
      
      
      rm(censored)
      fit2 = glm(
        family = mode,
        ifelse(Death_Censor_Age <= nx &
                 Death_censor_Any == 2, 1, 0) ~ cad.prs,
        data = atrisk
      )
      coef_cont[i, "death", "threeRF"] = coef(fit2)[2]
      se_cont[i, "death", "threeRF"] = sqrt(diag(vcov(fit2)))[2]
      
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

## how to run
# df_final = data.table(readRDS("~/Dropbox/pheno_dir/output/merged_pheno_censor_final.rds"))
# reg=fread("~/Dropbox/ukb+gp/gp_registrations.txt")
# df_final=df_final[identifier%in%reg$eid,]
# # restrict to regsitered in primary care
# ages = c(40:80)
# 
# nstates=c("Health","oneRF","twoRF","threeRF","Cad","death")
# 
# abinom = arrayrf(df_frame = df_final,
#                     ages = ages,
#                     nstates = nstates,mode = "binomial")
