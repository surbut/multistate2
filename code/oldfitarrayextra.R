
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




compute_empiricalrisk=function(age,age2,atrisk){
  
  rate=mean(atrisk$Cad_0_Any==2&atrisk$Cad_0_censor_age>age&atrisk$Cad_0_censor_age<age2)
  return(rate)
  
}


## pce func


## grab predicited ascvd 

compute_pce_predictedrisk=function(age,atrisk){
  
  lb=age-2
  ub=age+2
  atrisk=atrisk[phenos.enrollment>lb&phenos.enrollment<ub,]
  rate=mean(na.omit(atrisk$ascvd_10y_accaha))
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
