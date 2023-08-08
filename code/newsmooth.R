
### extract average coefficients from plot

coefplotsmooth=function(ages,start,stop,modelfit){
  require(ggplot2)
agenames=as.character(c(ages))
s=sapply(agenames,function(x){modelfit$model_list[[x]][[stop]][[start]][,"Estimate"]})
s=data.frame(s)
colnames(s)=agenames
s=t(s)
melt=melt(s)
g=ggplot(melt,aes(Var1,value,col=Var2))+stat_smooth()+geom_point()
m=ggplot_build(g)$data[[1]]
return(list("mat"=m,"plot"=g))}

### return P list of smoothed per age coefficient
coefsmooth=
  function(start,stop,ages,modelfit){
  m=coefplotsmooth(ages = ages,start = start,stop = stop,modelfit = modelfit)$m
  ##grab mean for each age from smoothed
  nterms=as.numeric(levels(as.factor(m$group)))
  returnlist=lapply(nterms,function(t){
  df=data.frame(m[m$group==t,]%>%group_by(round(x,0))%>%summarise(mean(y)))
  rownames(df)=as.character(c(ages))
  colnames(df)=c("age","avgcoef")
  return(df)
  })
  agenames=as.character(c(ages))
  s=sapply(agenames,function(x){modelfit$model_list[[x]][[stop]][[start]][,"Estimate"]})
  s=data.frame(s)
names(returnlist)=rownames(s)
return(returnlist)
  }


### return matrix of smoothed coefficeints
#library(purrr)

convertlistmat=function(smoothedlist){
  require(purrr)
list_of_dfs <- lapply(smoothedlist, as.data.frame)
final_df <- reduce(list_of_dfs, function(x, y) merge(x, y, by="age"))
final_matrix <- as.matrix(final_df)
names(final_df)[-1]=names(smoothedlist)
return(final_df)}


## now compute preidction for an indiviaul patient for a given year
## TODO: make generic for any coefficients

# compute_prediction <- function(coefficients_matrix, age, prs_quant, sex, smoking) {
#   # Find the row that corresponds to the given age
#   age_row <- coefficients_matrix[coefficients_matrix[, 'age'] == age, ]
#   
#   # Extract coefficients
#   intercept_coef <- age_row[,"(Intercept)"]
#   genetics_coef <- age_row[, "cad.prs"]
#   sex_coef <- age_row[, "f.31.0.01"]
#   smoking_coef <- age_row[, "smoke"]
#   sexnum=ifelse(sex=="male",1,0)
#   # Compute the prediction using the linear model
#   prediction <- intercept_coef + genetics_coef * prs_quant + sex_coef * sexnum + smoking_coef * smoking
#   p=exp(prediction)/(1+exp(prediction))
#   return(p)
# }

## now compute preidction for an indiviaul patient for a given year nesting the coefficient matrix production

compute_prediction_product <- function(modelfit,start,stop, agesmooth,agepredinterval, covariatevector) {
  
  smoothedlist=coefsmooth(start = start,stop = stop,ages = agesmooth,modelfit = modelfit)
  coefficients_matrix=convertlistmat(smoothedlist = smoothedlist)
  # Subset the matrix to include only the rows within the given age interval
  age_rows <- coefficients_matrix[coefficients_matrix[, 'age'] %in% agepredinterval, ]
  ## make sure the covariate vector has terms in the same order as the model
  xmat=matrix(as.numeric(c(1,covariatevector)))

 # Compute the prediction using the linear model
  logodds <-as.matrix(age_rows[,-1])%*%xmat
  prediction=exp(logodds)/(1+exp(logodds))
  prediction_not=1-prediction
  # Compute the product of the predictions
  prediction_product <- prod(prediction_not)
  
  return(1-prediction_product)
}


### compute prediction using smoothed coefficients

absrisksmoothedcoef=function(age,start,stop,cad.prs,sex,smoke,smoothedlist){
  agename=as.character(age)
  mod=sapply(seq(1:length(smoothedlist)),function(x){
    smoothedlist[[x]][agename,"avgcoef"]
  })
  names(mod)=names(smoothedlist)
  mod=data.frame(mod)
  
  if (start=="Health") {
    new.frame=as.matrix(data.frame(1,cad.prs,sex,smoke))
  } else {
    new.frame=as.matrix(data.frame(1,cad.prs,sex,smoke))
  }
  ## x%*%B to get log(p/1-p)
  pover1minp=exp(new.frame%*%mod[,1])
  return(pover1minp/(1+pover1minp))
}



### produce state array for enumerated PRS quants


stateriskfunc_smoking_smoothedcoef=function(ages,prs_quants,start,stop,modelfit){
  ## extract smoothed coefficients
  smoothedlist=coefsmooth(start = start,stop = stop,ages = ages,modelfit = modelfit)
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
          risk=absrisksmoothedcoef(age = age,start = start,stop = stop,cad.prs = prs,sex = sex,smoke = smoker,smoothedlist = smoothedlist)
          riskmat[i,j,r,s]=risk
        }
      }
    }}
  return(riskmat)
}

# s=stateriskfunc_smoking_smoothedcoef(ages = ages,prs_quants = prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
# p1=multipleprsfunc(s = s[,,,"none"],prsprobs = pnorm(prs_quants))
# smoothedlist=coefsmooth(start = start,stop = stop,ages = ages,modelfit = modelfit)
# 

### smoothing occurs beofre with coefficients so no need to extract plotbuilder
matriskfunc_coef=function(smoothedplot,ages,quantiles){
  
  ## extract fits
  ggp_data <- smoothedplot$data
  ggp_data$group=interaction(as.factor(ggp_data$Var3),as.factor(ggp_data$Var2))
  levels(ggp_data$group)=c(1:(length(quantiles)*2))
  cats=as.numeric(levels(as.factor(ggp_data$group)))
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
      #print(i)
      #print(c)
      category=cats[c]
      fit=gdt[Var1==age&group==category,]
      
      yearlyriskmat[i,c]=fit$value
      yearlyriskredmat[i,c]=0.8*fit$value
      conmat[i,c]=1-fit$value
      conredmat[i,c]=1-0.8*fit$value
    }
  }
  rownames(yearlyriskmat)=rownames(conmat)=rownames(yearlyriskredmat)=rownames(conredmat)=ages
  return(list("yearlyrisk"=yearlyriskmat,"yearlynotrisk"=conmat,"yearlyreducedmat"=yearlyriskredmat,"yearlyreducednotmad"=conredmat))
}


projection_with_plotcoef=function(plot,ages,quantiles,agestart,agestop){
  m=matriskfunc_coef(smoothedplot = plot,ages=ages,quantiles = quantiles)
  p=projection_withmat(m$yearlynotrisk,agestart = agestart,agestop = agestop)
  p=data.frame(p)
  g=as.factor(c("Female","Male"))
  pround=as.factor(round(pnorm(quantiles),1))
  rownames(p)=levels(interaction(g,pround))
  return(p)
}

# ## create states with smoothed coefficients
# 
# s=stateriskfunc_smoking_smoothedcoef(ages,prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
# p1=multipleprsfunc(s[,,,1],prsprobs = pnorm(prs_quants))
# pc=projection_with_plotcoef(p1,ages,quantiles = prs_quants,agestart =40,agestop = 80)
# 
# 
# s2=stateriskfunc_smoking(ages,prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
# p2=multipleprsfunc(s2[,,,1],prsprobs = pnorm(prs_quants))
# po=projection_with_plot(p2,ages = ages,quantiles = prs_quants,agestart = 40,agestop = 80)
# plot(pc$p,po$p,xlim=c(0,1),ylim=c(0,1),xlab="Coef Smoothing",ylab="Level Smoothing",main="Ht")
# abline(0,1)
# 
# 
# s=stateriskfunc_smoking_smoothedcoef(ages,prs_quants,start = "Ht",stop = "Cad",modelfit = fixedsmoke)
# p1=multipleprsfunc(s[,,,1],prsprobs = pnorm(prs_quants))
# pc=projection_with_plotcoef(p1,ages,quantiles = prs_quants,agestart =40,agestop = 80)
# 
# s2=stateriskfunc_smoking(ages,prs_quants,start = "Ht",stop = "Cad",modelfit = fixedsmoke)
# p2=multipleprsfunc(s2[,,,1],prsprobs = pnorm(prs_quants))
# po=projection_with_plot(p2,ages = ages,quantiles = prs_quants,agestart = 40,agestop = 80)
# plot(pc$p,po$p,xlim=c(0,1),ylim=c(0,1),xlab="Coef Smoothing",ylab="Level Smoothing",main="Ht")
# abline(0,1)
# 
# 
# s=stateriskfunc_smoking_smoothedcoef(ages,prs_quants,start = "HyperLip",stop = "Cad",modelfit = fixedsmoke)
# p1=multipleprsfunc(s[,,,1],prsprobs = pnorm(prs_quants))
# pc=projection_with_plotcoef(p1,ages,quantiles = prs_quants,agestart =40,agestop = 80)
# 
# s2=stateriskfunc_smoking(ages,prs_quants,start = "HyperLip",stop = "Cad",modelfit = fixedsmoke)
# p2=multipleprsfunc(s2[,,,1],prsprobs = pnorm(prs_quants))
# po=projection_with_plot(p2,ages = ages,quantiles = prs_quants,agestart = 40,agestop = 80)
# plot(pc$p,po$p,xlim=c(0,1),ylim=c(0,1),xlab="Coef Smoothing",ylab="Level Smoothing",main="HyperLip")
# abline(0,1)
# 
# s=stateriskfunc_smoking_smoothedcoef(ages,prs_quants,start = "Dm",stop = "Cad",modelfit = fixedsmoke)
# p1=multipleprsfunc(s[,,,1],prsprobs = pnorm(prs_quants))
# pc=projection_with_plotcoef(p1,ages,quantiles = prs_quants,agestart =40,agestop = 80)
# 
# s2=stateriskfunc_smoking(ages,prs_quants,start = "Dm",stop = "Cad",modelfit = fixedsmoke)
# p2=multipleprsfunc(s2[,,,1],prsprobs = pnorm(prs_quants))
# po=projection_with_plot(p2,ages = ages,quantiles = prs_quants,agestart = 40,agestop = 80)
# plot(pc$p,po$p,xlim=c(0,1),ylim=c(0,1),xlab="Coef Smoothing",ylab="Level Smoothing",main="Dm")
# abline(0,1)
# 
# 
# s=stateriskfunc_smoking_smoothedcoef(ages,prs_quants,start = "Ht&Dm",stop = "Cad",modelfit = fixedsmoke)
# p1=multipleprsfunc(s[,,,1],prsprobs = pnorm(prs_quants))
# pc=projection_with_plotcoef(p1,ages,quantiles = prs_quants,agestart =40,agestop = 80)
# 
# s2=stateriskfunc_smoking(ages,prs_quants,start = "Ht&Dm",stop = "Cad",modelfit = fixedsmoke)
# p2=multipleprsfunc(s2[,,,1],prsprobs = pnorm(prs_quants))
# po=projection_with_plot(p2,ages = ages,quantiles = prs_quants,agestart = 40,agestop = 80)
# plot(pc$p,po$p,xlim=c(0,1),ylim=c(0,1),xlab="Coef Smoothing",ylab="Level Smoothing",main="Ht&Dm")
# abline(0,1)
# 
# 
# ## 
# 
# ages=40:80
# s=stateriskfunc_smoking(ages,prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
# p=multipleprsfunc(s[,,,1],prsprobs = pnorm(prs_quants))
# m=matriskfun(p,ages,quantiles = prs_quants)
# agesint=seq(40,70,by=1)
# ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
# lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
# 
# for(i in 1:length(agesint)){
#   age=agesint[i]
#   ten.year[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
#   lifetime[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
# }
# 
# 
# ages=40:80
# 
# s2=stateriskfunc_smoking_smoothedcoef(ages,prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
# p2=multipleprsfunc(s2[,,,1],prsprobs = pnorm(prs_quants))
# m=matriskfunc_coef(p2,ages,quantiles = prs_quants)
# agesint=seq(40,70,by=1)
# ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
# lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
# 
# for(i in 1:length(agesint)){
#   age=agesint[i]
#   ten.year[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
#   lifetime[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
# }
# 
# 
# ## How to use
# 
# s=stateriskfunc_smoking(ages = c(20:80),prs_quants = prs_quants,start = "HyperLip",
#                         stop = "Cad",modelfit = fixedsmoke)
# 
# 
# p=multipleprsfunc(s = s[,,,1],prsprobs = pnorm(prs_quants))
# ages=c(40:80)
# m=matriskfun(p,ages,quantiles = prs_quants)
# agesint=seq(40,70,by=1)
# ten.yearhl=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
# lifetimehl=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
# 
# for(i in 1:length(agesint)){
#   age=agesint[i]
#   ten.yearhl[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
#   lifetimehl[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
# }
# 
# ten.yearhl=data.frame(ten.yearhl)
# lifetimehl=data.frame(lifetimehl)
# 
# rownames(ten.yearhl)=agesint
# colnames(ten.yearhl)=c(1:10)
# ten.yearhl$age=agesint
# 
# rownames(lifetimehl)=agesint
# colnames(lifetimehl)=c(1:10)
# lifetimehl$age=agesint
# 
# 
# 
# lookup_table <- data.frame(melt(ten.yearhl,id.vars = c("age")))
# names(lookup_table)[3]="ten.yearhls"
# ggplot(lookup_table,aes(age,y = ten.yearhls,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")
# 
# lookup_table2 <- data.frame(melt(lifetimehl,id.vars = c("age")))
# names(lookup_table2)[3]="lifetimehls"
# 
# ggplot(lookup_table2,aes(age,y = lifetimehls,color=as.factor(variable)))+
#   stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")
# 
# ###################
# ## DO it with coef method
# ###
# ages=c(40:80)
# 
# s2=stateriskfunc_smoking_smoothedcoef(ages,prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
# p2=multipleprsfunc(s2[,,,1],prsprobs = pnorm(prs_quants))
# m=matriskfunc_coef(p2,ages,quantiles = prs_quants)
# agesint=seq(40,70,by=1)
# ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
# lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
# 
# for(i in 1:length(agesint)){
#   age=agesint[i]
#   ten.year[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
#   lifetime[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
# }
# 
# ten.year=data.frame(ten.year)
# lifetime=data.frame(lifetime)
# 
# rownames(ten.year)=agesint
# colnames(ten.year)=c(1:10)
# ten.year$age=agesint
# 
# rownames(lifetime)=agesint
# colnames(lifetime)=c(1:10)
# lifetime$age=agesint
# 
# 
# 
# lookup_table <- data.frame(melt(ten.year,id.vars = c("age")))
# names(lookup_table)[3]="ten.year"
# ggplot(lookup_table,aes(age,y = ten.year,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")
# 
# lookup_table2 <- data.frame(melt(lifetime,id.vars = c("age")))
# names(lookup_table2)[3]="lifetime"
# 
# g=ggplot(lookup_table2,aes(age,y = lifetime,color=as.factor(variable)))+
#   stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")
# 
# make sure that atrisk has the dimension of regression

coefinput=function(start,stop,ages,modelfit){
  sl=coefsmooth(start = start,stop = stop,ages = ages,modelfit = modelfit)
  mat=convertlistmat(sl)
  return(mat)
}

compute_prediction_product_matrix=function(atrisk,agepredinterval,coefmat){
#require(dplyr) 
atrisk_for_regression=data.frame(atrisk)
## ensure that every column is numeric
atriskmat=sapply(1:ncol(atrisk_for_regression),function(x){as.numeric(atrisk_for_regression[,x])})

##YxP
age_rows <- as.matrix(coefmat[rownames(coefmat) %in% agepredinterval, ])
###PxY
tar=t(age_rows)
### now calc X %*% beta which should b NxY years (i.e., NxP * PXY)
logoddsall=atriskmat%*%tar

prediction=exp(logoddsall)/(1+exp(logoddsall))
prediction_not=1-prediction
# Compute the product of the predictions which is over the colums since XB is NY
prediction_product <- apply(prediction_not,1,prod)
risk=1-prediction_product
return(list("PredictedIntervalrisk"=risk,"Survival"=prediction_not,"Yearly Risk"=prediction))}


## return smoothedmatrix
return_smoothedmatrix=function(start,stop,ages,modelfit){
smoothedlist=coefsmooth(start = start,stop = stop,ages = ages,modelfit = fixedsmoke)
mat=convertlistmat(smoothedlist = smoothedlist)
return(mat)}



###

tenlifeplotting=function(start,stop,modelfit,agesmooth,agepred,prs_quants,agesint){
  

mat=return_smoothedmatrix(start = start,stop = stop,ages = agesmooth,modelfit = modelfit)
atrisk=data.frame(cad.prs=rep(prs_quants,2),f.31.0.0=c(rep(0,length(prs_quants)),rep(1,length(prs_quants))),smoke=rep(0,length(prs_quants)))
ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.year[i,]=compute_prediction_product_matrix(atrisk = atrisk,agepredinterval = c(age:(age+10)),coefmat = mat)
  
  lifetime[i,]=compute_prediction_product_matrix(atrisk = atrisk,agepredinterval = c(age:(80)),coefmat = mat)
}

ten.year=data.frame(ten.year)
lifetime=data.frame(lifetime)

e=length(prs_quants)
rownames(ten.year)=agesint
colnames(ten.year)=paste0(rep(round(prs_quants,1),2),":",rep(c("female","male"),each=e))
ten.year$age=agesint

rownames(lifetime)=agesint
colnames(lifetime)=paste0(rep(round(prs_quants,1),2),":",rep(c("female","male"),each=e))
lifetime$age=agesint



lookup_table <- data.frame(melt(ten.year,id.vars = c("age")))
names(lookup_table)[3]="ten.year"
gten=ggplot(lookup_table,aes(age,y = ten.year,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y=paste0("Ten Year Risk from ",start," to ",stop),col="PRS:Sex")

lookup_table2 <- data.frame(melt(lifetime,id.vars = c("age")))
names(lookup_table2)[3]="lifetime"

glife=ggplot(lookup_table2,aes(age,y = lifetime,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y=paste0("Lifetme Risk from ",start," to ",stop),col="PRS:Sex")
return(list("tenplot"=gten,"lifeplot"=glife))}


# ### to apply to an entire data.frame ###
# sapply(seq(1:nrow(test)),function(x){compute_prediction_product_matrix(mpce[x,],agepredinterval = c(mpce[x,round(phenos.enrollment,0)]:(mpce[x,round(phenos.enrollment,0)+10])),coefmat = mat)})
