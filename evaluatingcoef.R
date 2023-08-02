#with coef smooth



mpce=readRDS("output/mpcecomplete.rds")
ages=c(40:80)
nstates=c("Health", "Ht","HyperLip","Dm","Cad","death","Ht&HyperLip","HyperLip&Dm","Ht&Dm","Ht&HyperLip&Dm")
#fixedsmoke = fitfunc(df_frame = mpce,nstates = nstates,ages = ages,mode = "binomial",covariates ="cad.prs+f.31.0.0+smoke")
# saveRDS(fixedsmoke,"output/fixedsmokefullmpce.rds")
fixedsmoke=readRDS("output/fixedsmokefullmpce.rds")
empiric.quants=data.frame(mpce%>%group_by(newint)%>%summarize(median(cad.prs),mean(cad.prs)))[c(2,4,6,8,10,12,14,16,18,20),3]

## do the mean
prs_quants=c(data.frame(mpce%>%group_by(newint)%>%summarize(median(cad.prs),mean(cad.prs)))[c(2,4,6,8,10,12,14,16,18,20),3])

# snew=stateriskfunc_smoking(ages=c(40:80),prs_quants = prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
# 
# a=multipleprsfunc(s = snew[,,,1],prsprobs = pnorm(prs_quants))
# projection_with_plot(a,40:80,quantiles = prs_quants,agestart = 40,agestop = 80)

snew=stateriskfunc_smoking_smoothedcoef(ages=c(40:80),prs_quants = prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)

a=multipleprsfunc(s = snew[,,,1],prsprobs = pnorm(prs_quants))
projection_with_plotcoef(a,40:80,quantiles = prs_quants,agestart = 40,agestop = 80)


prsprobs=pnorm(prs_quants)
#a=multipleprsfunc(s = snew[,,,2],prsprobs = prsprobs)
m=matriskfunc_coef(a,ages=c(40:80),quantiles = prs_quants)


set.seed(456)
enrollments=c(41:70)
aucmat=matrix(NA,nrow=length(enrollments),ncol=4)
prcmat=matrix(NA,nrow=length(enrollments),ncol=3)
ages=40:80
fixedsmoke=readRDS("~/multistate/output/fixedsmoke.rds")
enrollments=c(41:70)
## return a matrix of coefficients over all ages for a given state to state transition
mat=return_smoothedmatrix(start = "Health",stop = "Cad",ages = ages,modelfit = fixedsmoke)
mpce=readRDS("~/multistate/output/mpcecomplete.rds")

for(z in 1:length(enrollments)){
  age=enrollments[z]
  start=age
  stop=80
  
  df_frame=mpce
  atrisk = df_frame[age < Cad_0_censor_age &
                      age < Ht_0_censor_age &
                      age < HyperLip_0_censor_age &
                      age < Dm_0_censor_age &smoke==0 , ]
  
 
  df_updated=atrisk
  
  df_updated$ms=compute_prediction_product_matrix(coefmat = mat,atrisk = atrisk,agepredinterval = c(start:stop))
  ### return matrix of smoothed coefficeints
  #library(purrr)

  rm(atrisk)

  #require(pROC)
  df_updated$outcome=ifelse(df_updated$Cad_0_Any==2&df_updated$Cad_0_censor_age<stop,1,0)
  d=df_updated[round(phenos.enrollment,0)==age&smoke==0,]
  
  #d=d[!is.na(d$ascvd_10y_accaha),]
  aucmat[z,1]=roc(d$outcome~d$ms)$auc
  aucmat[z,2]=roc(d$outcome~d$ascvd_10y_accaha)$auc
  aucmat[z,4]=roc(d$outcome~d$cad.prs)$auc
  d=d[!is.na(d$ascvd_10y_accaha),]
  h=d$ascvd_10y_accaha+d$ms
  aucmat[z,3]=roc(d$outcome~h)$auc
  
  
  
  require(PRROC)
  fg <- d$ms[d$outcome == 1]
  bg <- d$ms[d$outcome == 0]
  
  prcmat[z,1]=pr.curve(scores.class0 = fg, scores.class1 = bg)$auc.integral
  
  #semat[i,1]=roc(d$outcome~d$ms)$se
  
  #require(PRROC)
  fg <- na.omit(d$ascvd_10y_accaha[d$outcome == 1])
  bg <- na.omit(d$ascvd_10y_accaha[d$outcome == 0])
  prcmat[z,2]=pr.curve(scores.class0 = fg, scores.class1 = bg)$auc.integral
  
  fg <- na.omit(d$cad.prs[d$outcome == 1])
  bg <- na.omit(d$cad.prs[d$outcome == 0])
  prcmat[z,3]=pr.curve(scores.class0 = fg, scores.class1 = bg)$auc.integral
  
  print(paste0("Completedforage",age))
  
  
  
}



improv=mean(aucmat[,1]-aucmat[,2])*100
rownames(aucmat)=enrollments
m=melt(aucmat,id.vars="Age")
names(m)=c("Age","Model","AUC")
m$Model=as.factor(m$Model)

levels(m$Model)[1]="MSGene"
levels(m$Model)[2]="PCE"
levels(m$Model)[3]="Combined"
levels(m$Model)[4]="PRS"
m=m[m$Model%in%c("MSGene","PCE","PRS"),]
aucplot <- ggplot(m,aes(x = Age,y = AUC,color = Model,ymin=AUC,ymax=AUC))+geom_point()+geom_line(aes(group=Model,color =Model),linewidth=3)+geom_pointrange()+ylim(0.5,1)+theme_classic()+ggtitle(paste0("10-year risk prediction"))


improv=mean(prcmat[,1]-prcmat[,2])*100
rownames(prcmat)=enrollments
m=melt(prcmat,id.vars="Age")
names(m)=c("Age","Model","PRauc")
m$Model=as.factor(m$Model)
levels(m$Model)[1]="MSGene"
levels(m$Model)[2]="PCE"
levels(m$Model)[3]="Combined"
levels(m$Model)[4]="PRS"
prplot <- ggplot(m,aes(x = Age,y = PRauc,color = Model,ymin=PRauc,ymax=PRauc))+geom_point()+geom_line(aes(group=Model,color =Model),linewidth=3)+geom_pointrange()+ylim(0,0.5)+theme_classic()+ggtitle(paste0("10-year risk prediction"))

ga=ggarrange(aucplot,prplot,common.legend = T)

ggsave(ga,file="plots/png/aucprc.png",height = 5,width=10)


# now do a massive RMSE to compute average prediction across strata

```{r}
mpce$cad.prs=scale(mpce$cad.prs)
mpce$cad.prs.lec=cut(mpce$cad.prs,breaks = c(-5,-0.84,0.84,5),labels = c("low","mid","high"))
mpce$int=interaction(mpce$f.31.0.0,mpce$cad.prs.lec)
levels(mpce$int) <- c(1,2,3,4,5,6)

mpce%>%group_by(int)%>%summarise(mean(cad.prs),median(cad.prs))
## grab means

prs_quants=c(data.frame(mpce%>%group_by(int)%>%summarize(median(cad.prs),mean(cad.prs)))[c(2,4,6),3])
prsprobs= pnorm(prs_quants)


## generate new arrays

# s=stateriskfunc_smoking(ages = c(40:80),prs_quants = prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
# p=multipleprsfunc(s = s[,,,1],prsprobs = pnorm(prs_quants))
# m=matriskfun(p,ages,quantiles = prs_quants)
# 
ages = c(40:80)
s2=stateriskfunc_smoking_smoothedcoef(ages = c(40:80),prs_quants = prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
p=multipleprsfunc(s = s2[,,,1],prsprobs = pnorm(prs_quants))
m=matriskfunc_coef(p,ages,quantiles = prs_quants)




agesint=seq(40,70,by=5)
ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.year[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
  lifetime[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
}


emp.ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
emp.lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length((agesint))){
  age=agesint[i]
  for(j in 1:length(levels(mpce$int))){
    
    cat=levels(mpce$int)[j]
    atrisk = mpce[age < Cad_0_censor_age &
                        age < Ht_0_censor_age &
                        age < HyperLip_0_censor_age &
                        age < Dm_0_censor_age &smoke==0 , ]
    
    emp.ten.year[i,j]=compute_empiricalrisk(age=age,age2 = age+10,df_frame = mpce,cat = cat)
    emp.lifetime[i,j]=compute_empiricalrisk(age=age,age2 = 100,df_frame = mpce,cat = cat)
  }}



ascvd.ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length((agesint))){
  age=agesint[i]
  for(j in 1:length(levels(mpce$int))){
    cat=levels(mpce$int)[j]
   atrisk=mpce[age < Cad_0_censor_age &
                               age < Ht_0_censor_age &
                              age < HyperLip_0_censor_age &
                              age < Dm_0_censor_age &smoke==0 , ]
    ascvd.ten.year[i,j]=compute_pce_predictedrisk(age=age,df_frame =mpce,cat = cat)
    #ascvdriskmat[i,2]=compute_empiricalrisk(age=40,df_frame = mpce,cat = cat)
  }
}





library(ehaGoF)

pcten=gofRMSE(Obs = as.vector(emp.ten.year*100),Prd = as.vector(ascvd.ten.year))
msten=gofRMSE(Obs = as.vector(emp.ten.year*100),Prd = as.vector(ten.year*100))
mslif=gofRMSE(Obs = as.vector(emp.lifetime*100),Prd = as.vector(lifetime*100))

pcr=gofRRMSE(Obs = as.vector(emp.ten.year),Prd = as.vector(ascvd.ten.year))
msr=gofRRMSE(Obs = as.vector(emp.ten.year),Prd = as.vector(ten.year))
mslifr=gofRRMSE(Obs = as.vector(emp.lifetime),Prd = as.vector(lifetime))

print(xtable(data.frame("RMSE"=c(pcten,msten,mslif),"RRMSE"=c(pcr,msr,mslifr))),row.names=FALSE)

diff.ascvd=abs(data.frame(ascvd.ten.year/100-emp.ten.year))
d=as.matrix(diff.ascvd)
sqrt(mean(d^2))

diff.mstate=abs(data.frame(ten.year-emp.ten.year))
d=as.matrix(diff.mstate)
sqrt(mean(d^2))

diff.mstate.life=abs(data.frame(lifetime-emp.lifetime))
d=as.matrix(diff.mstate.life)
sqrt(mean(d^2))

# mean(as.matrix(emp.ten.year))
# mean.es=mean(as.matrix(emp.ten.year))
# rel.diff.mstate=sqrt(mean(d^2))/mean.es
# 
# diff.mstate=data.frame(lifetime-emp.lifetime)

# ##  sqrt(mean(as.matrix(diff.mstate)))
# [1] 0.1522733
# 
# sd(as.matrix(diff.mstate^2))
# Error: unexpected ')' in "sd(as.matrix(diff.mstate)))"
# # > sd(as.matrix(diff.mstate))
# # [1] 0.01173952


# rdiff=diff.mstate/emp.lifetime
# d=as.matrix(rdiff)
# sqrt(mean(d^2))/mean(as.matrix(emp.lifetime))
# 
# rel.diff.mstate=sqrt(mean((diff.mstate/(emp.lifetime*100))^2))

diff.ascvd$se=sd(as.matrix(sqrt(diff.ascvd^2)))
diff.ascvd$score=rep("PCE",7)
diff.ascvd$age=agesint

diff.mstate$se=sd(as.matrix(sqrt(diff.mstate^2)))
diff.mstate$score=rep("MSGene",7)
diff.mstate$age=agesint
r=rbind(diff.ascvd,diff.mstate)
rownames(r)=NULL
rf=r[,c(1,3,5,7,8,9)]
rf$sex=rep("female",nrow(rf))

rm=r[,c(2,4,6,7,8,9)]
rm$sex=rep("male",nrow(rm))
names(rf)[1:3]=names(rm)[1:3]=c("low","medium","high")



r=rbind(diff.ascvd,diff.mstate)
rownames(r)=NULL
rf=r[,c(1,3,5,7,8,9)]
rm=r[,c(2,4,6,7,8,9)]

#t.test(x = rm[c(1:7),c(1:3)],r[c(8:14),c(1:3)])

colnames(rm)=c("Low","Intermediate","High","se","score","age")
m=melt(rm,id.vars=c("age","score","se"))

m$se=m$se/sqrt(1000)
m$interaction=interaction(m$variable,m$score)
m$value=abs(m$value)
#interaction_colors=c(brewer.pal(n = 6, name = "RdBu"))
interaction_colors <- c(brewer.pal(n = 3, name = "Reds")[1:3], brewer.pal(n = 3, name = "Blues"))

r2_male=ggplot(data = m,
               aes(x=age,
                   y= value,
                   ymin=value-se,
                   ymax=value+se,
                   fill=interaction)) +scale_fill_manual(values=interaction_colors)+
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar( position = position_dodge(), colour="black") +labs(y="RMSE 10 year risk",x="Age",fill="Genomic Level: Score")+theme_classic(base_size = 20)+ylim(c(0,0.20))
#geom_point(position=position_dodge(.9), aes(y=value, colour=interaction))


ggsave(r2_male,file = "~/multistate/rmseten.png",dpi = 300,height = 4,width = 8)

## to do use mpce long
### do for female


colnames(rf)=c("Low","Intermediate","High","se","score","age")
m=melt(rf,id.vars=c("age","score","se"))
m$se=m$se/sqrt(1000)
m$interaction=interaction(m$variable,m$score)
m$value=abs(m$value)
#interaction_colors=c(brewer.pal(n = 6, name = "RdBu"))


r2_female=ggplot(data = m,
                 aes(x=age,
                     y= value,
                     ymin=value-se,
                     ymax=value+se,
                     fill=interaction)) +scale_fill_manual(values=interaction_colors)+
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar( position = position_dodge(), colour="black") +ylim(c(0,0.20))+labs(y="RMSE 10 year risk, Female",x="Age",fill="Genomic Level: Score")+theme_classic(base_size = 15)

#geom_point(position=position_dodge(.9), aes(y=value, colour=interaction))


ggsave(ggarrange(r2_male,r2_female,nrow = 1,common.legend = T,legend = "right"),file = "~/multistate/rmseten_both.png",dpi = 300,height = 4,width = 12)


############

r=rbind(rf,rm)
r$sex=c(rep("female",length(agesint),rep("male",length(agesint))))
m=melt(r,id.vars=c("age","sex","se","score"))
m$interaction=interaction(m$variable,m$sex)
m$se=c(rf$se,rm$se)
m$se=m$se/sqrt(1000)
m$value=as.numeric(m$value)
#interaction_colors=c(brewer.pal(n = 6, name = "RdBu"))

r2=ggplot(data = m,
          aes(x=age,
              y= value,
              ymin=value-se,
              ymax=value+se,
              fill=interaction)) +scale_fill_manual(values=interaction_colors)+
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar( position = position_dodge(), colour="black") 
labs(y="RMSE 10 year risk, Health to CAD",x="Age",fill="Genomic Level: Score")+theme_classic(base_size = 15)
#geom_point(position=position_dodge(.9), aes(y=value, colour=interaction))


diff.mstate.life$se=sd(as.matrix(sqrt(diff.mstate.life^2)))
diff.mstate.life$score=rep("mstate",7)
diff.mstate.life$age=agesint
lf=diff.mstate.life[,c(1,3,5,7,8,9)]
lf$sex=rep("female",nrow(lf))
lm=diff.mstate.life[,c(2,4,6,7,8,9)]
lm$sex=rep("male",nrow(lm))
names(lf)[1:3]=names(lm)[1:3]=c("low","medium","high")
l2=rbind(lf,lm)

l2$sex=as.factor(l2$sex)
m=melt(l2,id.vars=c("age","sex","se","score"))
m$interaction=interaction(m$variable,m$sex)
m$se=m$se/sqrt(1000)
m$value=as.numeric(m$value)
#interaction_colors=c(brewer.pal(n = 6, name = "RdBu"))

r2=ggplot(data = m,
          aes(x=age,
              y= value,
              ymin=value-se,
              ymax=value+se,
              fill=interaction)) +scale_fill_manual(values=interaction_colors)+
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar( position = position_dodge(), colour="black") +
  labs(y="RMSE Lifetime risk, Health to CAD",x="Age",fill="Genomic Level: Sex")+theme_classic(base_size = 15)
#geom_point(position=position_dodge(.9), aes(y=value, colour=interaction))

ggsave(r2,file = "~/multistate/rmselife.png",dpi = 300,height = 4,width = 8)
```

###
## return a matrix of coefficients over all ages for a given state to state transition
mat=return_smoothedmatrix(start = "Health",stop = "Cad",ages = ages,modelfit = fixedsmoke)

enrollments=c(41:70)
thresh_l=10
thresh_t=5

differencelist=list()
overa_list=list()


saved=matrix(NA,nrow=length(enrollments),ncol=7)
dfascvd=readRDS("~/multistate/output/dfascvd_newbp.rds")
m2=merge(mpce,dfascvd[,c("sample_id","sbp","hdladj","choladj","Race","as2")],by.x="identifier",by.y="sample_id")
differencelist=list()
overa_list=list()

for(z in 1:length(enrollments)){
  age=enrollments[z]
  start=age
  stop=80
   
  df_frame=m2
  

  
  atrisk = df_frame[age < Cad_0_censor_age &
                      age < Ht_0_censor_age &
                      age < HyperLip_0_censor_age &
                      age < Dm_0_censor_age &smoke==0 , ]
  atrisk$msten=100*compute_prediction_product_matrix(atrisk,agepredinterval = c(start:(start+10)),coefmat = mat)
  atrisk$mslife=100*compute_prediction_product_matrix(atrisk,agepredinterval = c(start:stop),coefmat = mat)
    
  df_updated=atrisk
  #require(pROC)
  df_updated$outcome=ifelse(df_updated$Cad_0_Any==2,1,0)
  d=df_updated
  d$Race=ifelse(d$Race=="white",0,1)
  #d=df_updated[round(phenos.enrollment,0)==age&smoke==0,]## here do not restrict to those at a particular age of enrolle
  #d=d[!is.na(d$ascvd_10y_accaha),]
  
  print(age)
  # d$ldladj=scale(d$ldladj)
  # d$hdladj=scale(d$hdladj)
  # d$choladj=scale(d$choladj)
  # d$sbp=scale(d$sbp)
  # d$dbp=scale(d$dbp)
  over=d[d$msten>thresh_t,]
  
  coolmat=over[,c("cad.prs","as2","hdladj","choladj","sbp","Race","f.31.0.0","phenos.enrollment")]
  coolmat$f.31.0.0=as.numeric(ifelse(coolmat$f.31.0.0=="1",1,0))
  overa_list[[z]]=coolmat
  under=d[d$msten<thresh_t&d$mslife>thresh_l,]
  undermat=under[,c("cad.prs","as2","hdladj","choladj","sbp","Race","f.31.0.0","phenos.enrollment")]
  undermat$f.31.0.0=as.numeric(ifelse(undermat$f.31.0.0=="1",1,0))
  differencelist[[z]]=undermat
  under_low=mean(d$msten<thresh_t&d$mslife>thresh_l&d$cad.prs.lev=="low")
  under_mid=mean(d$msten<thresh_t&d$mslife>thresh_l&d$cad.prs.lev=="mid")
  under_high=mean(d$msten<thresh_t&d$mslife>thresh_l&d$cad.prs.lev=="high")
  saved[z,1]=sum(d$msten<thresh_t&d$mslife>thresh_l)
  saved[z,2]=mean(d$msten<thresh_t&d$mslife>thresh_l)
  saved[z,3]=age
  saved[z,4]=nrow(d)
  saved[z,5]=under_low
  saved[z,6]=under_mid
  saved[z,7]=under_high
  
  
}

## at age 40, 48% (SEM 0.005) of people have a lifetime risk over 10% but an ascvd risk <5%, whil at age 60 only 2.6% (SEM 0.0015) do   saved[,1]/saved[,4]


c=cbind(enrollments,saved[,1]/saved[,4])
sd=sqrt(((1-saved[,1]/saved[,4])*(saved[,1]/saved[,4]))/saved[,5])

overmat=t(sapply(overa_list,function(x){
  colMeans(x)
}))

rownames(overmat)=rownames(saved)=enrollments
overmat=data.frame(overmat)
overmat$f=rep("over",nrow(overmat))

meanmat=t(sapply(differencelist,function(x){
  colMeans(x)
}))
rownames(meanmat)=enrollments
meanmat=data.frame(meanmat)
meanmat$f=rep("under",nrow(meanmat))

mass=rbind(overmat,meanmat)
colnames(mass)=c("CAD.prs","PCE 10 year","HDL-C","Total Cholesterol","Systolic BP","Race","Proportion Male","Age at Enroll","Score")
mass$age=enrollments
ma=melt(mass,id.vars=c("age","Score"))
#ma=ma[ma$variable!="PCE 10 year",]
#ma=ma[ma$variable!="Total Cholesterol",]
ma$Score=factor(ma$Score,levels = c("over","under"),labels = c("MSTen > 5%","MSLife > 10% & MSTen < 5%"))
wrap=ggplot(ma,aes(age,y=value,group=Score,color=Score))+stat_smooth()+labs(y="Value at Enrollment",x="Age at Calc",)+facet_wrap(~variable,scales="free")+theme_classic2()

####
###Now low ascvd but don't restrict

saved=matrix(NA,nrow=length(enrollments),ncol=7)
dfascvd=readRDS("~/multistate/output/dfascvd_newbp.rds")
m2=merge(mpce,dfascvd[,c("sample_id","sbp","hdladj","choladj","Race","as2")],by.x="identifier",by.y="sample_id")
differencelist=list()
overa_list=list()

for(z in 1:length(enrollments)){
  age=enrollments[z]
  start=age
  stop=80
  
  df_frame=m2
  
  
  
  atrisk = df_frame[age < Cad_0_censor_age &
                      age < Ht_0_censor_age &
                      age < HyperLip_0_censor_age &
                      age < Dm_0_censor_age &smoke==0 , ]
  atrisk$msten=100*compute_prediction_product_matrix(atrisk,agepredinterval = c(start:(start+10)),coefmat = mat)
  atrisk$mslife=100*compute_prediction_product_matrix(atrisk,agepredinterval = c(start:stop),coefmat = mat)
  
  df_updated=atrisk
  #require(pROC)
  df_updated$outcome=ifelse(df_updated$Cad_0_Any==2,1,0)
  ### now restrict to those a particular age at enrollment
  #d=df_updated[round(phenos.enrollment,0)==age&smoke==0,]
  d=df_updated
  d$Race=ifelse(d$Race=="white",0,1)
  
  d=d[!is.na(d$ascvd_10y_accaha),]
  
  print(age)
  # d$ldladj=scale(d$ldladj)
  # d$hdladj=scale(d$hdladj)
  # d$choladj=scale(d$choladj)
  # d$sbp=scale(d$sbp)
  # d$dbp=scale(d$dbp)
  over=d[d$ascvd_10y_accaha>thresh_t,]
  
  coolmat=over[,c("cad.prs","ascvd_10y_accaha","hdladj","choladj","sbp","Race","f.31.0.0","phenos.enrollment")]
  coolmat$f.31.0.0=as.numeric(ifelse(coolmat$f.31.0.0=="1",1,0))
  overa_list[[z]]=coolmat
  under=d[d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l,]
  undermat=under[,c("cad.prs","ascvd_10y_accaha","hdladj","choladj","sbp","Race","f.31.0.0","phenos.enrollment")]
  undermat$f.31.0.0=as.numeric(ifelse(undermat$f.31.0.0=="1",1,0))
  differencelist[[z]]=undermat
  under_low=mean(d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l&d$cad.prs.lev=="low")
  under_mid=mean(d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l&d$cad.prs.lev=="mid")
  under_high=mean(d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l&d$cad.prs.lev=="high")
  saved[z,1]=sum(d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l)
  saved[z,2]=mean(d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l)
  saved[z,3]=age
  saved[z,4]=nrow(d)
  saved[z,5]=under_low
  saved[z,6]=under_mid
  saved[z,7]=under_high
  
  
}

## at age 40, 48% (SEM 0.005) of people have a lifetime risk over 10% but an ascvd risk <5%, whil at age 60 only 2.6% (SEM 0.0015) do   saved[,1]/saved[,4]


c=cbind(enrollments,saved[,1]/saved[,4])
sd=sqrt(((1-saved[,1]/saved[,4])*(saved[,1]/saved[,4]))/saved[,5])

overmat=t(sapply(overa_list,function(x){
  colMeans(x)
}))

rownames(overmat)=rownames(saved)=enrollments
overmat=data.frame(overmat)
overmat$f=rep("over",nrow(overmat))

meanmat=t(sapply(differencelist,function(x){
  colMeans(x)
}))
rownames(meanmat)=enrollments
meanmat=data.frame(meanmat)
meanmat$f=rep("under",nrow(meanmat))

mass=rbind(overmat,meanmat)
colnames(mass)=c("CAD.prs","PCE 10 year","HDL-C","Total Cholesterol","Systolic BP","Race","Proportion Male","Age at Enrollment","Score")
mass$age=enrollments
ma=melt(mass,id.vars=c("age","Score"))
#ma=ma[ma$variable!="PCE 10 year",]
#ma=ma[ma$variable!="Total Cholesterol",]
ma$Score=factor(ma$Score,levels = c("over","under"),labels = c("PCE > 5%","MSLife > 10% & PCE < 5%"))
wrap=ggplot(ma,aes(age,y=value,group=Score,color=Score))+stat_smooth()+labs(y="Value at Enrollment",x="Age at Enrollment",)+facet_wrap(~variable,scales="free")+theme_classic2()



## now low ASCVD and restirction
mat=return_smoothedmatrix(start = "Health",stop = "Cad",ages = ages,modelfit = fixedsmoke)
saved=matrix(NA,nrow=length(enrollments),ncol=7)
dfascvd=readRDS("~/multistate/output/dfascvd_newbp.rds")
m2=merge(mpce,dfascvd[,c("sample_id","sbp","hdladj","choladj","Race","as2")],by.x="identifier",by.y="sample_id")
differencelist=list()
overa_list=list()

for(z in 1:length(enrollments)){
  age=enrollments[z]
  start=age
  stop=80
  
  df_frame=m2
  
  
  
  atrisk = df_frame[age < Cad_0_censor_age &
                      age < Ht_0_censor_age &
                      age < HyperLip_0_censor_age &
                      age < Dm_0_censor_age &smoke==0 , ]
  atrisk$msten=100*compute_prediction_product_matrix(atrisk,agepredinterval = c(start:(start+10)),coefmat = mat)
  atrisk$mslife=100*compute_prediction_product_matrix(atrisk,agepredinterval = c(start:stop),coefmat = mat)
  
  df_updated=atrisk
  #require(pROC)
  df_updated$outcome=ifelse(df_updated$Cad_0_Any==2,1,0)
  ### now restrict to those a particular age at enrollment
  d=df_updated[round(phenos.enrollment,0)==age&smoke==0,]

  d$Race=ifelse(d$Race=="white",0,1)
  
  d=d[!is.na(d$ascvd_10y_accaha),]
  
  print(age)
  # d$ldladj=scale(d$ldladj)
  # d$hdladj=scale(d$hdladj)
  # d$choladj=scale(d$choladj)
  # d$sbp=scale(d$sbp)
  # d$dbp=scale(d$dbp)
  over=d[d$ascvd_10y_accaha>thresh_t,]
  
  coolmat=over[,c("cad.prs","ascvd_10y_accaha","hdladj","choladj","sbp","Race","f.31.0.0","phenos.enrollment")]
  coolmat$f.31.0.0=as.numeric(ifelse(coolmat$f.31.0.0=="1",1,0))
  overa_list[[z]]=coolmat
  under=d[d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l,]
  undermat=under[,c("cad.prs","ascvd_10y_accaha","hdladj","choladj","sbp","Race","f.31.0.0","phenos.enrollment")]
  undermat$f.31.0.0=as.numeric(ifelse(undermat$f.31.0.0=="1",1,0))
  differencelist[[z]]=undermat
  under_low=mean(d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l&d$cad.prs.lev=="low")
  under_mid=mean(d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l&d$cad.prs.lev=="mid")
  under_high=mean(d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l&d$cad.prs.lev=="high")
  saved[z,1]=sum(d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l)
  saved[z,2]=mean(d$ascvd_10y_accaha<thresh_t&d$mslife>thresh_l)
  saved[z,3]=age
  saved[z,4]=nrow(d)
  saved[z,5]=under_low
  saved[z,6]=under_mid
  saved[z,7]=under_high
  
  
}

## at age 40, 48% (SEM 0.005) of people have a lifetime risk over 10% but an ascvd risk <5%, whil at age 60 only 2.6% (SEM 0.0015) do   saved[,1]/saved[,4]


c=cbind(enrollments,saved[,1]/saved[,4])
sd=sqrt(((1-saved[,1]/saved[,4])*(saved[,1]/saved[,4]))/saved[,5])

overmat=t(sapply(overa_list,function(x){
  colMeans(x)
}))

rownames(overmat)=rownames(saved)=enrollments
overmat=data.frame(overmat)
overmat$f=rep("over",nrow(overmat))

meanmat=t(sapply(differencelist,function(x){
  colMeans(x)
}))
rownames(meanmat)=enrollments
meanmat=data.frame(meanmat)
meanmat$f=rep("under",nrow(meanmat))

mass=rbind(overmat,meanmat)
colnames(mass)=c("CAD.prs","PCE 10 year","HDL-C","Total Cholesterol","Systolic BP","Race","Proportion Male","Age at Enrollment","Score")
mass$age=enrollments
ma=melt(mass,id.vars=c("age","Score"))
#ma=ma[ma$variable!="PCE 10 year",]
#ma=ma[ma$variable!="Total Cholesterol",]
ma$Score=factor(ma$Score,levels = c("over","under"),labels = c("PCE > 5%","MSLife > 10% & PCE < 5%"))
wrap=ggplot(ma,aes(age,y=value,group=Score,color=Score))+stat_smooth()+labs(y="Value at Enrollment",x="Age at Enrollment",)+facet_wrap(~variable,scales="free")+theme_classic2()


####
prs_quants=c(-5,-4,-3,-2,-1,0,1,2,3,4,5)
s2=stateriskfunc_smoking_smoothedcoef(ages,prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
p2=multipleprsfunc(s2[,,,1],prsprobs = pnorm(prs_quants))
m=matriskfunc_coef(p2,ages,quantiles = prs_quants)
agesint=seq(40,70,by=1)
ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.year[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
  lifetime[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
}

ten.year=data.frame(ten.year)
lifetime=data.frame(lifetime)

rownames(ten.year)=agesint
colnames(ten.year)=1:(length(prs_quants)*2)
ten.year$age=agesint

rownames(lifetime)=agesint
colnames(lifetime)=1:(length(prs_quants)*2)
lifetime$age=agesint



lookup_table <- data.frame(melt(ten.year,id.vars = c("age")))
names(lookup_table)[3]="ten.year"
ggplot(lookup_table,aes(age,y = ten.year,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")

lookup_table2 <- data.frame(melt(lifetime,id.vars = c("age")))
names(lookup_table2)[3]="lifetime"

g=ggplot(lookup_table2,aes(age,y = lifetime,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")

###now do with compute_prediction fun


mat=return_smoothedmatrix(start = "HyperLip",stop = "Cad",ages = c(40:80),modelfit = fixedsmoke)
atrisk=data.frame(cad.prs=rep(c(-5:5),2),f.31.0.0=c(rep(0,11),rep(1,11)),smoke=rep(0,22))

#compute_prediction_product_matrix(atrisk = atrisk,agepredinterval = 40:80,coefmat = mat)
agesint=seq(40,70,by=1)
ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.year[i,]=compute_prediction_product_matrix(atrisk = atrisk,agepredinterval = c(age:(age+10)),coefmat = mat)

  lifetime[i,]=compute_prediction_product_matrix(atrisk = atrisk,agepredinterval = c(age:(80)),coefmat = mat)
}

ten.year=data.frame(ten.year)
lifetime=data.frame(lifetime)

rownames(ten.year)=agesint
colnames(ten.year)=paste0(rep(prs_quants,2),":",rep(c("female","male"),each=11))
ten.year$age=agesint

rownames(lifetime)=agesint
colnames(lifetime)=paste0(rep(prs_quants,2),":",rep(c("female","male"),each=11))
lifetime$age=agesint



lookup_table <- data.frame(melt(ten.year,id.vars = c("age")))
names(lookup_table)[3]="ten.year"
ggplot(lookup_table,aes(age,y = ten.year,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")

lookup_table2 <- data.frame(melt(lifetime,id.vars = c("age")))
names(lookup_table2)[3]="lifetime"

g=ggplot(lookup_table2,aes(age,y = lifetime,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")


###

tenlifeplotting(start = "HyperLip",stop = "Cad",modelfit = fixedsmoke,agesmooth = 40:80,agepred = 40:80,prs_quants = prs_quants,agesint = 41:70)

phenos=c("Health","Ht","Dm","HyperLip")
atrisk=data.frame(cad.prs=rep(c(-5:5),2),f.31.0.0=c(rep(0,11),rep(1,11)),smoke=rep(0,22))
ten.year=matrix(NA,nrow = length(phenos),ncol=length(prs_quants)*2)

for(i in 1:length("Health","Ht","Dm","HyperLip")){
  
mat=return_smoothedmatrix(start = x,stop = "Cad",ages = ages,modelfit = fixedsmoke)
compute_prediction_product_matrix(atrisk = atrisk,agepredinterval = c(40:50),coefmat = mat)

#compute_prediction_product_matrix(atrisk = atrisk,agepredinterval = 40:80,coefmat = mat)
agesint=seq(40,70,by=1)
ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.year[i,]=compute_prediction_product_matrix(atrisk = atrisk,agepredinterval = c(age:(age+10)),coefmat = mat)
  
