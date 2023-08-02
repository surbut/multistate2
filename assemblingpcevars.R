##using the results from Denaxas at el, get BP readings

#https://github.com/spiros/ukb-biomarker-phenotypes/tree/master#implementation-8

# IF data_provider = England Vision (1)
# IF read_2 = “246..00 O/E - blood pressure reading”' 
#         SBP = value1
#         DBP = value2

gp_event


engV=gp_event[gp_event$data_provider==1 & gp_event$read_2=="246..",c("eid","read_2","read_3","event_dt","value1")]


# ELSE IF data_provider = Scotland (2)
#     IF read_2 = “246..00 O/E - blood pressure reading”
#         DBP = value1
#         SBP = value2

s=gp_event[gp_event$data_provider==2 & gp_event$read_2=="246..",c("eid","read_2","read_3","event_dt","value2")]

#     ELSE IF read_2 = “2469.00 O/E - Systolic BP reading”
#         SBP = value1

s2=gp_event[gp_event$data_provider==2 & gp_event$read_2=="2469.",c("eid","read_2","read_3","event_dt","value1")]

#     ELSE IF read_2 = “246A.00 O/E - Diastolic BP reading”
#         DBP = value1
# 
# IF data_provider = England TPP (3)
#     IF read_3 = “2469.00 O/E - Systolic BP reading” 
#         SBP = value1

et=gp_event[gp_event$data_provider==3&gp_event$read_3=="2469.",c("eid","read_2","read_3","event_dt","value1")]

#     ELSE IF read_3 = “246A.00 O/E - Diastolic BP reading” 
#         DBP = value1


# 
# IF data_provider = Wales (4)
#     IF read_2 = “246..00 O/E - blood pressure reading”
#         SBP = value1

ew1=gp_event[gp_event$data_provider==4 & gp_event$read_2=="246..",c("eid","read_2","read_3","event_dt","value1")]

#         DBP = value2
#     ELSE IF read_2 = “2469.00 O/E - Systolic BP reading”
#         SBP = value1

ew=gp_event[gp_event$data_provider==4 & gp_event$read_2=="2469.",c("eid","read_2","read_3","event_dt","value1")]

#     ELSE IF read_2 = “246A.00 O/E - Diastolic BP reading”
#         DBP = value1


sbp=rbind(rbind(rbind(rbind(rbind(engV,s,use.names=FALSE),s2),et),ew1),ew)
sbp$value1=as.numeric(sbp$value1)
sbp=sbp[sbp$value1>50&sbp$value1<220,]

saveRDS(sbp,"~/Dropbox/denaxasvalues/sbp.rds")

### now do for HDL from Denaxas et al Supplementary Figure 1: Flowchart description of the main steps involved in data extraction for each biomarker. In England and Wales data sources, measurements are recorded in the value1 field whereas in Scotland measurements are recoded in value2 and units in value3.
# HDL	read2	44P5.00	Serum HDL cholesterol level
# HDL	read2	44PB.00	Serum fasting HDL cholesterol level
# HDL	read2	44PC.00	Serum random HDL cholesterol level
# HDL	ctv3	44P5.	Serum HDL cholesterol level
# HDL	ctv3	44PB.	Serum fasting HDL cholesterol level
# HDL	ctv3	44PC.	Serum random HDL cholesterol level

HDL=gp_clinical[gp_clinical$read_2%in%c("44P5.","44PB.","44PC.","44P5.","44PB.","44PC.")|gp_clinical$read_3%in%c("44P5.","44PB.","44PC."),]
HDL$value1=as.numeric(HDL$value1)*38.66976

HDL=HDL[!is.na(HDL$value1),]
HDL=HDL[HDL$value1>10&HDL$value1<200,]
saveRDS(HDL,"~/Dropbox/denaxasvalues/hdl.rds")

# ####
# Cholesterol	read2	44P..00	Serum cholesterol
# Cholesterol	read2	44PJ.00	Serum total cholesterol level
# Cholesterol	read2	44PH.00	Total cholesterol measurement
# Cholesterol	read2	44P3.00	Serum cholesterol raised
# Cholesterol	read2	44P1.00	Serum cholesterol normal
# Cholesterol	read2	44PZ.00	Serum cholesterol NOS
# Cholesterol	read2	44P9.00	Serum cholesterol studies
# Cholesterol	read2	44P2.00	Serum cholesterol borderline
# Cholesterol	read2	44PK.00	Serum fasting total cholesterol
# Cholesterol	read2	44P4.00	Serum cholesterol very high
# Cholesterol	ctv3	XE2eD	Serum cholesterol
# Cholesterol	ctv3	XE2eD	Serum cholesterol level
# Cholesterol	ctv3	XaJe9	Serum total cholesterol level
# Cholesterol	ctv3	XSK14	Total cholesterol measurement
# Cholesterol	ctv3	44P3.	Serum cholesterol raised
# Cholesterol	ctv3	44P1.	Serum cholesterol normal
# Cholesterol	ctv3	44PZ.	Serum cholesterol NOS
# Cholesterol	ctv3	44P9.	Serum cholesterol studies
# Cholesterol	ctv3	44P2.	Serum cholesterol borderline
# Cholesterol	ctv3	XaLux	Serum fasting total cholesterol
# Cholesterol	ctv3	44P4.	Serum cholesterol very high


TC=gp_clinical[gp_clinical$read_2%in%c("44P..","44PJ.","44PH.","44P3.","44P1.","44PZ.","44P9.","44P2.","44PK.","44P4.")|
                 gp_clinical$read_3%in%c("XE2eD","XaJe9","XSK14","44P3.","44P1.","44PZ.","44P9.","44P2.","XaLux","44P4")]
TC$value1=as.numeric(TC$value1)*38.66976
TC=TC[!is.na(TC$value1),]

TC=TC[which(TC$value1>50&TC$value1<300),]

saveRDS(TC,"~/Dropbox/denaxasvalues/tc.rds")

###


i=intersect(intersect(sbp$eid,HDL$eid),TC$eid)

###
load("~/Dropbox (Personal)/pheno_dir/output/merged_pheno_censor_final_withdrugs.rds")

order_biomarker=function(biomarker_file){
  biomarker=biomarker_file
  m=merge(biomarker,dfh[,c("identifier","Birthdate")],by.x = "eid",by.y = "identifier")
  m2=m[with(m, order(event_dt)), ]
  m3=m2[!duplicated(m2$eid), ]
  m3$age=as.numeric(difftime(time1 = m3$event_dt,time2 = m3$Birthdate,units = "days")/365.25)
  m3=m3[m3$age>18,]
  
  ## show that biomarker often followed later even though first reading earlier
  m3$followup=as.numeric(difftime(time1 ="2021-03-31",time2 = m3$event_dt,units = "days")/365.25)
  return(m3)
 
}

HDLo=order_biomarker(HDL)
rm(HDL)

SBPo=order_biomarker(sbp)
rm(sbp)

TCo=order_biomarker(TC)
rm(TC)

### how many have greater than 20 years of follow up since first measurement
sum(SBPo$followup>20)
sum(HDLo$followup>20)
sum(TCo$followup>20)

### folks with all
m=merge(SBPo[,c("eid","event_dt","age","value1","followup")],
        merge(
TCo[,c("eid","event_dt","age","value1","followup")],HDLo[,c("eid","event_dt","age","value1","followup","Birthdate")],
            by="eid"),by="eid")
names(m)[4]="SBP"
names(m)[8]="TC"
names(m)[12]="HDL"

# Load packages
library(lcsm)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

m=data.frame(m)

m$maxfollow=apply(m,1,function(x){max(x["followup"],max(x["followup.x"],x["followup.y"]))})

m$minfollow=apply(m,1,function(x){min(x["followup"],min(x["followup.x"],x["followup.y"]))})
m=m[order(m$minfollow,decreasing=T),]

#`%out%` <- function(a,b) ! a %in% b
mlong=m

mlong$agecompute=apply(mlong,1,function(x){min(x["age"],min(x["age.x"],x["age.y"]))})

mlong=merge(mlong,dfh[,c("identifier","f.31.0.0","Dm_0_censor_age","Dm_0_Any","Cad_0_Any","Cad_0_censor_age",
                         "antihtn","htn_age","statin","statin_age")],by.x = "eid",by.y = "identifier")

mlong$agecompute=as.numeric(mlong$agecompute)


##biased sample
mean(mlong$antihtn==1)
mlong$htn_age=as.numeric(mlong$htn_age)
mlong$statin_age=as.numeric(mlong$statin_age)

##Only count med if started more than a year before
mlong[mlong$antihtn==1&(mlong$agecompute-mlong$htn_age)>1,"antihtn_compute"]=1
mlong[mlong$antihtn==1&(mlong$agecompute-mlong$htn_age)<1,"antihtn_compute"]=0
mlong$antihtn_compute[mlong$antihtn==0]=0 

mlong[mlong$statin==1&(mlong$agecompute-mlong$statin_age)>1,"statin_compute"]=1
mlong[mlong$statin==1&(mlong$agecompute-mlong$statin_age)<1,"statin_compute"]=0
mlong$statin_compute[mlong$statin==0]=0

##Only count DM if DM occured before compute time 

mlong[mlong$Dm_0_Any==2&(mlong$agecompute-mlong$Dm_0_censor_age)>1,"Dm_compute"]=1
mlong[mlong$Dm_0_Any==2&(mlong$agecompute-mlong$Dm_0_censor_age)<1,"Dm_compute"]=0
mlong$Dm_compute[mlong$Dm_0_Any==1]=0
saveRDS(mlong,"~/Dropbox/pheno_dir/mlong.rds")
### Biased sample

t.test(mshort$SBP,mlong$SBP)$conf.int
t.test(mshort$TC,mlong$TC)$conf.int
t.test(mshort$HDL,mlong$HDL)$conf.int

#####
ascvddf=readRDS("~/multistate/output/dfascvd_newbp.rds")
colnames(ascvddf)[1]="eid"
colnames(ascvddf)[4]="ageenroll"
mlong=readRDS("~/Dropbox/pheno_dir/mlong.rds")
dfukb_baseline=readRDS("~/Dropbox/pheno_dir/dfukbaseline_7.rds")

mlong=merge(mlong,ukb[,c("id","in_white_British_ancestry_subset")],by.x ="eid",by.y="id")
mlong=merge(mlong,dfukb_baseline[,c("identifier","f.20116.0.0")],by.x="eid",by.y="identifier",all.x = T)
mlong$Race=mlong$in_white_British_ancestry_subset
mlong$smoke=ifelse(mlong$f.20116.0.0==2,1,0)
mlong$smoke[is.na(mlong$smoke)]=0
mlong$gender=ifelse(mlong$f.31.0.0==0,"female","male")
mlong$Race[is.na(mlong$Race)]="white"
mlong$Race=ifelse(mlong$Race==1,"white","other")

## compute based on first SBP/HDL/TC
c=compute_CVrisk2(scores =c("as2"),age="agecompute",df = mlong,gender = "gender",sbp ="SBP",
                  hdl = "HDL",totchol = "TC",race = "Race",
                  bp_med = "antihtn_compute",lipid_med = "statin_compute",diabetes ="Dm_compute",
               smoker ="smoke",fh_heartattack = NULL, cac = NULL)

# 
# mp=mlong
# mp$bpmed1=ifelse(mp$antihtn_compute==1,"yes","no")
# mp$smoke2=ifelse(mp$smoke==1,"yes","no")
# mp$dm2=ifelse(mp$Dm_compute==1,"yes","no")
# 
# 
# c2=with(mp,
#         predict_10yr_ascvd_risk(age_years = agecompute,
#                            race = Race,sex = gender,smoke_current = smoke2,
#                            chol_total_mgdl = TC,bp_sys_mmhg = SBP,
#                            bp_meds = bpmed1,
#                            chol_hdl_mgdl = HDL,diabetes = dm2,equation_version = "Goff_2013",
#                            override_boundary_errors = TRUE,
#                            race_levels = list(black = 'AA', white = c('white', 'other')),
#                             smoke_current_levels = list(no = c('no', 'former', 'never'), yes = 'yes'),
#                         bp_meds_levels = list(no = 'no', yes = 'yes'),
#                         diabetes_levels = list(no = 'no', yes = 'yes')))
# 



###It's the same                           
                          
# d=merge(c[,c("eid","minfollow","agecompute","as2")],
#         dfh[,c("identifier","f.21003.0.0","Cad_0_Any","Cad_0_censor_age")],
#         by.x = "eid",by.y="identifier")
# d$followlevel=cut(as.numeric(c$minfollow),breaks = c(0,10,20,100))
c$followlevel=cut(as.numeric(c$minfollow),breaks = c(0,10,20,100))


a=ggplot(c) + 
aes(x = TC,fill=followlevel) +
geom_density() +
geom_density( 
    
    position="identity",
    alpha = 0.4,
   linetype = "dashed")+labs(fill="Minimum Followup",y="Density")+theme_classic()


b=ggplot(c) + 
  aes(x = HDL,fill=followlevel) +
  geom_density() +
  geom_density( 
    
    position="identity",
    alpha = 0.4,
    linetype = "dashed")+labs(fill="Minimum Followup",y="Density")+theme_classic()


c=ggplot(c) + 
  aes(x = SBP,fill=followlevel) +
  geom_density() +
  geom_density( 
    
    position="identity",
    alpha = 0.4,
    linetype = "dashed")+labs(fill="Minimum Followup",y="Density")+theme_classic()

library(ggpubr)
skew=ggarrange(a,b,c,nrow = 1,common.legend = T)
saveRDS(skew,"skew.rds")
ggplot(c[,c("statin","statin_compute")],aes(x="identity",y=count)+geom_hist(stat="identity")

d=merge(d,ascvddf,by= "eid")
qqplot(c$as2,ascvddf$as2,xlab="First Longitudinal Measurement 10-yr PCE",ylab="All 10 yr PCE At enrollment")
qqplot(c$agecompute,ascvddf$ageenroll,xlab="First Longitudinal Measurement AGE",ylab="Age at Enrollment General Cohort")
plot(d$as2.x,d$as2.y,xlab="10 Year Score First Long Measure",ylab="10 year score at enrollment")
d$as=ifelse(d$agecompute>d$ageenroll,d$as2.y,d$as2.x)
plot(d$as,d$as2.y,xlab="10 Year Score First Long Measure or if Enroll Before",ylab="10 year score at enrollment")
## some people have higher BP etcbut in general...

#c$followlevel=cut(as.numeric(c$minfollow),breaks = c(0,10,20,100))

cnew=merge(c,ascvddf,by= "eid")
#cnew=merge(c,dfh[,c("identifier","cad.prs")],by.x="eid",by.y="identifier")

plotas=ggplot(cnew[cnew$agecompute<cnew$ageenroll,],mapping = aes(as2.x,as2.y,color=as.factor(cad.prs.lec)))+
geom_point()+
  labs(x="10 Year Score First Long Measure or if Enroll Before",y="10 year score at enrollment",color="PRS Level")+
  theme_classic()
saveRDS(plotas,file = "~/multistate/plotas.rds")

sumframe=data.frame(cnew%>%group_by(followlevel)%>%summarise(antiHtnEver=mean(antihtn),
                                                          antiHtn_compute=mean(antihtn_compute),
                                                          statinEver=mean(statin),
                                                          statinCompute=mean(statin_compute),
                                                          agecompute=mean(agecompute),
                                                          ageenroll=mean(ageenroll),
                                                          prs=mean(cad.prs))
                    )


##plot ascvd at enrollment versus all

## how do we compare to MSTATE


ml=c[which(c$minfollow>20),]
ml$fac="long"
ms=c[which(c$minfollow<20),]
ms$fac=rep("short")



r=rbind(ml,ms)

set.seed(456)
enrollments=c(41:70)
aucmat=matrix(NA,nrow=length(enrollments),ncol=4)
prcmat=matrix(NA,nrow=length(enrollments),ncol=3)
ages=40:80
fixedsmoke=readRDS("~/multistate/output/fixedsmoke.rds")
enrollments=c(41:70)

s2=stateriskfunc_smoking_smoothedcoef(ages = c(40:80),prs_quants = prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
p=multipleprsfunc(s = s2[,,,1],prsprobs = pnorm(prs_quants))
m=matriskfunc_coef(p,ages,quantiles = prs_quants)

c3=merge(c,dfh[,c("identifier","cad.prs","Ht_0_censor_age","Ht_0_Any","HyperLip_0_Any","HyperLip_0_censor_age")],by.x = "eid",by.y="identifier")


c3$cad.prs=as.numeric(scale(c3$cad.prs))
c3$cad.prs.lec=cut(c3$cad.prs,breaks = c(-5.02,-0.84,0.84,5.02),labels = c("low","mid","high"))
c3$int=interaction(c3$f.31.0.0,c3$cad.prs.lec)
levels(c3$int) <- c(1,2,3,4,5,6)

c3%>%group_by(int)%>%summarise(mean(cad.prs),median(cad.prs))
## grab means
saveRDS(c3,file="~/Dropbox/c3longitudinal.rds")

###

```
c3=readRDS("~/Dropbox/c3longitudinal.rds")
c3$cad.prs=as.numeric(scale(c3$cad.prs))
c3$cad.prs.lec=cut(c3$cad.prs,breaks = c(-5.02,-0.84,0.84,5.02),labels = c("low","mid","high"))
c3$int=interaction(c3$f.31.0.0,c3$cad.prs.lec)
levels(c3$int) <- c(1,2,3,4,5,6)

prs_quants=c(data.frame(c3%>%group_by(int)%>%summarize(median(cad.prs),mean(cad.prs)))[c(2,4,6),3])
prsprobs= pnorm(prs_quants)
m=matriskfun(p,ages,quantiles = prs_quants)
c3=c3[c3$minfollow>20,]
agesint=seq(40,60,by=2)
ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
c3=data.table(c3)
for(i in 1:length(agesint)){
  age=agesint[i]
  ten.year[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
  lifetime[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
}


emp.ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
emp.lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length((agesint))){
  age=agesint[i]
  for(j in 1:length(levels(c3$int))){
    
    cat=levels(c3$int)[j]
    atrisk = c3[age < Cad_0_censor_age &
                        age < Ht_0_censor_age &
                        age < HyperLip_0_censor_age &
                        age < Dm_0_censor_age &smoke==0 , ]
    
    emp.ten.year[i,j]=compute_empiricalrisk(age=age,age2 = age+10,df_frame = atrisk,cat = cat)
    emp.lifetime[i,j]=compute_empiricalrisk(age=age,age2 = 100,df_frame = atrisk,cat = cat)
  }}



ascvd.ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

c3$ascvd_10y_accaha=c3$as2
c3$phenos.enrollment=c3$agecompute

for(i in 1:length((agesint))){
  age=agesint[i]
  for(j in 1:length(levels(c3$int))){
    cat=levels(c3$int)[j]
    atrisk=c3[age < Cad_0_censor_age &
                      age < Ht_0_censor_age &
                      age < HyperLip_0_censor_age &
                      age < Dm_0_censor_age &smoke==0 , ]
    ascvd.ten.year[i,j]=compute_pce_predictedrisk(age=age,df_frame =atrisk,cat = cat)
    #ascvdriskmat[i,2]=compute_empiricalrisk(age=40,df_frame = mpce,cat = cat)
  }
}





library(ehaGoF)

pcten=gofRMSE(Obs = as.vector(emp.ten.year*100),Prd = as.vector(ascvd.ten.year))
msten=gofRMSE(Obs = as.vector(emp.ten.year*100),Prd = as.vector(ten.year*100))
mslif=gofRMSE(Obs = as.vector(emp.lifetime*100),Prd = as.vector(lifetime*100))
pclif=gofRMSE(Obs = as.vector(emp.lifetime*100),Prd = as.vector(ascvd.ten.year))

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
diff.ascvd$score=rep("PCE",length(agesint))
diff.ascvd$age=agesint

diff.mstate$se=sd(as.matrix(sqrt(diff.mstate^2)))
diff.mstate$score=rep("MSGene",length(agesint))
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
#interaction_colors=c(brewer.pal(n = 6, name = "RdBu"))
interaction_colors <- c(brewer.pal(n = 3, name = "Reds")[1:3], brewer.pal(n = 3, name = "Blues"))

r2_male=ggplot(data = m,
               aes(x=age,
                   y= value,
                   ymin=value-se,
                   ymax=value+se,
                   fill=interaction)) +scale_fill_manual(values=interaction_colors)+
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar( position = position_dodge(), colour="black") +labs(y="RMSE 10 year risk",x="Age",fill="Genomic Level: Score")+theme_classic(base_size = 20)
#geom_point(position=position_dodge(.9), aes(y=value, colour=interaction))


ggsave(r2_male,file = "~/multistate/rmseten.png",dpi = 300,height = 4,width = 8)


# Create life RMSE
diff.mstate.life=abs(data.frame(lifetime-emp.lifetime))
d=as.matrix(diff.mstate.life)
sqrt(mean(d^2))

diff.ascvd.life=abs(data.frame(ascvd.ten.year/100-emp.lifetime))
d=as.matrix(diff.ascvd.life)
sqrt(mean(d^2))

diff.ascvd.life$se=sd(as.matrix(sqrt(diff.ascvd.life^2)))
diff.ascvd.life$score=rep("PCE",length(agesint))
diff.ascvd.life$age=agesint

diff.mstate.life$se=sd(as.matrix(sqrt(diff.mstate.life^2)))
diff.mstate.life$score=rep("MSGene",length(agesint))
diff.mstate.life$age=agesint
r=rbind(diff.ascvd.life,diff.mstate.life)
rownames(r)=NULL
rf=r[,c(1,3,5,7,8,9)]
rf$sex=rep("female",nrow(rf))

rm=r[,c(2,4,6,7,8,9)]
rm$sex=rep("male",nrow(rm))
names(rf)[1:3]=names(rm)[1:3]=c("low","medium","high")



r=rbind(diff.ascvd.life,diff.mstate.life)
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
  geom_errorbar( position = position_dodge(), colour="black") +labs(y="RMSE Lifetme risk",x="Age",fill="Genomic Level: Score")+theme_classic(base_size = 20)
#geom_point(position=position_dodge(.9), aes(y=value, colour=interaction))


ggsave(r2_male,file = "~/multistate/rmselife_fullfollow.png",dpi = 300,height = 4,width = 8)

colnames(rf)=c("Low","Intermediate","High","se","score","age")
m=melt(rf,id.vars=c("age","score","se"))

m$se=m$se/sqrt(1000)
m$interaction=interaction(m$variable,m$score)
m$value=abs(m$value)
#interaction_colors=c(brewer.pal(n = 6, name = "RdBu"))
interaction_colors <- c(brewer.pal(n = 3, name = "Reds")[1:3], brewer.pal(n = 3, name = "Blues"))

r2_female=ggplot(data = m,
               aes(x=age,
                   y= value,
                   ymin=value-se,
                   ymax=value+se,
                   fill=interaction)) +scale_fill_manual(values=interaction_colors)+
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar( position = position_dodge(), colour="black") +labs(y="RMSE Lifetme risk",x="Age",fill="Genomic Level: Score")+theme_classic(base_size = 20)
#geom_point(position=position_dodge(.9), aes(y=value, colour=interaction))


grmse=ggarrange(r2_male,r2_female,nrow = 1,common.legend = T,legend = "right",labels = c("Male","Female"))
saveRDS(grmse,"mstate/grmselife.rds")


##3 and now lifetime AUG
#c3=data.table(readRDS("~/Dropbox/c3longitudinal.rds"))
#c3=merge(c3,ascvddf[,c("eid","ageenroll","as2")],by="eid")
#c3$as2=ifelse(c3$ageenroll<c3$agecompute,c3$as2.y,c3$as2.x)
c3=c3old
c3=data.table(readRDS("~/Dropbox/c3longitudinal.rds"))
c3$cad.prs=scale(c3$cad.prs)
c3$ascvd_10y_accaha=c3$as2
set.seed(456)
enrollments=c(40:70)
aucmat=matrix(NA,nrow=length(enrollments),ncol=4)
prcmat=matrix(NA,nrow=length(enrollments),ncol=3)
ages=20:80
fixedsmoke=readRDS("~/multistate/output/fixedsmoke.rds")
#

## return a matrix of coefficients over all ages for a given state to state transition
mat=return_smoothedmatrix(start = "Health",stop = "Cad",ages = ages,modelfit = fixedsmoke)
c3=c3[,-c("age")]
#c3=c3[c3$minfollow>20,]


for(z in 1:length(enrollments)){
  age=enrollments[z]
  start=age
  stop=80
  
  df_frame=data.table(c3)
  atrisk = df_frame[age < Cad_0_censor_age &
                      age < Ht_0_censor_age &
                      age < HyperLip_0_censor_age &
                      age < Dm_0_censor_age &smoke==0 , ]
  print(dim(atrisk))
  
  df_updated=atrisk
  
  df_updated$ms=compute_prediction_product_matrix(coefmat = mat,atrisk = atrisk,agepredinterval = c(start:stop))
  ### return matrix of smoothed coefficeints
  #library(purrr)
  
  rm(atrisk)
  
  #require(pROC)
  df_updated$outcome=ifelse(df_updated$Cad_0_Any==2&df_updated$Cad_0_censor_age<stop,1,0)
  d=df_updated[round(agecompute,0)==age&smoke==0,]
  #d=d[!is.na(d$ascvd_10y_accaha),]
  aucmat[z,1]=roc(d$outcome~d$ms)$auc
  aucmat[z,2]=roc(d$outcome~d$as2)$auc
  aucmat[z,4]=roc(d$outcome~d$cad.prs)$auc
  d=d[!is.na(d$as2),]
  print(dim(d))
  h=d$as2+d$ms
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
aucplot <- ggplot(m,aes(x = Age,y = AUC,color = Model,ymin=AUC,ymax=AUC))+geom_point()+geom_line(aes(group=Model,color =Model),linewidth=3)+geom_pointrange()+theme_classic()+ggtitle(paste0("Lifetime Risk Prediction Using First Score"))


improv=mean(prcmat[,1]-prcmat[,2])*100
rownames(prcmat)=enrollments
m=melt(prcmat,id.vars="Age")
names(m)=c("Age","Model","PRauc")
m$Model=as.factor(m$Model)
levels(m$Model)[1]="MSGene"
levels(m$Model)[2]="PCE"
levels(m$Model)[3]="Combined"
levels(m$Model)[4]="PRS"
prplot <- ggplot(m,aes(x = Age,y = PRauc,color = Model,ymin=PRauc,ymax=PRauc))+geom_point()+geom_line(aes(group=Model,color =Model),linewidth=3)+geom_pointrange()+ylim(0,0.5)+theme_classic()+ggtitle(paste0("Lifetime risk prediction using First Score"))

ga=ggarrange(aucplot,prplot,common.legend = T)
saveRDS(ga,file="lifetimeriskusingfirstscore.rds")
ggsave(ga,file="lifetimeriskusingfirstscore.png")

mat=return_smoothedmatrix(start = "Health",stop = "Cad",ages = 20:80,modelfit = fixedsmoke)

cnew$f.31.0.0=as.numeric(as.character(cnew))
cnew$smoke=cnew$smoke
s=sapply(seq(1:nrow(cnew)),function(x){compute_prediction_product_matrix(cnew[x,],agepredinterval = c(cnew[x,round(agecompute,0)]:(cnew[x,round(agecompute,0)+10])),coefmat = mat)})

cnew$mstate10=s

sl=sapply(seq(1:nrow(cnew)),function(x){compute_prediction_product_matrix(cnew[x,],agepredinterval = c(cnew[x,round(agecompute,0)]:80),coefmat = mat)})

cnew$mstatelife=sl

Create the histogram for data frame 1
p1 <- ggplot(r, aes(x = as2)) +
  geom_histogram(bins = 100, fill = "blue", alpha = 0.7,aes(y=..density..)) +
  geom_density(aes(y=..density..))+xlim(c(0,20))+
  labs(title = "Distribution of Scores - Data Frame 1",
       x = "Score",
       y = "Frequency") +
  theme_minimal()+facet_wrap(vars(fac),nrow = 2)

# Create the histogram for data frame 2
p2 <- ggplot(ms, aes(x = as2)) +
  geom_histogram(binwidth = 5, fill = "green", alpha = 0.7) +
  geom_density()+
  labs(title = "Distribution of Scores - Data Frame 2",
       x = "Score",
       y = "Frequency") +
  theme_minimal()

# Combine and display the histograms side by side using facet_wrap
plot_combined <- p1 + p2 + facet_wrap(~ ., ncol = 2)

# Display the combined plot
print(plot_combined)


compute_CVrisk2=function (df, scores = c("ascvd_10y_accaha", "as2", 
                                         "ascvd_10y_frs_simple", "chd_10y_mesa", "chd_10y_mesa_cac"), 
                          age, gender, race, sbp = NULL, bmi = NULL, hdl = NULL, totchol = NULL, 
                          bp_med = NULL, smoker = NULL, diabetes = NULL, lipid_med = NULL, 
                          fh_heartattack = NULL, cac = NULL) 
{
  all_args <- as.list(environment())
  valid_pred <- c("age", "gender", "race", "sbp", "bmi", "hdl", 
                  "totchol", "bp_med", "smoker", "diabetes", "lipid_med", 
                  "fh_heartattack", "cac")
  pred_args <- list()
  for (var in valid_pred) {
    if (!is.null(all_args[[var]])) 
      pred_args[[var]] <- df[[all_args[[var]]]]
  }
  results <- sapply(scores, function(x) do.call(x, pred_args))
  row.names(results) <- NULL
  return(cbind(df, results))
}

as2=function (race = "white", gender = c("male", "female"), age, 
              totchol, hdl, sbp, bp_med, smoker, diabetes, ...) 
{
  if (!all(gender %in% c("male", "female")) | missing(gender)) {
    stop("gender must be either 'male' or 'female'")
  }
  ascvd_pooled_coef <- NULL
  utils::data(ascvd_pooled_coef, envir = environment())
  race <- ifelse(race %in% c("white", "aa"), race, "white")
  race_sex <- data.frame(race, gender)
  race_sex$id <- as.numeric(row.names(race_sex))
  pooled_coef <- merge(race_sex, ascvd_pooled_coef)
  pooled_coef <- pooled_coef[order(pooled_coef$id), ]
  sbp_treated <- ifelse(bp_med == 1, sbp, 1)
  sbp_untreated <- ifelse(bp_med == 0, sbp, 1)
  indv_sum <- log(age) * pooled_coef$ln_age + log(age)^2 * 
    pooled_coef$ln_age_squared + log(totchol) * pooled_coef$ln_totchol + 
    log(age) * log(totchol) * pooled_coef$ln_age_totchol + 
    log(hdl) * pooled_coef$ln_hdl + log(age) * log(hdl) * 
    pooled_coef$ln_age_hdl + log(sbp_treated) * pooled_coef$ln_treated_sbp + 
    log(sbp_treated) * log(age) * pooled_coef$ln_age_treated_sbp + 
    log(sbp_untreated) * pooled_coef$ln_untreated_sbp + log(sbp_untreated) * 
    log(age) * pooled_coef$ln_age_untreated_sbp + smoker * 
    pooled_coef$smoker + smoker * log(age) * pooled_coef$ln_age_smoker + 
    diabetes * pooled_coef$diabetes
  risk_score <- round((1 - (pooled_coef$baseline_survival^exp(indv_sum - 
                                                                pooled_coef$group_mean))) * 100, 2)
  ifelse(risk_score < 1, 1, risk_score)
}


NRI_mat=matrix(NA,nrow=2,ncol=6)
num_mat=matrix(NA,nrow=2,ncol=6)
ages=c(40,45,50,55,60,65,70)
cnew=readRDS("~/Dropbox/c3longitudinal_withms.rds")
oldf=cnew
#for(q in 1:(length(ages)-1)){
for(q in 1:(length(ages)-1)){
  is=ages[q]
  i=ages[q+1]
  atrisk = oldf[i < Cad_0_censor_age &
  i < Ht_0_censor_age &
  i < HyperLip_0_censor_age &
  i < Dm_0_censor_age &smoke==0 , ]

  #print(i)
  df=oldf[oldf$agecompute>is&oldf$agecompute<i,]
  #df=oldf[oldf$agecompute<is,]
  print(nrow(df))
  (nocad_PCEpos_PRSneg=dim(df[df$Cad_0_Any==1&df$as2>7.5&df$msl<0.10,])[1])
  (nocad_PCEneg_PRSpos=dim(df[df$Cad_0_Any==1&df$as2<7.5&df$msl>0.10,])[1])
  
  
  (cad_PCEpos_PRSneg=dim(df[df$Cad_0_Any==2&df$as2>7.5&df$msl<0.10,])[1])
  (cad_PCEneg_PRSpos=dim(df[df$Cad_0_Any==2&df$as2<7.5&df$msl>0.10,])[1])
  
  
  NRI_e=(cad_PCEneg_PRSpos-cad_PCEpos_PRSneg)/(sum(df$Cad_0_Any==2))
  NRI_ne=(nocad_PCEpos_PRSneg-nocad_PCEneg_PRSpos)/(sum(df$Cad_0_Any==1))
  
  #print(NRI_e)
  #print(NRI_ne)
  
  NRI_mat[1,q]=NRI_e
  NRI_mat[2,q]=NRI_ne
  num_mat[c(1,2),q]=nrow(df)
}

errors=sqrt((abs(NRI_mat)*(1-abs(NRI_mat)))/(num_mat))
rownames(NRI_mat)=rownames(errors)=c("NRI_event","NRI_nonevent")
t=data.frame(t(NRI_mat))
e=data.frame(t(errors))
t$ages=e$ages=c("40-45","45-50","50-55","55-60","60-65","65-70")
m=melt(t,id.vars = "ages")
m2=melt(e,id.vars = "ages")
m$error=m2$value
m$variable=factor(m$variable,levels = c("NRI_event","NRI_nonevent"),labels = c("Net Reclassification Event","Net Reclassification Non-event"))

n=ggplot(m,aes(ages,value,fill=variable,ymin=value-error,ymax=value+error))+geom_bar(stat="identity",position = "dodge")+theme_classic(base_size = 20)+
  labs(x="Ages of Enrollment, years",y="Net Reclassification Index",fill="Metric Considered")+
  scale_fill_nejm()+ylim(c(-1,1))+
  geom_errorbar(aes(ymin=value-0.02,ymax=value+0.02),width=.2,position=position_dodge(.9))
saveRDS(n,"~/multistate/nri.rds")
