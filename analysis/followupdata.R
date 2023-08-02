library(data.table)
library(lubridate)
library('dplyr')
source("~/dynamichr/code/utils.R")
###

fos=read.table("~/Desktop/sraempty/phs000007.v33.pht006027.v4.p14.c1.vr_wkthru_ex09_1_1001s.HMB-IRB-MDS.txt.gz",header = T,skip = 1,sep="\t")
fhs=read.table("~/Desktop/sraempty/phs000007.v33.pht007777.v3.p14.c1.vr_wkthru_ex32_0_0997s.HMB-IRB-MDS.txt.gz",header = T,skip = 1,sep="\t")
#surv=read.table("/Users/sarahurbut/Desktop/Framingham/85598/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v32.p13.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v32.pht003316.v9.p13.c1.vr_survcvd_2018_a_1267s.HMB-IRB-MDS.txt.gz",skip=1,sep="\t",header=T)
surv=read.table("/Users/sarahurbut/Desktop/sraempty/phs000007.v33.pht003316.v10.p14.c1.vr_survcvd_2019_a_1334s.HMB-IRB-MDS.txt.gz",skip=1,sep="\t",header=T)

m=merge(fhs,surv,by = "shareid")

SEX=grep(x = names(m2),pattern="SEX")

HDL=m[,c(1,grep(x = names(m),pattern="HDL"))]
melt(HDL,id.vars = "shareid")

TC=m[,c(1,grep(x = names(m),pattern="TC"))]
Trig=m[,c(1,grep(x = names(m),pattern="TRIG"))]
SMOKE=m[,c(1,grep(x = names(m),pattern="CURR"))]
LIP=m[,c(1,grep(x = names(m),pattern="LIPRX"))]
HRX=m[,c(1,grep(x = names(m),pattern="HRX"))]
LDL=m[,c(1,grep(x = names(m),pattern="LDL"))]



# 
# for(i in 1:32){
#   df=m[,"dbGaP_Subject_ID.x","SEX",paste0("AGE",i),"CURRSMK1","DMRX1","TC1","HDL]
#   compute_CVrisk2(df = m,scores = "as2",age = paste0("AGE",i),gender = )
# }

m2=merge(fos,surv,by="shareid")
HDL=m[,c(1,grep(x = names(m),pattern="HDL"))]
sum(rowSums(!is.na(HDL[,-1]))==0)

#######
HDL=m2[,c(1,grep(x = names(m2),pattern="HDL"))]
sum(rowSums(!is.na(HDL[,-1]))==0)
HDL=HDL%>%
  pivot_longer(
    cols = starts_with("HDL"),
    names_to = "HDL",
    names_prefix = "HDL",
    values_to = "Value",
    values_drop_na = FALSE
  )

#######
HRX=m2[,c(1,grep(x = names(m2),pattern="HRX"))]
sum(rowSums(!is.na(HRX[,-1]))==0)
HRX=HRX%>%
  pivot_longer(
    cols = starts_with("HRX"),
    names_to = "HRX",
    names_prefix = "HRX",
    values_to = "Value",
    values_drop_na = FALSE
  )

#######
LIPRX=m2[,c(1,grep(x = names(m2),pattern="LIPRX"))]
sum(rowSums(!is.na(LIPRX[,-1]))==0)
LIPRX=LIPRX%>%
  pivot_longer(
    cols = starts_with("LIPRX"),
    names_to = "LIPRX",
    names_prefix = "LIPRX",
    values_to = "Value",
    values_drop_na = FALSE
  )


###
AGE=m2[,c(1,grep(x = names(m2),pattern="AGE"))]
sum(rowSums(!is.na(AGE[,-1]))==0)
AGE=AGE%>%
  pivot_longer(
    cols = starts_with("AGE"),
    names_to = "AGE",
    names_prefix = "AGE",
    values_to = "Value",
    values_drop_na = FALSE
  )

###
SMOKE=m2[,c(1,grep(x = names(m2),pattern="SMK"))]
sum(rowSums(!is.na(SMOKE[,-1]))==0)

SMOKE=SMOKE%>%
  pivot_longer(
    cols = starts_with("CURR"),
    names_to = "SMOKE",
    names_prefix = "SMOKE",
    values_to = "Value",
    values_drop_na = FALSE
  )

###
DM=m2[,c(1,grep(x = names(m2),pattern="DMRX"))]
sum(rowSums(!is.na(DM[,-1]))==0)
DM=DM%>%
  pivot_longer(
    cols = starts_with("DM"),
    names_to = "DM",
    names_prefix = "DM",
    values_to = "Value",
    values_drop_na = FALSE
  )

###
TC=m2[,c(1,grep(x = names(m2),pattern="TC"))]
sum(rowSums(!is.na(TC[,-1]))==0)
TC=TC%>%
  pivot_longer(
    cols = starts_with("TC"),
    names_to = "TC",
    names_prefix = "TC",
    values_to = "Value",
    values_drop_na = FALSE
  )

###
SBP=m2[,c(1,grep(x = names(m2),pattern="SBP"))]
sum(rowSums(!is.na(SBP[,-1]))==0)
SBP=SBP%>%
  pivot_longer(
    cols = starts_with("SBP"),
    names_to = "SBP",
    names_prefix = "SBP",
    values_to = "Value",
    values_drop_na = FALSE
  )



# Takes too long
# f=merge(SMOKE,
#               merge(AGE,
#                         merge(SBP,
#                                   merge(HDL,DM,by = "shareid",suffixes = c("HDL","DM")),
#                                         by="shareid"),
#                               by = "shareid",suffixes = c("AGE","SBP")),
#                   by="shareid")

f=cbind(AGE,SBP[,-1],HDL[,-1],DM[,-1],SMOKE[,-1],HRX[,-1],TC[,-1],LIPRX[,-1])
colnames(f)[c(3,5,7,9,11,13,15,17)]=c("AGEV","SBPV","HDLV","DMV","SMOKEV","HRXV","TCV","LIPRXV")



f2=f%>%fill(c("AGEV","SBPV","HDLV","DMV","SMOKEV","HRXV","TCV","LIPRXV"),.direction = "downup")
f2$SEX=rep((fos$SEX),each=9)
f2$SEX=ifelse(f2$SEX==2,"female","male")

sbu=SBP%>%group_by(shareid)%>%summarise(length(unique(Value)))
u=sample(sbu$shareid[sbu$`length(unique(Value))`==9],9)
s2=f2[f2$shareid%in%u,]
ggplot(s2, aes(AGEV, SBPV, col = factor(shareid))) +
geom_point() +
geom_line() +
facet_wrap(~ shareid) +
labs(x = "Age,y", y = "SBP,mmHg", col = "Pt id") +
guides(col = guide_legend(nrow = 3))


cframe=compute_CVrisk2(f2,scores = "as2",age = "AGEV",gender = "SEX",
                       race = "white",sbp = "SBPV",hdl = "HDLV",
                bp_med ="HRXV",totchol = "TCV",diabetes = "DMV",smoker = "SMOKEV",lipid_med = "LIPRXV")

r=reshape(cframe[,c("shareid","AGE","as2")],idvar = "shareid",v.names = "as2",timevar = "AGE",direction = "wide")
sum(m2$chd==1)
sum(m$chd==1)

coolfhsplot=ggplot(cframe[cframe$AGEV>20,],aes(AGEV,as2,col=as.factor(AGEV)))+geom_point()+labs(y="10-year Risk",x="AGE")
saveRDS(coolfhsplot,file="~/Dropbox/sequentialPCE.rds")

sbu=SBP%>%group_by(shareid)%>%summarise(length(unique(Value)))
u=sample(sbu$shareid[sbu$`length(unique(Value))`==9],9)
s2=cframe[cframe$shareid%in%u,]
p=ggplot(s2, aes(AGEV, as2, col = factor(shareid))) +
  geom_point() +
  geom_line() +
  facet_wrap(~ shareid) +
  labs(x = "Age,y", y = "PCE", col = "Pt id") +
  guides(col = guide_legend(nrow = 3))+theme_classic()

s=ggplot(s2, aes(AGEV, SBPV, col = factor(shareid))) +
  geom_point() +
  geom_line() +
  facet_wrap(~ shareid) +
  labs(x = "Age,y", y = "SBPV", col = "Pt id") +
  guides(col = guide_legend(nrow = 3))+theme_classic()

t=ggplot(s2, aes(AGEV, TCV, col = factor(shareid))) +
  geom_point() +
  geom_line() +
  facet_wrap(~ shareid) +
  labs(x = "Age,y", y = "TC", col = "Pt id") +
  guides(col = guide_legend(nrow = 3))+theme_classic()

> gf=ggarrange(p,s,t,common.legend = T,nrow=1)
> saveRDS(gf,file = "gf.rds")
### compare with ours

## plots
##How to use
prs=fread("~/Dropbox/Fram_allchr_CAD_c1.profile")
prs$score=scale(prs$V2)
m3=merge(cframe,prs,by.x="shareid",by.y="V1")


## for more levels
mpce=data.frame(readRDS("~/multistate/output/mpcecomplete.rds"))
prs_quants=qnorm(c(seq(0.1,0.9,by=0.1)))
mpce$newlev=cut(mpce$cad.prs,breaks = c(-5,prs_quants,5),labels = c(1:10))
m3$newlev=cut(m3$score,breaks = c(-5,prs_quants,5),labels = c(1:10))
mpce$newint=interaction(mpce$f.31.0.0,mpce$newlev)
levels(mpce$newint)=c(1:20)
m3$newint=interaction(ifelse(m3$SEX=="male",1,0),m3$newlev)
levels(m3$newint)=c(1:20)

mpce%>%group_by(newint)%>%summarise(mean(cad.prs),median(cad.prs))
## grab means

prs_quants=c(data.frame(mpce%>%group_by(newint)%>%summarize(median(cad.prs),mean(cad.prs)))[c(2,4,6,8,10,12,14,16,18,20),3])
prsprobs= pnorm(prs_quants)


## generate new arrays
source("~/multistate/code/fitarray.R")

fixedsmoke=readRDS("~/multistate/output/fixedsmoke.rds")
s=stateriskfunc_smoking(ages = c(20:80),prs_quants = prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)


p=multipleprsfunc(s = s[,,,1],prsprobs = pnorm(prs_quants))
ages=c(20:80)
m=matriskfun(p,ages,quantiles = prs_quants)
agesint=seq(20,70,by=1)
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
colnames(ten.year)=c(1:20)
ten.year$age=agesint

rownames(lifetime)=agesint
colnames(lifetime)=c(1:20)
lifetime$age=agesint



lookup_table <- data.frame(melt(ten.year,id.vars = c("age")))
names(lookup_table)[3]="tenyear"
ggplot(lookup_table,aes(age,y = tenyear,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")

lookup_table$variable=as.factor(lookup_table$variable)
# Join the lookup table with the original data frame based on the factor variable
df_updated <- m3 %>% 
  left_join(lookup_table, by = c("newint" = "variable","AGEV"="age"))

lookup_table2 <- data.frame(melt(lifetime,id.vars = c("age")))
names(lookup_table2)[3]="lifetime"

ggplot(lookup_table2,aes(age,y = lifetime,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")

df_updated <- df_updated %>% 
  left_join(lookup_table2, by = c("newint" = "variable","AGEV"="age"))



####
### for Hypertensive

s=stateriskfunc_smoking(ages = c(20:80),prs_quants = prs_quants,start = "Ht",
                        stop = "Cad",modelfit = fixedsmoke)


p=multipleprsfunc(s = s[,,,1],prsprobs = pnorm(prs_quants))
ages=c(20:80)
m=matriskfun(p,ages,quantiles = prs_quants)
agesint=seq(20,70,by=1)
ten.yearh=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetimeh=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.yearh[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
  lifetimeh[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
}

ten.yearh=data.frame(ten.yearh)
lifetimeh=data.frame(lifetimeh)

rownames(ten.yearh)=agesint
colnames(ten.yearh)=c(1:20)
ten.yearh$age=agesint

rownames(lifetimeh)=agesint
colnames(lifetimeh)=c(1:20)
lifetimeh$age=agesint



lookup_table <- data.frame(melt(ten.yearh,id.vars = c("age")))
names(lookup_table)[3]="tenyearh"
ggplot(lookup_table,aes(age,y = tenyearh,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")

lookup_table$variable=as.factor(lookup_table$variable)
# Join the lookup table with the original data frame based on the factor variable
df_updated <- df_updated %>% 
  left_join(lookup_table, by = c("newint" = "variable","AGEV"="age"))

lookup_table2 <- data.frame(melt(lifetimeh,id.vars = c("age")))
names(lookup_table2)[3]="lifetimeh"

ggplot(lookup_table2,aes(age,y = lifetimeh,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")

df_updated <- df_updated %>% 
  left_join(lookup_table2, by = c("newint" = "variable","AGEV"="age"))

#################
## Hyperlip

####
s=stateriskfunc_smoking(ages = c(20:80),prs_quants = prs_quants,start = "HyperLip",
                        stop = "Cad",modelfit = fixedsmoke)


p=multipleprsfunc(s = s[,,,1],prsprobs = pnorm(prs_quants))
ages=c(20:80)
m=matriskfun(p,ages,quantiles = prs_quants)
agesint=seq(20,70,by=1)
ten.yearhl=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetimehl=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.yearhl[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
  lifetimehl[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
}

ten.yearhl=data.frame(ten.yearhl)
lifetimehl=data.frame(lifetimehl)

rownames(ten.yearhl)=agesint
colnames(ten.yearhl)=c(1:20)
ten.yearhl$age=agesint

rownames(lifetimehl)=agesint
colnames(lifetimehl)=c(1:20)
lifetimehl$age=agesint



lookup_table <- data.frame(melt(ten.yearhl,id.vars = c("age")))
names(lookup_table)[3]="ten.yearhl"
ggplot(lookup_table,aes(age,y = ten.yearhl,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")

lookup_table$variable=as.factor(lookup_table$variable)
# Join the lookup table with the original data frame based on the factor variable
df_updated <- df_updated %>% 
  left_join(lookup_table, by = c("newint" = "variable","AGEV"="age"))

lookup_table2 <- data.frame(melt(lifetimehl,id.vars = c("age")))
names(lookup_table2)[3]="lifetimehl"

ggplot(lookup_table2,aes(age,y = lifetimehl,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")

df_updated <- df_updated %>% 
  left_join(lookup_table2, by = c("newint" = "variable","AGEV"="age"))



################
### DM

####
s=stateriskfunc_smoking(ages = c(20:80),prs_quants = prs_quants,start = "Dm",
                        stop = "Cad",modelfit = fixedsmoke)


p=multipleprsfunc(s = s[,,,1],prsprobs = pnorm(prs_quants))
ages=c(20:80)
m=matriskfun(p,ages,quantiles = prs_quants)
agesint=seq(20,70,by=1)
ten.yeardm=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetimedm=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.yeardm[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
  lifetimedm[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
}

ten.yeardm=data.frame(ten.yeardm)
lifetimedm=data.frame(lifetimedm)

rownames(ten.yeardm)=agesint
colnames(ten.yeardm)=c(1:20)
ten.yeardm$age=agesint

rownames(lifetimedm)=agesint
colnames(lifetimedm)=c(1:20)
lifetimedm$age=agesint



lookup_table <- data.frame(melt(ten.yeardm,id.vars = c("age")))
names(lookup_table)[3]="ten.yeardm"
ggplot(lookup_table,aes(age,y = ten.yeardm,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")

lookup_table$variable=as.factor(lookup_table$variable)
# Join the lookup table with the original data frame based on the factor variable
df_updated <- df_updated %>% 
  left_join(lookup_table, by = c("newint" = "variable","AGEV"="age"))

lookup_table2 <- data.frame(melt(lifetimedm,id.vars = c("age")))
names(lookup_table2)[3]="lifetimedm"

ggplot(lookup_table2,aes(age,y = lifetimedm,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")

df_updated <- df_updated %>% 
  left_join(lookup_table2, by = c("newint" = "variable","AGEV"="age"))



###
saveRDS(df_updated,"~/Dropbox/framinghamlookup.rds")
dfu=readRDS("~/Dropbox/framinghamlookup.rds")

## only do predictions for those not on statins and non smokers
dfu=df_updated[df_updated$LIPRXV==0&df_updated$SMOKEV==0,]



###
dfu$hypertension=ifelse(dfu$SBPV > 150,1,0)
dfu$hyperlip=ifelse(dfu$TCV > 220,1,0)
dfu$dm=ifelse(dfu$DMV==1,1,0)

one_cond=dfu[rowSums(dfu[,c("hypertension","hyperlip","dm")])<2,]
one_cond=data.table(one_cond)

one_cond[hyperlip==1,"risk"]=one_cond[hyperlip==1,ten.yearhl]
one_cond[hypertension==1,"risk"]=one_cond[hypertension==1,tenyearh]
one_cond[dm==1,"risk"]=one_cond[dm==1,ten.yeardm]
one_cond=data.frame(one_cond)
one_cond[which(rowSums(one_cond[,c("hypertension","hyperlip","dm")])==0),"risk"]=one_cond[which(rowSums(one_cond[,c("hypertension","hyperlip","dm")])==0),"tenyear"]

saveRDS(one_cond,file = "~/Dropbox/one_cond_nosmoke.rds")
###

## add smokers

fixedsmoke=readRDS("~/multistate/output/fixedsmoke.rds")
s=stateriskfunc_smoking(ages = c(20:80),prs_quants = prs_quants,
                        start = "Health",stop = "Cad",modelfit = fixedsmoke)


p=multipleprsfunc(s = s[,,,"smoke"],prsprobs = pnorm(prs_quants))
ages=c(20:80)
m=matriskfun(p,ages,quantiles = prs_quants)
agesint=seq(20,70,by=1)
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
colnames(ten.year)=c(1:20)
ten.year$age=agesint

rownames(lifetime)=agesint
colnames(lifetime)=c(1:20)
lifetime$age=agesint



lookup_table <- data.frame(melt(ten.year,id.vars = c("age")))
names(lookup_table)[3]="tenyear.s"
ggplot(lookup_table,aes(age,y = tenyear.s,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")

lookup_table$variable=as.factor(lookup_table$variable)
# Join the lookup table with the original data frame based on the factor variable
df_updated <- m3 %>% 
  left_join(lookup_table, by = c("newint" = "variable","AGEV"="age"))

lookup_table2 <- data.frame(melt(lifetime,id.vars = c("age")))
names(lookup_table2)[3]="lifetime.s"

ggplot(lookup_table2,aes(age,y = lifetime.s,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")

df_updated <- df_updated %>% 
  left_join(lookup_table2, by = c("newint" = "variable","AGEV"="age"))



####
### for Hypertensive

s=stateriskfunc_smoking(ages = c(20:80),prs_quants = prs_quants,start = "Ht",
                        stop = "Cad",modelfit = fixedsmoke)


p=multipleprsfunc(s = s[,,,2],prsprobs = pnorm(prs_quants))
ages=c(20:80)
m=matriskfun(p,ages,quantiles = prs_quants)
agesint=seq(20,70,by=1)
ten.yearh=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetimeh=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.yearh[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
  lifetimeh[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
}

ten.yearh=data.frame(ten.yearh)
lifetimeh=data.frame(lifetimeh)

rownames(ten.yearh)=agesint
colnames(ten.yearh)=c(1:20)
ten.yearh$age=agesint

rownames(lifetimeh)=agesint
colnames(lifetimeh)=c(1:20)
lifetimeh$age=agesint



lookup_table <- data.frame(melt(ten.yearh,id.vars = c("age")))
names(lookup_table)[3]="tenyearh.s"
ggplot(lookup_table,aes(age,y = tenyearh.s,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")

lookup_table$variable=as.factor(lookup_table$variable)
# Join the lookup table with the original data frame based on the factor variable
df_updated <- df_updated %>% 
  left_join(lookup_table, by = c("newint" = "variable","AGEV"="age"))

lookup_table2 <- data.frame(melt(lifetimeh,id.vars = c("age")))
names(lookup_table2)[3]="lifetimeh.s"

ggplot(lookup_table2,aes(age,y = lifetimeh.s,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")

df_updated <- df_updated %>% 
  left_join(lookup_table2, by = c("newint" = "variable","AGEV"="age"))

#################
## Hyperlip

####
s=stateriskfunc_smoking(ages = c(20:80),prs_quants = prs_quants,start = "HyperLip",
                        stop = "Cad",modelfit = fixedsmoke)


p=multipleprsfunc(s = s[,,,2],prsprobs = pnorm(prs_quants))
ages=c(20:80)
m=matriskfun(p,ages,quantiles = prs_quants)
agesint=seq(20,70,by=1)
ten.yearhl=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetimehl=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.yearhl[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
  lifetimehl[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
}

ten.yearhl=data.frame(ten.yearhl)
lifetimehl=data.frame(lifetimehl)

rownames(ten.yearhl)=agesint
colnames(ten.yearhl)=c(1:20)
ten.yearhl$age=agesint

rownames(lifetimehl)=agesint
colnames(lifetimehl)=c(1:20)
lifetimehl$age=agesint



lookup_table <- data.frame(melt(ten.yearhl,id.vars = c("age")))
names(lookup_table)[3]="ten.yearhls"
ggplot(lookup_table,aes(age,y = ten.yearhls,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")

lookup_table$variable=as.factor(lookup_table$variable)
# Join the lookup table with the original data frame based on the factor variable
df_updated <- df_updated %>% 
  left_join(lookup_table, by = c("newint" = "variable","AGEV"="age"))

lookup_table2 <- data.frame(melt(lifetimehl,id.vars = c("age")))
names(lookup_table2)[3]="lifetimehls"

ggplot(lookup_table2,aes(age,y = lifetimehls,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")

df_updated <- df_updated %>% 
  left_join(lookup_table2, by = c("newint" = "variable","AGEV"="age"))



################
### DM

####
s=stateriskfunc_smoking(ages = c(20:80),prs_quants = prs_quants,start = "Dm",
                        stop = "Cad",modelfit = fixedsmoke)


p=multipleprsfunc(s = s[,,,2],prsprobs = pnorm(prs_quants))
ages=c(20:80)
m=matriskfun(p,ages,quantiles = prs_quants)
agesint=seq(20,70,by=1)
ten.yeardm=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetimedm=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.yeardm[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
  lifetimedm[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
}

ten.yeardm=data.frame(ten.yeardm)
lifetimedm=data.frame(lifetimedm)

rownames(ten.yeardm)=agesint
colnames(ten.yeardm)=c(1:20)
ten.yeardm$age=agesint

rownames(lifetimedm)=agesint
colnames(lifetimedm)=c(1:20)
lifetimedm$age=agesint



lookup_table <- data.frame(melt(ten.yeardm,id.vars = c("age")))
names(lookup_table)[3]="ten.yeardms"
ggplot(lookup_table,aes(age,y = ten.yeardms,color=as.factor(variable)))+stat_smooth()+labs(x="Age",y="Ten Year Risk",col="PRS:Sex")

lookup_table$variable=as.factor(lookup_table$variable)
# Join the lookup table with the original data frame based on the factor variable
df_updated <- df_updated %>% 
  left_join(lookup_table, by = c("newint" = "variable","AGEV"="age"))

lookup_table2 <- data.frame(melt(lifetimedm,id.vars = c("age")))
names(lookup_table2)[3]="lifetimedms"

ggplot(lookup_table2,aes(age,y = lifetimedms,color=as.factor(variable)))+
  stat_smooth()+labs(x="Age",y="Lifetime Risk",col="PRS:Sex")

df_updated <- df_updated %>% 
  left_join(lookup_table2, by = c("newint" = "variable","AGEV"="age"))


saveRDS(df_updated,"~/Dropbox/framinghamlookup_smoking.rds")

## only do predictions for those not on statins and non smokers
dfu=readRDS("~/Dropbox/framinghamlookup_smoking.rds")
dfu=df_updated[df_updated$LIPRXV==0&df_updated$SMOKEV==1,]



###
dfu$hypertension=ifelse(dfu$SBPV > 150,1,0)
dfu$hyperlip=ifelse(dfu$TCV > 220,1,0)
dfu$dm=ifelse(dfu$DMV==1,1,0)

one_cond=dfu[rowSums(dfu[,c("hypertension","hyperlip","dm")])<2,]
one_cond=data.table(one_cond)
one_cond[hyperlip==1,"risk"]=one_cond[hyperlip==1,ten.yearhls]
one_cond[hypertension==1,"risk"]=one_cond[hypertension==1,tenyearh.s]
one_cond[dm==1,"risk"]=one_cond[dm==1,ten.yeardms]
one_cond=data.frame(one_cond)
one_cond[which(rowSums(one_cond[,c("hypertension","hyperlip","dm")])==0),"risk"]=one_cond[which(rowSums(one_cond[,c("hypertension","hyperlip","dm")])==0),"tenyear.s"]


saveRDS(one_cond,file = "~/Dropbox/one_cond_smoke.rds")

ocsmoke=readRDS("~/Dropbox/one_cond_smoke.rds")
oc_nosmoke=readRDS("~/Dropbox/one_cond_nosmoke.rds")

data_all <- rbind(ocsmoke,         # Rename columns & rbind
                  setNames(oc_nosmoke, names(ocsmoke)))
da=data_all[order(data_all$shareid, data_all$AGEV),]


data_all                         # Print row-binded 
###
pheno=f2%>%group_by(shareid)%>%summarise(sum(SMOKEV)>0,sum(HRXV)>0,sum(LIPRXV)>0,max(TCV),max(SBPV),sum(DMV)>0)

g=merge(df_updated,surv[,c("shareid","chd","chddate")],by="shareid")



### merge FOS with genetics
### 
c2=with(f2,
        predict_10yr_ascvd_risk(age_years = AGEV,
                                race = rep("white",34524),
                                sex = SEX,smoke_current = SMOKEV,
                                chol_total_mgdl = TCV,bp_sys_mmhg = SBPV,
                                bp_meds = HRXV,
                                chol_hdl_mgdl = HDLV,diabetes = DMV,equation_version = "Goff_2013",
                                override_boundary_errors = TRUE,
                                race_levels = list(black = 'AA', white = c('white', 'other')),
                                smoke_current_levels = list(no = '0', yes = '1'),
                                bp_meds_levels = list(no = '0', yes = '1'),
                                diabetes_levels = list(no = '0', yes = '1')))


cframe=cframe[-which(cframe$SBPV>200),]

cframe=cframe[-which(cframe$TCV>300),]
r=reshape(cframe[,c("shareid","AGE","as2")],idvar = "shareid",v.names = "as2",timevar = "AGE",direction = "wide")
ggplot(cframe,aes(AGEV,as2,col=SEX))+geom_point()
ggplot(cframe[cframe$AGEV>20,],aes(AGEV,as2,col=as.factor(AGEV)))+geom_point()+labs(y="10-year Risk",x="AGE")

tolerance=cframe
names(tolerance)[1]="id"
names(tolerance)[2]="time"

p <- ggplot(data = tolerance, aes(x = time, y = SBPV, group = id))

p + geom_line()


p + geom_line() + stat_smooth(aes(group = 1)) + stat_summary(aes(group = 1),
                                                             geom = "point", fun.y = mean, shape = 17, size = 3) + 
  facet_grid(. ~ male)
melted=melt(cframe[,c("shareid","AGE","SBPV","HDLV","TCV")],id.vars = c("shareid","AGE"))
ggplot(melted,aes(x=AGE,y=value,group=shareid))+geom_line()+facet_wrap(~variable,nrow = 3,scales ="free")


## now find individuals who we can also compute PRS on


m4=merge(m2,prs,by.x="shareid",by.y="V1")

###
sum(m4$chd==1)
505
sum(m3$chd==1)
247

sum(m4$AGE1>40&m4$chd==1)


sum(m3$AGE1>40&m3$chd==1)
### Now ask about SBP and LDL

ldl_gp=readRDS("~/Dropbox (Personal)/pheno_dir/output/ldl_gp.rds")
dfh=readRDS("~/Dropbox (Personal)/pheno_dir/output/merged_pheno_censor_final.rds")
m=merge(ldl_gp,dfh[,c("identifier","Birthdate")],by.x = "eid",by.y = "identifier")

m2=m[with(m, order(event_dt)), ]
m3=m2[!duplicated(m2$eid), ]
lc=ldl_gp%>%group_by(eid)%>%summarise(n=length(event_dt))

m3$age=as.numeric(difftime(time1 = m3$event_dt,time2 = m3$Birthdate,units = "days")/365.25)
m3=m3[m3$age>18,]

## show that LDL often followed later even though first reading earlier
m3$followup=as.numeric(difftime(time1 ="2021-03-31",time2 = m3$event_dt,units = "days")/365.25)
hist(m3$followup[m3$followup>20],xlab="Years Follow up from first measurement in UKB EHR",ylab="Frequency", main="Lipids first measurement")


#### SBPP

sbp_gp=readRDS("~/Dropbox (Personal)/pheno_dir/output/sbp_gp.rds")
dfh=readRDS("~/Dropbox (Personal)/pheno_dir/output/merged_pheno_censor_final.rds")
m=merge(sbp_gp,dfh[,c("identifier","Birthdate")],by.x = "eid",by.y = "identifier")

m2=m[with(m, order(event_dt)), ]
m3=m2[!duplicated(m2$eid), ]
lc=sbp_gp%>%group_by(eid)%>%summarise(n=length(event_dt))

m3$age=as.numeric(difftime(time1 = m3$event_dt,time2 = m3$Birthdate,units = "days")/365.25)
m3=m3[m3$age>18,]

## show that SBP often followed later even though first reading earlier
m3$followup=as.numeric(difftime(time1 ="2021-03-31",time2 = m3$event_dt,units = "days")/365.25)
hist(m3$followup[m3$followup>20],xlab="Years Follow up from first measurement in UKB EHR",ylab="Frequency", main="SBP first measurement")

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7717266/#sup1

## For RMSE, sort by first age and PCE


###
saveRDS(df_updated,"~/Dropbox/framinghamlookup.rds")
dfu=readRDS("~/Dropbox/framinghamlookup.rds")

## only do predictions for those not on statins and non smokers
dfu=dfu[dfu$LIPRXV==0&dfu$SMOKEV==0,]

firstmeasure=data.frame(dfu%>%group_by(shareid)%>%summarise(min(AGEV),as2[1],score[1],SEX[1]))
names(firstmeasure)=c("shareid","Age","ASCVD","cad.prs","sex")
firstmeasure=merge(firstmeasure,surv[,c("shareid","chd","chddate")],by="shareid")

firstmeasure$SEX=ifelse(firstmeasure$sex=="male",1,0)
firstmeasure$cad.prs.lec=cut(firstmeasure$cad.prs,breaks = c(-5.02,-0.84,0.84,5.02),labels = c("low","mid","high"))
firstmeasure$int=interaction(firstmeasure$SEX,firstmeasure$cad.prs.lec)
levels(firstmeasure$int) <- c(1,2,3,4,5,6)

firstmeasure%>%group_by(int)%>%summarise(mean(cad.prs),median(cad.prs))
## grab means
saveRDS(firstmeasure,file="~/Dropbox/firstmeasure.rds")

prs_quants=c(data.frame(firstmeasure%>%group_by(int)%>%summarize(median(cad.prs),mean(cad.prs)))[c(2,4,6),3])
prsprobs= pnorm(prs_quants)

s2=stateriskfunc_smoking_smoothedcoef(ages = c(20:80),prs_quants = prs_quants,start = "Health",stop = "Cad",modelfit = fixedsmoke)
p=multipleprsfunc(s = s2[,,,1],prsprobs = pnorm(prs_quants))
ages=c(20:80)
m=matriskfunc_coef(p,ages,quantiles = prs_quants)


agesint=seq(30,50,by=5)
ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)

for(i in 1:length(agesint)){
  age=agesint[i]
  ten.year[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart =age,agestop = age+10)
  lifetime[i,]=projection_withmat(survivalmat = m$yearlynotrisk,agestart = age,agestop = 80)
}

library(plyr)

emp.ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
emp.lifetime=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)
atrisk=firstmeasure
firstmeasure$round5=round_any(firstmeasure$Age,5,f = round)
for(i in 1:length((agesint))){
  age=agesint[i]
  for(j in 1:length(levels(firstmeasure$int))){
    
    cat=levels(firstmeasure$int)[j]
    atrisk=firstmeasure[firstmeasure$int==cat&firstmeasure$round5==age,]
    print(dim(atrisk))
    print(age)
    print(head(atrisk))
    emp.ten.year[i,j]=mean(atrisk$chd==1&atrisk$chddate<3650)
    
    emp.lifetime[i,j]=mean(atrisk$chd==1)
  }}



ascvd.ten.year=matrix(NA,nrow = length(agesint),ncol=length(prs_quants)*2)


for(i in 1:length((agesint))){
  age=agesint[i]
  for(j in 1:length(levels(firstmeasure$int))){
    cat=levels(firstmeasure$int)[j]
    atrisk=firstmeasure[firstmeasure$int==cat&firstmeasure$round5==age,]
    print(dim(atrisk))
    
    ascvd.ten.year[i,j]=mean(atrisk$ASCVD)
    #ascvdriskmat[i,2]=mean(atrisk$chd==1)
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

diff.mstate.life=abs(data.frame(lifetime-emp.lifetime))
d=as.matrix(diff.mstate.life)
sqrt(mean(d^2))

diff.ascvd.life=abs(data.frame(ascvd.ten.year/100-emp.lifetime))
d=as.matrix(diff.ascvd.life)
sqrt(mean(d^2))


mat=cbind(pcten,msten,pclif,mslif)
colnames(mat)=c("PCE.10","MSgene.10","PCE.life","MSgene.life")
d=data.frame(mat)
d=t(d)
d=data.frame(d)
d$score=rownames(d)
d$score=as.factor(d$score)
levels(d$score)=c("PCE.10","MSgene.10","PCE.life","MSgene.life")

p <- d %>%
  mutate(score = fct_relevel(score, 
                             "PCE.10","MSgene.10","PCE.life","MSgene.life")) %>%
  ggplot(aes(x=score,y= d,
          ymin=d-0.12,
          ymax=d+0.12,
          fill=as.factor(score))) +
  geom_bar(stat="identity") + geom_errorbar(colour="black")+labs(y="RMSE, FOS",x="Age of first Interaction",
                                                                 fill="Score")+theme_classic(base_size = 15)
                                                              

  saveRDS(p,file = "fosrmse.rds")


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
diff.ascvd$score=rep("PCE",5)
diff.ascvd$age=agesint

diff.mstate$se=sd(as.matrix(sqrt(diff.mstate^2)))
diff.mstate$score=rep("MSGene",5)
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



#t.test(x = rm[c(1:7),c(1:3)],r[c(8:14),c(1:3)])

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
  geom_errorbar( position = position_dodge(), colour="black") +labs(y="RMSE 10 year risk",x="Age",fill="Genomic Level: Score")+theme_classic(base_size = 20)+ylim(c(0,0.20))
#geom_point(position=position_dodge(.9), aes(y=value, colour=interaction))


# Create life RMSE

diff.ascvd.life$se=sd(as.matrix(sqrt(diff.ascvd.life^2)))
diff.ascvd.life$score=rep("PCE",5)
diff.ascvd.life$age=agesint

diff.mstate.life$se=sd(as.matrix(sqrt(diff.mstate.life^2)))
diff.mstate.life$score=rep("MSGene",5)
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
  geom_errorbar( position = position_dodge(), colour="black") +labs(y="RMSE Lifetme risk",x="Age",fill="Genomic Level: Score")+theme_classic(base_size = 20)+ylim(c(0,0.20))
#geom_point(position=position_dodge(.9), aes(y=value, colour=interaction))



### 
