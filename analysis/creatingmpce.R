#rm(list=ls())
library('data.table')
install.packages("CVrisk")
library(CVrisk)
library("dplyr")
#######
d2=fread("~/Dropbox/output/Diabetes_Type_2.tab.tsv.gz")
d2[d2$prevalent_disease==1,"incident_disease"]=0
d3=data.frame(phenos.enrollment=d2$enroll_age,sample_id=d2$sample_id,prev_disease_dm=d2$prevalent_disease,inc_disease_dm=d2$incident_disease)




dfukb_baseline=readRDS("~/dynamichr/output/dfukbaseline.rds")
dat2=merge(x = d3,y = dfukb_baseline,by.x = "sample_id",by.y = "identifier",all.x = TRUE)
dat2$bp_med2=ifelse(dat2$f.6153.0.0==2|dat2$f.6153.0.1==2|dat2$f.6153.0.2==2|dat2$f.6177.0.0==2|dat2$f.6177.0.1==2|dat2$f.6177.0.2==2,1,0)
dat2$bp_med2[which(is.na(dat2$bp_med2))]=0

### lipids info
lipids=fread("~/Dropbox/ukbb-lipids-meds.txt")
sum(rowSums(is.na(lipids[,c("eid","hdladj","choladj","anylipidmed0")]))>0)
# 72721

lip=na.omit(lipids[,c("eid","statin0","ldladj","hdladj","choladj","anylipidmed0")])

d3=merge(dat2,lip,by.x="sample_id",by.y="eid")

#ukb=fread("~/Dropbox/big_ukb_file.txt",header=T,sep="\t")
cov=data.frame("age"=ukb$age,"iid"=ukb$id,"sex"=ukb$Sex,"sbp"=ukb$SBP_adjMeds,"dbp"=ukb$DBP_adjMeds,"smoke"=ukb$SmokingStatus,"Race"=ukb$in_white_British_ancestry_subset)
cov$Race=ifelse(cov$Race==1,"white","other")

sum(rowSums(is.na(cov))>0)
## exclude 41192 without baseling info

dim(cov)
cov=na.omit(cov)
dim(na.omit(cov))


df=merge(d3,cov,by.x="sample_id",by.y="iid")
dim(df)

#401116  



## compute ascvd scores using better BP metric

df$sex=ifelse(as.factor(df$sex)=="Male","male","female")
df$old_smoke=df$smoke
df$smoke=ifelse(df$smoke==2,1,0)
rownames(df)=as.character(df$sample_id)

df2=compute_CVrisk(df,scores = c("ascvd_10y_accaha"),
                   age = "phenos.enrollment", race = "Race", gender = "sex", bmi = "bmi", sbp = "sbp",hdl = "hdladj", totchol = "choladj", bp_med = "bp_med2", smoker = "smoke",
                   diabetes = "prev_disease_dm", lipid_med = "anylipidmed0",
                   fh_heartattack = NULL, cac = NULL)

saveRDS(df2,"~/multistate/output/dfascvd_newbp.rds")




###

library('data.table')
library("CVrisk")
library("dplyr")
#######
d2=fread("~/Dropbox/output/Diabetes_Type_2.tab.tsv.gz")
d2[d2$prevalent_disease==1,"incident_disease"]=0
d3=data.frame(phenos.enrollment=d2$enroll_age,sample_id=d2$sample_id,prev_disease_dm=d2$prevalent_disease,inc_disease_dm=d2$incident_disease)



dfukb_baseline=readRDS("~/dynamichr/output/dfukbaseline.rds")
dat2=merge(x = d3,y = dfukb_baseline,by.x = "sample_id",by.y = "identifier",all.x = TRUE)
dat2$bp_med2=ifelse(dat2$f.6153.0.0==2|dat2$f.6153.0.1==2|dat2$f.6153.0.2==2|dat2$f.6177.0.0==2|dat2$f.6177.0.1==2|dat2$f.6177.0.2==2,1,0)
dat2$bp_med2[which(is.na(dat2$bp_med2))]=0
dat2$sbp=dat2$f.4080.0.0
dat2$sbp[which(is.na(dat2$f.4080.0.0&!is.na(dat2$f.4080.0.1)))]=dat2$f.4080.0.1[which(is.na(dat2$f.4080.0.0&!is.na(dat2$f.4080.0.1)))]
dat2$sex=ifelse(dat2$f.31.0.0==1,"male","female")
dat2$age=dat2$f.21003.0.0
dat2$smoke=ifelse(dat2$f.20116.0.0==2,1,0)
### lipids info
lipids=fread("~/Dropbox/ukbb-lipids-meds.txt")
sum(rowSums(is.na(lipids[,c("eid","hdladj","choladj","anylipidmed0")]))>0)
# 72721

lip=na.omit(lipids[,c("eid","statin0","ldladj","hdladj","choladj","anylipidmed0")])

d3=merge(dat2,lip,by.x="sample_id",by.y="eid")

#ukb=fread("~/Dropbox/big_ukb_file.txt",header=T,sep="\t")
cov=data.frame("iid"=ukb$id,"Race"=ukb$in_white_British_ancestry_subset)
cov$Race=ifelse(cov$Race==1,"white","other")

sum(rowSums(is.na(cov))>0)
## exclude 14255 without baseling info

dim(cov)
cov=na.omit(cov)
dim(na.omit(cov))


df=merge(d3,cov,by.x="sample_id",by.y="iid")
dim(df)
df=df[,c("sample_id","prev_disease_dm","sex","age","bp_med2","sbp","smoke","hdladj","choladj","anylipidmed0","Race")]


df=na.omit(df)
#400969     10

## compute ascvd scores using better BP metric

# df$sex=ifelse(as.factor(df$sex)=="Male","male","female")
# df$old_smoke=df$smoke
# df$smoke=ifelse(df$smoke==2,1,0)
# rownames(df)=as.character(df$sample_id)

df2=compute_CVrisk2(df,scores = "as2",
                    age = "age", race = "Race", gender = "sex",
                    sbp = "sbp",hdl = "hdladj", totchol = "choladj", bp_med = "bp_med2",
                    smoker = "smoke",
                    diabetes = "prev_disease_dm", lipid_med = "anylipidmed0",
                    fh_heartattack = NULL, cac = NULL)

saveRDS(df2,"~/multistate/output/dfascvd_newbp.rds")

df_final = load("~/Dropbox/pheno_dir/output/merged_pheno_censor_final_withdrugs.rds")$dfh
#gpc=readRDS("~/Dropbox/ukbb-ehr-data/data/gp_event.rds")
#g2=gpc[!duplicated(gpc$eid),]
g2=readRDS("output/g2.rds")
a=merge(g2,df_final,by.x="eid",by.y="identifier")


a=data.frame(a)
a$agerec=difftime(a$event_dt,a$Birthdate,units="day")/365.25
a$Durationfollowed=as.numeric(a$Death_Censor_Age)-as.numeric(a$agerec)
head(a)
a=merge(g2,df_final,by.x="eid",by.y="identifier")
a=data.frame(a)
a$agerec=difftime(a$event_dt,a$Birthdate,units="day")/365.25
a$Durationfollowed=as.numeric(a$Death_Censor_Age)-as.numeric(a$agerec)
saveRDS(a,"output/fortable1.rds")



library("ukbpheno")
library("ggplot2")
library("survival")
library("dplyr")
library("tidyr")
# the directory with datafiles
pheno_dir="~/Dropbox (Personal)/pheno_dir/"
fukbtab <- paste(pheno_dir,"ukb47823.tab",sep="")
# meta data file
fhtml <- paste(pheno_dir,"ukb47823.html",sep="")
# hospital inpatient data
fhesin <- paste(pheno_dir,"ukb_showcase_9.01/hesin.txt",sep="")
fhesin_diag <- paste(pheno_dir,"ukb_showcase_9.01/hesin_diag.txt",sep="")
fhesin_oper <- paste(pheno_dir,"ukb_showcase_9.01/hesin_oper.txt",sep="")
dfhtml <- read_ukb_metadata(fhtml)
#
# # Rename the identifier column in the metadata
#
dfhtml[which(dfhtml$field.tab=="f.eid"),]$field.tab<-"identifier"
#
baseline_fields<-c(21003,31,6177,6153,20116,20160,4080)
#
dfukb_baseline <- read_ukb_tabdata(fukbtab,dfhtml,fields_to_keep = baseline_fields)
#
head(dfukb_baseline)
sum(dfukb_baseline$f.20116.0.0==2)
summary(as.factor(dfukb_baseline$f.20116.0.0))
table(as.factor(dfukb_baseline$f.20116.0.0))
table(as.factor(dfukb_baseline$f.20160.0.0))
saveRDS(dfukb_baseline,"~/Dropbox/pheno_dir/dfukbaseline_7.rds")

