## for table1

nstates=c("Health", "Ht","HyperLip","Dm","Cad","death","Ht&HyperLip","HyperLip&Dm","Ht&Dm","Ht&HyperLip&Dm")
gpc=readRDS("~/Library/CloudStorage/Dropbox-Personal//ukbb-ehr-data/data/gp_event.rds")
g2=gpc[!duplicated(gpc$eid),]
saveRDS(g2,"~/Library/CloudStorage/Dropbox-Personal/g2.rds")
pdf("output/barplotyears.pdf")
barplot(table(year(g2$event_dt)),las=2,xlab="First Measurement Year",ylab="Number of Individuals")
dev.off()

load("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/output/merged_pheno_censor_final_withdrugs_smoke.rds")
df_final=dfh
dim(dfh)
#dim(lst.harmonized.data$dfukb)[1]-dim(dfh)[1]
a=merge(g2,df_final,by.x="eid",by.y="identifier")

a=data.frame(a)
a$agerec=difftime(a$event_dt,a$Birthdate,units="day")/365.25
#a$Durationfollowed=as.numeric(a$Death_Censor_Age)-as.numeric(a$agerec)
a$Durationfollowed=as.numeric(a$Death_Censor_Age)-as.numeric(a$agerec)
saveRDS(a,"../output/fortable1.rds")
summary(a$Durationfollowed)

df_final$cad.prs=scale(df_final$cad.prs)


library(table1)
dat=readRDS("../output/fortable1.rds")
df_ascvd=readRDS("~/multistate2//output/dfascvd_newbp.rds")
dat=merge(dat,df_ascvd[,c("sample_id","Race")],by.x="eid",by.y = "sample_id",all.x=T)

dat$gpmemb=ifelse(dat$eid%in%g2$eid,1,0)
dat$gpmemb=factor(dat$gpmemb,levels = c(0,1),labels = c("Not Member","Member"))
dat$race=as.numeric(ifelse(dat$Race=="white",1,0))
dat$smoke=as.factor(dat$smoke)
dat$gpm=as.factor(dat$gpmemb)

dat$cad.prs=scale(dat$cad.prs)
dat$cad.prs.lev=cut(dat$cad.prs,breaks = c(-10,qnorm(0.20),qnorm(0.80),10),labels=c("Low","Intermediate","High"))
dat$Durationfollowed=as.numeric(dat$Death_Censor_Age)-as.numeric(40)
#dat=merge(dat,df_ascvd[,c("sample_id","as2")],by.x="eid",by.y = "sample_id",all.x=TRUE)
dat=dat[-which(dat$Cad_0_censor_age<40),]
dat$agerec=as.numeric(dat$agerec)
dat$agerec[dat$agerec<18]=18
dat$f.31.0.0=ifelse(dat$f.31.0.0=="1","Male","Female")
#dat$as2[is.na(dat$as2)]=median(!is.na(dat$as2))
dat$Birthdate=year(dat$Birthdate)
dat$Ht_0_Any=as.factor(ifelse(dat$Ht_0_Any==2,1,0))
dat$HyperLip_0_Any=as.factor(ifelse(dat$HyperLip_0_Any==2,1,0))
dat$Dm_0_Any=as.factor(ifelse(dat$Dm_0_Any==2,1,0))
dat$Cad_0_Any=as.factor(ifelse(dat$Cad_0_Any==2,1,0))
dat$smoke=as.factor(dat$smoke)
dat
library(table1)
label(dat$f.31.0.0) <- "Sex"
label(dat$antihtn) <- "Start an anti-Hypertensive"
label(dat$smoke) <- "Current Smoker"
label(dat$Ht_0_Any) <- "Develop Hypertension"
label(dat$Dm_0_Any) <- "Develop Diabetes"
label(dat$HyperLip_0_Any) <- "Develop Hyperlipidemia"
label(dat$Cad_0_Any) <- "Develop Coronary Disease"
#label(dat$as2) <- "ASCVD at Enrollment"
label(dat$agerec) <- "Age First Enrolled in NHS"
label(dat$Durationfollowed) <- "Years Followed"



f=table1(~ f.31.0.0 + agerec+Durationfollowed+Birthdate+Ht_0_Any+Cad_0_Any+HyperLip_0_Any+antihtn+smoke+race+gpmemb|cad.prs.lev, data=dat)

library(xtable)

install.packages(c("officer", "flextable"))
library(officer)
library(flextable)
doc <- read_docx()

ft <- flextable(as.data.frame(f))

doc <- body_add_flextable(doc, ft)

print(doc, target = "../output/output_file.docx")
print(xtable(as.data.frame(f)),include.rownames=F)


###


library(khsmisc)
library(khsverse)


i=intersect(g2$eid,dfh$identifier)
nrow(dfh)-length(i)


dfascvd=readRDS("~/multistate2//output/dfascvd_newbp.rds")
test=dfh[!(dfh$identifier%in%train$identifier),]
test=merge(test,dfascvd[,-which(names(dfascvd)%in%c("age","anylipidmed0","bp_med2","smoke"))],by.x="identifier",by.y="sample_id")
  
traindim=221351-nrow(test)

  
length(intersect(dfh$identifier,dfascvd$sample_id))
#400552 
design <- tibble::tribble(
  ~left,~n_left, ~right,~n_right,
  "Study base with Outcome Data", 502461 ,"Lack QC Genotype, Sex, Birthdate or Smoking information", 20534,
  "Contain baseline covariates",  481927,    "Not present in GP records", 259930,
  "Individuals in GP Records", 221997,   "With CAD at baseline",   646,
  "PRS, Pheno, Covariate info", 221351,  "",  NA_integer_,
  "Training Set", 142234, "", NA_integer_,
  "Testing Set", 79117, "", NA_integer_
)



##

e=exclusion_flowchart(design, width = 2)
export_svg(
  
e %>%
  export_svg() %>%
  read_xml() %>%
  write_xml("../output/flowchart_msgene.svg")
###


length(intersect(dfh$identifier,dfascvd$sample_id))
#400552 

design <- tibble::tribble(
  ~left,~n_left, ~right,~n_right,
  "Study base with Outcome Data", 502461 ,"Lack QC Genotype, Sex, Birthdate or Smoking information", 20534,
  "Contain baseline covariates",  481927, "CAD at baseline", 1289,
  "CAD Free at 40 for model fit", 480638 , "80% for training", 384510,
  "Remain for testing", 96127,   "lacking TC, SBP or HDL for comparison with FRS",   17010,
  "Available for testing", 79117 ,  "",  NA_integer_,
 
)

e=exclusion_flowchart(design, width = 2)
export_svg(
  
e %>%
    export_svg() %>%
    read_xml() %>%
    write_xml("../output/flowchart_msgenereal.svg")
  
###

load("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/output/merged_pheno_censor_final_withdrugs_smoke.rds")
df_final=dfh
dim(dfh)
#dim(lst.harmonized.data$dfukb)[1]-dim(dfh)[1]
a=merge(g2,df_final,by.x="eid",by.y="identifier",all.y = T)

a=data.frame(a)
a$agerec=difftime(a$event_dt,a$Birthdate,units="day")/365.25
a$Durationfollowed=as.numeric(a$Death_Censor_Age)-as.numeric(a$agerec)
#saveRDS(a,"../output/fortable1.rds")

df_final$cad.prs=scale(df_final$cad.prs)


library(table1)
dat=a
df_ascvd=readRDS("~/multistate2//output/dfascvd_newbp.rds")
dat$cad.prs=scale(dat$cad.prs)
g2=readRDS("~/multistate2/output//g2.rds")
dat$cad.prs.lev=cut(dat$cad.prs,breaks = c(-10,qnorm(0.20),qnorm(0.80),10),labels=c("Low","Intermediate","High"))
dat$Durationfollowed=as.numeric(dat$Death_Censor_Age)-as.numeric(40)
#dat=merge(dat,df_ascvd[,c("sample_id","as2")],by.x="eid",by.y = "sample_id",all.x=TRUE)
dat=dat[-which(dat$Cad_0_censor_age<40),]
dat$agerec=as.numeric(dat$agerec)
dat$agerec[dat$agerec<18]=18
dat$f.31.0.0=ifelse(dat$f.31.0.0=="1","Male","Female")
#dat$as2[is.na(dat$as2)]=median(!is.na(dat$as2))
dat$Birthdate=year(dat$Birthdate)
dat$Ht_0_Any=as.factor(ifelse(dat$Ht_0_Any==2,1,0))
dat$HyperLip_0_Any=as.factor(ifelse(dat$HyperLip_0_Any==2,1,0))
dat$Dm_0_Any=as.factor(ifelse(dat$Dm_0_Any==2,1,0))
dat$Cad_0_Any=as.factor(ifelse(dat$Cad_0_Any==2,1,0))
dat$smoke=as.factor(dat$smoke)

dat=merge(dat,df_ascvd[,c("sample_id","Race")],by.x="eid",by.y = "sample_id",all.x=T)
#dat=dat[-which(dat$Cad_0_censor_age<40),]
dat$Race=as.factor(dat$Race)
dat$Race[is.na(dat$Race)]="white"
dat$race=as.numeric(ifelse(dat$Race=="white",1,0))
dat$race=as.factor(dat$race)
dat$gpmemb=ifelse(dat$eid%in%g2$eid,1,0)
dat$gpmemb=factor(dat$gpmemb,levels = c(0,1),labels = c("Not Member","Member"))
dat$smoke=as.factor(dat$smoke)
dat$gpm=as.factor(dat$gpmemb)

library(table1)
label(dat$f.31.0.0) <- "Sex"
label(dat$antihtn) <- "Start an anti-Hypertensive"
label(dat$smoke) <- "Current Smoker"
label(dat$Ht_0_Any) <- "Develop Hypertension"
label(dat$Dm_0_Any) <- "Develop Diabetes"
label(dat$HyperLip_0_Any) <- "Develop Hyperlipidemia"
label(dat$Cad_0_Any) <- "Develop Coronary Disease"
#label(dat$as2) <- "ASCVD at Enrollment"
label(dat$agerec) <- "Age First Enrolled in NHS"
label(dat$Durationfollowed) <- "Years Followed"


library(table1)
label(dat$f.31.0.0) <- "Sex"
label(dat$antihtn) <- "Start an anti-Hypertensive"
label(dat$smoke) <- "Current Smoker"
label(dat$Ht_0_Any) <- "Develop Hypertension"
label(dat$Dm_0_Any) <- "Develop Diabetes"
label(dat$HyperLip_0_Any) <- "Develop Hyperlipidemia"
label(dat$Cad_0_Any) <- "Develop Coronary Disease"

label(dat$Ht_0_censor_age) <- "Age HT"
label(dat$Dm_0_censor_age) <- "Age DM"
label(dat$HyperLip_0_censor_age) <- "Age HL"
label(dat$Cad_0_censor_age) <- "Age CAD"

label(dat$agerec) <- "Age First Enrolled in NHS"
label(dat$Durationfollowed) <- "Years Followed"
label(dat$gpmemb) <- "GPDR Members"
label(dat$race) <- "Proportion White"

f=table1(~ f.31.0.0 + Birthdate+agerec+Durationfollowed+Ht_0_Any+Cad_0_Any+HyperLip_0_Any+smoke+race+gpmemb|cad.prs.lev, data=dat)


f=table1(~ f.31.0.0 + Birthdate+agerec+Durationfollowed+Ht_0_Any+Cad_0_Any+Dm_0_Any+HyperLip_0_Any+smoke+race+Ht_0_censor_age+Cad_0_censor_age+Dm_0_censor_age+HyperLip_0_censor_age|gpmemb, data=dat)


