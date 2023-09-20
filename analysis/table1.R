## for table1

nstates=c("Health", "Ht","HyperLip","Dm","Cad","death","Ht&HyperLip","HyperLip&Dm","Ht&Dm","Ht&HyperLip&Dm")
gpc=readRDS("~/Library/CloudStorage/Dropbox-Personal//ukbb-ehr-data/data/gp_event.rds")
g2=gpc[!duplicated(gpc$eid),]

pdf("output/barplotyears.pdf")
barplot(table(year(g2$event_dt)),las=2,xlab="First Measurement Year",ylab="Number of Individuals")
dev.off()

load("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/output/merged_pheno_censor_final_withdrugs_smoke.rds")
df_final=dfh
dim(dfh)
dim(lst.harmonized.data$dfukb)[1]-dim(dfh)[1]
a=merge(g2,df_final,by.x="eid",by.y="identifier")

a=data.frame(a)
a$agerec=difftime(a$event_dt,a$Birthdate,units="day")/365.25
a$Durationfollowed=as.numeric(a$Death_Censor_Age)-as.numeric(a$agerec)
saveRDS(a,"../output/fortable1.rds")

df_final$cad.prs=scale(df_final$cad.prs)


library(table1)
dat=readRDS("../output/fortable1.rds")

dat=merge(dat,df_ascvd[,c("sample_id","as2")],by.x="eid",by.y = "sample_id",all.x=TRUE)
dat=dat[-which(dat$Cad_0_censor_age<40),]
dat$agerec=as.numeric(dat$agerec)
dat$agerec[dat$agerec<18]=18
dat$f.31.0.0=ifelse(dat$f.31.0.0=="1","Male","Female")
dat$as2[is.na(dat$as2)]=median(!is.na(dat$as2))
dat$Birthdate=year(dat$Birthdate)
dat$Ht_0_Any=as.factor(ifelse(dat$Ht_0_Any==2,1,0))
dat$HyperLip_0_Any=as.factor(ifelse(dat$HyperLip_0_Any==2,1,0))
dat$Dm_0_Any=as.factor(ifelse(dat$Dm_0_Any==2,1,0))
dat$Cad_0_Any=as.factor(ifelse(dat$Cad_0_Any==2,1,0))
dat$cad.prs.lev=factor(dat$cad.prs.lev,levels = c("low","mid","high"),labels=c("Low","Intermediate","High"))

library(table1)
label(dat$f.31.0.0) <- "Sex"
label(dat$antihtn) <- "Start an anti-Hypertensive"
label(dat$smoke) <- "Current Smoker"
label(dat$Ht_0_Any) <- "Develop Hypertension"
label(dat$Dm_0_Any) <- "Develop Diabetes"
label(dat$HyperLip_0_Any) <- "Develop Hyperlipidemia"
label(dat$Cad_0_Any) <- "Develop Coronary Disease"
label(dat$as2) <- "ASCVD at Enrollment"
label(dat$agerec) <- "Age First Enrolled in NHS"
label(dat$Durationfollowed) <- "Years Followed"

f=table1(~ f.31.0.0 + agerec+Durationfollowed+Birthdate+Ht_0_Any+Cad_0_Any+HyperLip_0_Any+antihtn+smoke|cad.prs.lev, data=dat)

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

design <- tibble::tribble(
  ~left,               ~n_left, ~right,              ~n_right,
  "Study base with Outcome Data", 502,461 , "Lack QC Genotype, Sex, Birthdate or smoking information", 20534,
  "Contain baseline covariates",  481,927,    "Excluded from gp clinical atlas", 259930,
  "Individuals in GP Clinical", 221997,   "With CAD at baseline",   646,
  "PRS, Pheno, Covariate info", 221351,  "",                  NA_integer_)


##

e=exclusion_flowchart(design, width = 2)
e
###

