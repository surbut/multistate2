## Sensitivity analysis 



mpce=readRDS("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/output/mpce.rds")

## show phentoyping good

## check out phenotyping 

mpce$HyperLipidemia=factor(mpce$HyperLip_0_Any,levels = c(1,2),labels = c("control","case"))
mpce$Hypertension=factor(mpce$Ht_0_Any,levels = c(1,2),labels = c("control","case"))

g1=ggplot(mpce,aes(mpce$ldladj,fill=HyperLipidemia))+
  geom_density()+labs(x="Baseline LDL Cholesterol (ng/dl)",y="Density"
)+theme_classic()

g2=ggplot(mpce,aes(sbp,fill=Hypertension))+geom_density()+theme_classic()+
  labs(x="Baseline Systolic Blood Pressure (mmHg)",y="Density"
)

library(ggpubr)

ggsave(ggarrange(g1,g2,nrow=1),file="~/multistate2/rebuttal/baseline_phenotype.png",width = 15)



## longitudinal measuremments


sbp=readRDS("~/Library/CloudStorage/Dropbox-Personal/denaxasvalues/sbp.rds")
ldl=readRDS("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/output/ldl_gp.rds")
tc=readRDS("~/Library/CloudStorage/Dropbox-Personal/denaxasvalues/tc.rds")
  
hyperlips=mpce[mpce$HyperLip_0_Any==2,"identifier"]
hts=mpce[mpce$Ht_0_Any==2,"identifier"]

ldl$case=ifelse(ldl$eid%in%hyperlips$identifier,"Case","Control")
sbp$case=ifelse(sbp$eid%in%hts$identifier,"Case","Control")
tc$case=ifelse(tc$eid%in%hyperlips$identifier,"Case","Control")



g1=ggplot(ldl,aes(value1,fill=as.factor(case)))+geom_density()+
  labs(x="LDL-C Longitudinal Values (mg/dl)",fill="Disease Status")+theme_classic()

g2=ggplot(tc,aes(value1,fill=as.factor(case)))+geom_density()+
  labs(x="Total Cholesterol Longitudinal Values (mg/dl)",fill="Disease Status")+theme_classic()

ggsave(ggarrange(g1,g2,nrow=1),file="~/multistate2/rebuttal/allmeas.png",width = 15)




### number of measuremets

lc=ldl%>%group_by(eid)%>%summarise(n=length(event_dt))
sc=sbp%>%group_by(eid)%>%summarise(n=length(event_dt))

lc$case=ifelse(lc$eid%in%hyperlips$identifier,1,0)
sc$case=ifelse(sc$eid%in%hts$identifier,1,0)

n1=ggplot(lc,aes(n,fill=as.factor(case)))+geom_density(adjust=5)+
  labs(x="Number of LDL Measurements",fill="Disease Status")+theme_classic()
n2=ggplot(sc,aes(n,fill=as.factor(case)))+geom_density(adjust=5)+
  labs(x="Number of SBP Measurements",fill="Disease Status")+theme_classic()

ggsave(ggarrange(n1,n2,nrow=1),file="~/multistate2/rebuttal/nummeas.png",width = 15)


### wshow that these measruements are similar depending on age of enrollemnet

load("~/Library/CloudStorage/Dropbox-Personal///pheno_dir/output/merged_pheno_censor_final_withdrugs_smoke.rds")
dfh$cad.prs.lec=cut(dfh$cad.prs,breaks = c(-5,-0.84,0.84,5),labels = c("low","mid","high"))
dfh$int=interaction(dfh$f.31.0.0,dfh$cad.prs.lec)
# Relabel the levels of the interaction variable
levels(dfh$int) <- c(1,2,3,4,5,6)
train=dfh[1:(nrow(dfh)*0.80),]
dfascvd=readRDS("~/multistate2//output/dfascvd_newbp.rds")
test=dfh[!(dfh$identifier%in%train$identifier),]
df=merge(test,dfascvd[,-which(names(dfascvd)%in%c("age","anylipidmed0","bp_med2","smoke"))],by.x="identifier",by.y="sample_id")

df$round_age=as.numeric(round(df$enrollage/10,0)*10)
# Assuming df is your dataframe with age, sbp, and choladj columns
# Filtering the dataframe for individuals at ages 40, 50, 60, and 70
df_filtered <- df[df$round_age %in% c(40, 50, 60, 70), ]

# Plotting cholesterol at baseline for ages 40, 50, 60, and 70
hist_total_chol_by_age=
  ggplot(df_filtered[df_filtered$statin==0,], aes(x = choladj,fill=as.factor(round_age))) +
  geom_histogram(bins = 100, color = "black") +
  facet_wrap(~ round_age, scales = "free_y") +
  labs(title = "",
       x = "Total Cholesterol (mg/dl) at Enrollment",
       y = "Count",fill="Age Group")+
  theme_minimal()

hist_sbp_by_age=
  ggplot(df_filtered[df_filtered$antihtn==0,], aes(x = sbp,fill=as.factor(round_age))) +
  geom_histogram(bins = 100, color = "black") +
  facet_wrap(~ round_age, scales = "free_y") +
  labs(title = "",
       x = "Systolic Blood Pressure (mmHg) at Enrollment",
       y = "Count",fill="Age Group")+
  theme_minimal()


ggsave(ggarrange(hist_total_chol_by_age,hist_sbp_by_age,nrow=1),
       file="~/multistate2/rebuttal/hist_by_age.png",width = 15)


### 

baseline=readRDS("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/output/dfukb_baseline_pheno.rds")
baseline$enrollage=difftime(baseline$f.53.0.0,baseline$reference_date,units = "days")/365.25
baseline$indata=ifelse(baseline$identifier%in%dfh$identifier,1,0)
gae=fread("~/Library/CloudStorage/Dropbox-Personal/ukb.kgp_projected.tsv.gz")
baseline$identifier=as.character(baseline$identifier)
baseline=merge(baseline,gae[,c("eid","rf80")],by.x="identifier",by.y="eid",all.x = T)
baseline$race=factor(baseline$rf80)
dat=baseline

library(table1)



dat$Sex=ifelse(dat$f.31.0.0=="1","Male","Female")
#dat$as2[is.na(dat$as2)]=median(!is.na(dat$as2))
dat$Birthdate=dat$f.34.0.0
dat$RecruitAge=as.numeric(dat$enrollage)
dat$Race=as.factor(dat$race)

library(table1)
label(dat$f.31.0.0) <- "Sex"


f=table1(~ Sex + Birthdate + RecruitAge+Race|indata, data=dat)


$### trait summary
  
df_definitions <- fread("~/Downloads//definitions_cardiometabolic_traits.tsv")

# Filter for specific traits
traits_of_interest <- c('Dm', 'Cad', 'Ht', 'HyperLip')
df=dfDefinitions_processed_expanded[dfDefinitions_processed_expanded$TRAIT%in%traits_of_interest,c("TRAIT","DESCRIPTION","TS","ICD10","ICD9","OPCS4","READ2","CTV3","f.20002","f.20004")] 
write.table(df ,"~/Library/CloudStorage/Dropbox-Personal/long_definitions.tsv",sep = "\t",quote = F,row.names = FALSE)

write.table(melt(df,id.vars=c("TRAIT","DESCRIPTION")),"~/Library/CloudStorage/Dropbox-Personal/long_definitions_melt.csv",sep = "\t",quote = F,row.names = FALSE)
# Select and rename columns for a comprehensive summary
library(dplyr)

# Print the summary table
print(df_summary)

### densityof prs 

# Create a new column for age group by decade
library(ggridges)
dfh$prs=scale(dfh$cad.prs)
dfh$age_group <- cut(as.numeric(dfh$enrollage),
                     breaks = c(39, 49, 59, 69, 79,85),
                     labels = c("40s", "50s", "60s", "70s","80s"),
                     include.lowest = TRUE)

# Plot the distribution of PRS by age group
a1=ggplot(dfh[!is.na(dfh$age_group),], aes(x = prs, y = age_group, fill = age_group)) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_brewer(palette = "Spectral") +
  labs(title = "Distribution of PRS by Age at Enrollment",
       x = "Polygenic Risk Score (PRS)",
       y = "Age Group") +
  theme_ridges() +
  theme(legend.position = "none")  # Remove legend if age_group is used in fill


a2=ggplot() +
  geom_qq(aes(sample = dfh$prs[dfh$age_group%in%"40s"]), color = "blue") +
  geom_qq_line(aes(sample = dfh$prs[dfh$age_group%in%"40s"]), color = "blue") +
  geom_qq(aes(sample = dfh$prs[dfh$age_group%in%"70s"]), color = "red") +
  geom_qq_line(aes(sample = dfh$prs[dfh$age_group%in%"70s"]), color = "red") +
  labs(title = "QQ Plot Comparing PRS for Enrollment at 40 vs 70",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")+theme_classic()


ggsave(ggarrange(a1,a2,nrow=1),
       file="~/multistate2/rebuttal/prs_by_age.png",width = 15)

### Ethnicity
library(tibble)
library(dplyr)

## see updated Fig4.html

extdata_dir <- paste0(system.file("extdata", package="ukbpheno"),"/")
fdata_setting <- paste0(extdata_dir,"data.settings.tsv")

fdefinitions="~/Library/CloudStorage/Dropbox-Personal/definitions_cardiometabolic_traits2.tsv"
dfDefinitions_processed_expanded<-read_definition_table(fdefinitions,fdata_setting,extdata_dir)

dfukb_baseline_pheno=readRDS("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/output/dfukb_baseline.rds")
diseases<-dfDefinitions_processed_expanded$TRAIT[c(1,2,3,7,8,9,13,14,15,24,41,42,c(43:81))]

lst.harmonized.data=readRDS("~/Library/CloudStorage/Dropbox-Personal//pheno_dir/output/lst.harmonized.data.updated.rds")

fdata_setting <- paste0(extdata_dir,"data.settings.tsv")
dfData.settings <- fread(fdata_setting)

d4<-lst.harmonized.data$dfukb[,c("identifier","f.52.0.0","f.34.0.0")]
# # # f.34.0.0 contains year of birth, f.52 is month of birth, create birthday on 15th of month
lst.harmonized.data$dfukb$Birthdate<-as.Date(with(d4,paste(f.34.0.0,f.52.0.0,15,sep="-")),"%Y-%m-%d")
#
#df_reference_dt_v0<-lst.harmonized.data$dfukb[,c("identifier","Birthdate")]
df_reference_dt_v0<-lst.harmonized.data$dfukb[,c("identifier","f.53.0.0")]

trait="Ht"
make_upsetplot(definition=dfDefinitions_processed_expanded%>%filter(TRAIT==trait),
               lst.harmonized.data$lst.data,dfData.settings,df.reference.dates = df_reference_dt_v0)

trait="Cad"
make_upsetplot(definition=dfDefinitions_processed_expanded%>%filter(TRAIT==trait),
               lst.harmonized.data$lst.data,dfData.settings,df.reference.dates = df_reference_dt_v0)

### 

baseline_data=readRDS("~/Library/CloudStorage/Dropbox-Personal//pheno_dir/output/baseline_sbp_ldl_hdl_tc.rds")
colnames(baseline_data)=c("eid","sbp","ldl","hdl","tc","baseline_date")
m=na.omit(melt(baseline_data,id.vars = c("eid","baseline_date")))


sbp=readRDS("~/Library/CloudStorage/Dropbox-Personal/denaxasvalues/sbp.rds")[,c("eid","event_dt","value1")]
ldl=readRDS("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/output/ldl_gp.rds")[,c("eid","event_dt","value1")]
tc=readRDS("~/Library/CloudStorage/Dropbox-Personal/denaxasvalues/tc.rds")[,c("eid","event_dt","value1")]
hdl=readRDS("~/Library/CloudStorage/Dropbox-Personal/denaxasvalues/hdl.rds")[,c("eid","event_dt","value1")]

# Add a new column 'test_type' to each dataset
sbp$test_type="sbp" 
ldl$test_type= "ldl"
hdl$test_type= "hdl"
tc$test_type= "tc"

# Combine the datasets
gp_data <- rbind(sbp, ldl, tc,hdl)

# Format the 'event_dt' column as a Date object
gp_data[, event_dt := as.Date(event_dt)]

library(dplyr)

gp_data <- bind_rows(sbp, ldl,tc, hdl)

library(data.table)
setDT(baseline_data)

m$baseline_date=as.Date(m$baseline_date)
# Function to process and compare a specific lab test
process_and_compare_lab_test <- function(lab, m, gp_data) {
  # Filter GP data for the specific test type
  specific_test_data <- gp_data[test_type==lab,]
  baseline_dat=m[m$variable==lab,]
  
  # For each patient and test type, calculate the average of GP measurements within one year of baseline
  specific_test_data[, event_dt := as.Date(event_dt)]
  baseline_dat$baseline_date=as.Date(baseline_dat$baseline_date)
 
  
  #Merge GP data with baseline data based on patient identifier
  merged_data <- na.omit(merge(specific_test_data,baseline_dat,  by.x = "eid", by.y = "eid", all.x = TRUE, suffixes = c("_baseline", "_gp")))
  merged_data$abstime=as.numeric(abs(merged_data$event_dt-merged_data$baseline_date))
  
  averages_within_year <- merged_data%>%group_by(eid)%>%filter(abstime<180)%>%summarise(amb=median(value1),baseline=mean(value))
  
  #averages_within_year <- merged_data%>%group_by(eid)%>%summarise(amb=median(value1),
  averages_within_year$diff=averages_within_year$amb-averages_within_year$baseline
  
  averages_within_year=averages_within_year[!is.na(averages_within_year$diff),]
   # Perform statistical analysis (e.g., paired t-test)
  t.test.result=t.test(averages_within_year$amb,averages_within_year$baseline,paired = T)
  
  # Return results
  list(mean_difference = mean(averages_within_year$diff, na.rm = TRUE), t_test_result = t.test.result$estimate)
}

# Example usage for LDL
sbp_results <- process_and_compare_lab_test("sbp", m, gp_data)
ldl_results <- process_and_compare_lab_test("ldl", m, gp_data)
hdl_results <- process_and_compare_lab_test("hdl", m, gp_data)
tc_results <- process_and_compare_lab_test("tc", m, gp_data)


## population sd
sd(na.omit(baseline_data$sbp))
sd(na.omit(baseline_data$hdl))
sd(na.omit(baseline_data$tc))


### differences between lab and ukb
cases_ph=lst.case_control$all_event_dt.Include_in_cases$identifier
i=setdiff(cases_ph,cadlab$sample_id[cadlab$has_disease==1])
summary(as.factor(lst.case_control$all_event_dt.Include_in_cases[lst.case_control$all_event_dt.Include_in_cases$identifier%in%i,code]))



### whoing BMI not sueful

train2=merge(dfh,dfascvd[,-which(names(dfascvd)%in%c("age","anylipidmed0","bp_med2","smoke"))],by.x="identifier",by.y="sample_id")
train2$phenos.enrollment=train2$f.21003.0.0
train2$bmi=train2$f.21001.0.0

mod_bmi=fitfunc2(data.table(train2),ages = ages,nstates = nstates,mode = "binomial",covariates ="cad.prs+f.31.0.0+smoke+antihtn_now+statin_now+bmi+choladj+hdladj")
b=coefplotsmooth2(ages = ages,start = "Health",stop = "Cad",modelfit = mod_bmi,window_width = 20,span = 0.75,degree = 2)
coefs_mod4=b$custom_smooth


#### reassurance with HR

# in practice
atrisk = expand.grid(intercept = intercept, 
                     CAD_PRS = cad.prs, 
                     Sex = sex, 
                     Smoke = smoke, 
                     AntiHTN_Now = antihtn_now, 
                     Statin_Now = statin_now)

# Display the first few rows to check the structure
head(atrisk)

predictedrisks=readRDS("~/output/predictedrisk_forfixedmatrix.rds")

mat=cbind(atrisk[c(1:20),c(1:4)],round(predictedrisks[c(1:20),"40",c("Health","Ht","Dm")],2))
