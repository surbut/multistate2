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

baseline=readRDS("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/output/dfukb_baseline.rds")
baseline$indata=ifelse(baseline$identifier%in%dfh$identifier,1,0)
ethnicity=readRDS("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/ethnicity.rds")

baseline%>%group_by(indata)%>%summary(mean(as.factor(f.31.0.0),mean(f.21003.0.0),