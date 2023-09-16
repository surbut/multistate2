source("~/dynamichr/code/utils.R")
source("~/multistate2//code/smoothtest.R")
source("~/multistate2//code/newsmooth.R")
source("~/multistate2/code/fitarray.R")
library("reshape2")
source("~/multistate2/code/arrayindicate.R")

load("~/Library/CloudStorage/Dropbox-Personal///pheno_dir/output/merged_pheno_censor_final_withdrugs_smoke.rds")
dfh$cad.prs.lec=cut(dfh$cad.prs,breaks = c(-5,-0.84,0.84,5),labels = c("low","mid","high"))
dfh$int=interaction(dfh$f.31.0.0,dfh$cad.prs.lec)
# Relabel the levels of the interaction variable
levels(dfh$int) <- c(1,2,3,4,5,6)


train=dfh[1:(nrow(dfh)*0.80),]

dfascvd=readRDS("~/multistate2//output/dfascvd_newbp.rds")
test=dfh[!(dfh$identifier%in%train$identifier),]
test=merge(test,dfascvd[,-which(names(dfascvd)%in%c("age","anylipidmed0","bp_med2","smoke"))],by.x="identifier",by.y="sample_id")
test$ascvd_10y_accaha=test$as2
test$phenos.enrollment=test$f.21003.0.0
test=data.table(test)
source("~/multistate2/code/frs30_URBUT/fun.frs_30ynew.R")
library(CVrisk)
ages=c(40:80)
nstates=c("Health", "Ht","HyperLip","Dm","Cad","death","Ht&HyperLip","HyperLip&Dm","Ht&Dm","Ht&HyperLip&Dm")


atrisk=test
agesint=c(40:79)
ascvd.ten.year=ascvd.30.year=emp.ten.year=emp.lifetime=ten.year.new=lifetime.new=pce.reverse.tenyear=ascvd.30.year.rc=matrix(NA,nrow=nrow(test),ncol=length(agesint))


for(i in 1:length(agesint)) {
  age = agesint[i]


    dat = data.frame(
      "id" = atrisk$identifier,
      "mysex" = as.factor(atrisk$sex),
      "myage" = rep(age, nrow(atrisk)),
      "mysbp" = atrisk$sbp,
      "mytreat" = ifelse(atrisk$antihtn == 1 &
                           atrisk$htn_age < age, 1, 0),
      "mysmoking" = atrisk$smoke,
      "mydiabetes" = ifelse(atrisk$Dm_0_Any == 2 &
                              atrisk$Dm_0_censor_age < age, 1, 0),
      "mytotalchol" = atrisk$choladj,
      "myhdl" = atrisk$hdladj,
      "Race" = atrisk$Race,
      "mystatnow" = ifelse(atrisk$statin == 1 & atrisk$statin_age < age, 1, 0)
    )
    #
    frs = fun.frs_30yn(
      dat,
      id = "id",
      sex = "mysex",
      age = "myage",
      sbp = "mysbp",
      treat = "mytreat",
      smoking = "mysmoking",
      diabetes = "mydiabetes",
      totalchol = "mytotalchol",
      hdl = "myhdl"
    )
    ascvd.30.year[, i] = frs$frs_orig
    ascvd.30.year.rc[, i] = frs$frs_recali
    pce.reverse.tenyear[, i] = 
      compute_CVrisk2(
        df = dat,
        scores = "as2",
        age = "myage",
        gender = "mysex",
        race = "Race",
        totchol = "mytotalchol",
        sbp = "mysbp",
        hdl = "myhdl",
        bp_med = "mytreat",
        diabetes = "mydiabetes",
        smoker = "mysmoking",
        lipid_med = "mystatnow"
      )$as2
      
      
    
}

saveRDS(pce.reverse.tenyear,"multistate2/output/pce.reverse.tenyear.rds")
#saveRDS(ascvd.30.year.rc,"multistate2/output/ascvd.30year.rc.rds")
saveRDS(ascvd.30.year,"multistate2/output/ascvd.30year.rds")
saveRDS(ascvd.30.year.rc,"~/multistate2/output/ascvd.30year.rcnew.rds")



## compare to 