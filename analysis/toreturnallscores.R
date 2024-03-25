# how to return the array that has the score for every person
## pm is the bootrapped risk under all models


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
source("~/multistate2/code/frs30_URBUT/fun.frs_30y.R")
library(CVrisk)
ages=c(40:80)
nstates=c("Health", "Ht","HyperLip","Dm","Cad","death","Ht&HyperLip","HyperLip&Dm","Ht&Dm","Ht&HyperLip&Dm")


B=50

## increase bootstrap

# for(i in 1:B){
#   
#   row=sample(nrow(train),size = nrow(train),replace = T)
#   modelfit = fitfunc2(
#     data.table(train[row,]),
#     ages = ages,
#     nstates = nstates,
#     mode = "binomial",
#     covariates = "cad.prs+f.31.0.0+smoke+antihtn_now+statin_now")
#   
#   saveRDS(object = modelfit,file = paste0("~/multistate2/output/bootstrap",i,".rds"))
# }


##

## return a matrix of coefficients over all ages for a given state to state transition

agesint=c(40:79)

predicted_risks = array(NA, dim = c(B, nrow(test),length(agesint), 8))
dimnames(predicted_risks) = list(1:B,
                                 as.numeric(as.character(test$identifier)),
                                 agesint,
                                 c(
                                   "Health",
                                   "Ht",
                                   "HyperLip",
                                   "Dm",
                                   "Ht&HyperLip",
                                   "HyperLip&Dm",
                                   "Ht&Dm",
                                   "Ht&HyperLip&Dm"
                                 )
)

time=proc.time()
for(t in c((1:4),(6:B))){
  
  modelfit=readRDS(paste0("~/multistate2/output/bootstrap",t,".rds"))
  
  a = coefplotsmooth2(
    ages = ages,
    start = "Health",
    stop = "Cad",
    modelfit = modelfit,
    window_width = 20,
    span = 0.75,
    degree = 2
  )
  healthy_coefs = a$custom_smooth
  
  b = coefplotsmooth2(
    ages = ages,
    start = "Ht",
    stop = "Cad",
    modelfit = modelfit,
    window_width = 30,
    span = 0.75,
    degree = 2
  )
  ht_coefs = b$custom_smooth
  
  c = coefplotsmooth2(
    ages = ages,
    start = "HyperLip",
    stop = "Cad",
    modelfit = modelfit,
    window_width = 20,
    span = 0.75,
    degree = 2
  )
  hyperlip_coefs = c$custom_smooth
  
  d = coefplotsmooth2(
    ages = ages,
    start = "Dm",
    stop = "Cad",
    modelfit = modelfit,
    window_width = 20,
    span = 0.75,
    degree = 2
  )
  dm_coefs = d$custom_smooth
  
  
  e = coefplotsmooth2(
    ages = ages,
    start = "Ht&HyperLip",
    stop = "Cad",
    modelfit = modelfit,
    window_width = 20,
    span = 0.75,
    degree = 2
  )
  HH_coefs = e$custom_smooth
  
  f = coefplotsmooth2(
    ages = ages,
    start = "HyperLip&Dm",
    stop = "Cad",
    modelfit = modelfit,
    window_width = 20,
    span = 0.75,
    degree = 2
  )
  HD_coefs = f$custom_smooth
  
  g = coefplotsmooth2(
    ages = ages,
    start = "Ht&Dm",
    stop = "Cad",
    modelfit = modelfit,
    window_width = 20,
    span = 0.75,
    degree = 2
  )
  TD_coefs = g$custom_smooth
  
  h = coefplotsmooth2(
    ages = ages,
    start = "Ht&HyperLip&Dm",
    stop = "Cad",
    modelfit = modelfit,
    window_width = 20,
    span = 0.75,
    degree = 2
  )
  
  HHD_coefs = h$custom_smooth
  
  coeflist = list(
    healthy_coefs,
    ht_coefs,
    hyperlip_coefs,
    dm_coefs,
    HH_coefs,
    HD_coefs,
    TD_coefs,
    HHD_coefs
  )
  
  agesint = seq(40, 79, by = 1)
  
  
  for (j in 1:length(coeflist)) {
    coefs = coeflist[[j]]
    
    for (i in 1:length(agesint)) {
      age = agesint[i]
      start = age
      #stop = 80
      #stop=min(age+30,80)
      stop=min(age+1,80)
      atrisk = test
      ar = data.frame(
        "intercept" = 1,
        "cad.prs" = atrisk$cad.prs,
        "sex" = atrisk$f.31.0.0,
        "smoke" = atrisk$smoke,
        "antihtn_now" = ifelse(atrisk$antihtn == 1 &
                                 atrisk$htn_age < age, 1, 0),
        "statin_now" = ifelse(atrisk$statin == 1 &
                                atrisk$statin_age < age, 1, 0)
        
      )
      
      
      mso = compute_prediction_product_matrix(
        coefmat = coefs,
        atrisk = ar,
        agepredinterval = c(start:stop)
      )
      
      predicted_risks[t,, i, j] = mso$PredictedIntervalrisk
      #predicted_risks[t,, i, j] = mso$Hazard_treated
      print(paste0("completeboot",t,"age",age,"state",j))
      
    }
    
  }
}

print(proc.time()-time)

pm=predicted_risks[c((1:4),(6:50)),,,]

saveRDS(pm,"~/multistate2/output/predictedrsiskboot_1year.rds")
#saveRDS(pm,"~/multistate2/output/predictedrsiskboot_30year.rds")

#saveRDS(pm,"../output/predictedrsiskboot.rds")

###


#pm=readRDS("../output/predictedrsiskboot.rds")
s=statusarray(df_frame = data.table(test),ages = ages,nstates = nstates)
s2=s[c(1:40),,c(1,2,3,4,7,8,9,10)]

get_person_cross_threshold <- function(pm, s2, person_idx) {
  # Convert the risk array for the given person to a data frame
  risk_df <- as.data.frame(as.table(pm[, person_idx, , ]))
  colnames(risk_df) <- c("Bootstrap", "Age", "State", "Risk")
  
  # Convert the state indicator array for the given person to a data frame
  state_df <- as.data.frame(as.table(s2[, person_idx, ]))
  colnames(state_df) <- c("Age", "State", "Indicator")
  
  # Extract the state info for the given person
  current_state_df <- state_df %>%
    filter(Indicator == 1) %>%
    select(-Indicator)
  
  # Merge the data to get risks based on the state the person is in
  merged_df <- risk_df %>%
    inner_join(current_state_df, by = c("Age", "State"))
  
  # Compute mean and SE for the person across ages
  summary_df <- merged_df %>%
    group_by(Age) %>%
    summarise(
      mean_risk = mean(Risk),
      se_risk = sd(Risk) / sqrt(n()),
      #lower_95 = mean_risk - 1.96 * se_risk,  # For 95% CI
      #upper_95 = mean_risk + 1.96 * se_risk,  # For 95% CI
      lower_95 = quantile(Risk,0.025),  # For 95% CI
      upper_95 = quantile(Risk,0.975),  # For 95% CI
      num_exceeded=sum(Risk>thresh),
      
      .groups = "drop"
    )
  
  return(list(mean_risk=summary_df$mean_risk,se=summary_df$se_risk,num_exceeded=summary_df$num_exceeded))}

states=sapply(c(1:nrow(test)),function(x){get_person_cross_threshold(pm,s2,as.character(test$identifier[x]))$mean_risk})

#saveRDS(states,"~/multistate2/output/state_occupancy_risk_oneyear.rds")
#saveRDS(states,"~/multistate2/output/state_occupancy_risk_30.rds")

saveRDS(states,"~/multistate2/output/state_occupancy_risk.rds")

states=sapply(c(1:nrow(test)),function(x){get_person_cross_threshold(pm,s2,as.character(test$identifier[x]))$mean_risk})
