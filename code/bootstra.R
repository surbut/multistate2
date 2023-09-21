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
    stop=min(age+10,80)
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

saveRDS(pm,"../output/predictedrsiskboot.rds")
saveRDS(pm,"../output/predictedrsiskboot_withtreat.rds")


pm=readRDS("../output/predictedrsiskboot.rds")

library(ggplot2)
library(viridis)

# Assuming `array1` is predicted risk and `array2` is risk under treatment
diff_array <- pmuntreated - pm

# Calculate mean and SD across bootstraps
mean_diff <- apply(diff_array, c(2,3,4), mean)
sd_diff <- apply(diff_array, c(2,3,4), sd)

# Melt the array to a dataframe for ggplot2
df <- as.data.frame(as.table(mean_diff))

ggplot(df, aes(x = Var2, y = Var3, fill = Freq)) + 
  geom_tile() +
  scale_fill_viridis() + 
  theme_minimal() +
  labs(x = "Years", y = "Risk States", fill = "Mean Difference")


oldtest=test

## now do for a fixed matrix
intercept=1
cad.prs=c(-2,-1,0,1,2)
sex=c(0,1)
smoke=c(0,1)
antihtn_now=c(0,1)
statin_now=c(0,1)
B=50
atrisk=expand.grid(intercept,cad.prs,sex,smoke,antihtn_now,statin_now)


test=data.frame(atrisk)
rownames(test)=c(1:nrow(test))
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

predicted_risks_ten=predicted_risks_ten_treat=predicted_risks
ages=c(40:80)

time=proc.time()
B=50
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
    #stop = 80 ## here we use lifetime
    stop=min((age+10),80)
   age=agesint[i]
    start=age
    #stop=80
    ar=test
    stop=min(age+10,80)

    
    mso = compute_prediction_product_matrix(
      coefmat = coefs,
      atrisk = ar,
      agepredinterval = c(start:stop)
    )
    
    predicted_risks_ten[t,, i, j] = mso$PredictedIntervalrisk
    predicted_risks_ten_treat[t,, i, j] = mso$Hazard_treated
    print(paste0("completeboot",t,"age",age,"state",j))

}
}}



pm2=predicted_risks_ten[c((1:4),(6:50)),,,]
saveRDS(pm2,"../output/predictedrsiskboot_fixed_ten.rds")

pm2=predicted_risks_ten_treat[c((1:4),(6:50)),,,]
saveRDS(pm2,"../output/predictedrsiskboot_fixedtentreat.rds")


pm2t=predicted_risks_ten_treat[c((1:4),(6:50)),,,]

saveRDS(pm2,"../output/predictedrsiskboot_fixed.rds")

saveRDS(pm2,"../output/predictedrsiskboot_fixed_benefit.rds")

untreated=readRDS("../output/predictedrsiskboot_fixed.rds")
treated=pm2
# Assuming `array1` is predicted risk and `array2` is risk under treatment
diff_array <- untreated[,c(6:10),,] - treated[,c(6:10),,] 

# Calculate mean and SD across bootstraps
mean_diff <- apply(diff_array, c(2,3,4), mean)
sd_diff <- apply(diff_array, c(2,3,4), sd)

# Melt the array to a dataframe for ggplot2
df <- as.data.frame(as.table(mean_diff))

ggsave(ggplot(df, aes(x = Var2, y = Var3, fill = Freq)) + 
  geom_tile() +
  scale_fill_viridis() + 
  theme_minimal() +
  labs(x = "Years", y = "Risk States", fill = "Mean Difference"),file="../output/benefit.pdf",dpi=600)


# Calculate mean and SD across bootstraps
mean_diff <- apply(diff_array, c(2,3,4), mean)
sd_diff <- apply(diff_array, c(2,3,4), sd)

# Melt the array to a dataframe for ggplot2
df <- as.data.frame(as.table(mean_diff))

ggplot(df, aes(x = Var2, y = Var3, fill = Freq)) + 
  geom_tile() +
  scale_fill_viridis() + 
  theme_minimal() +
  labs(x = "Years", y = "Risk States", fill = "Mean Difference")

library(ggplot2)
library(ggridges)

df <- as.data.frame(as.table(apply(diff_array, c(3,4), mean)))


benefit=ggplot(df, aes(x = Freq, y = as.factor(Var2), fill = ..x..)) + 
geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
theme_ridges() + scale_fill_viridis_b()+
labs(x = "Mean Difference in Risk", y = "Years",fill="Absolute Risk Deifference")

saveRDS(benefit,file = "../output/benefit.rds")
s=statusarray(df_frame = data.table(test),ages = ages,nstates = nstates)

pm2=readRDS("../output/predictedrsiskboot_fixed.rds")

library(ggplot2)
library(tidyr)
library(dplyr)

# Function to convert the data of one person to a data frame
person_to_df <- function(pm, person_idx) {
  # Extract the person's data
  person_data <- pm[, person_idx, , ]
  
  # State names from the 4th dimension of pm
  state_names <- dimnames(pm)[[4]]
  age_names=dimnames(pm)[[3]]
  # Convert to a data frame
  df <- expand.grid(Bootstrap = 1:dim(pm)[1], Age = age_names, State = state_names)
  df$Risk <- c(person_data)
  
  return(df)
}

# Function to get the summary statistics for one person
get_summary <- function(person_df) {
  summary_df <- person_df %>%
    group_by(Age, State) %>%
    summarise(
      mean_risk = mean(Risk),
      #lower_ci = quantile(Risk, 0.025),
      #upper_ci = quantile(Risk, 0.975),
      lower_ci = mean(Risk)-1.96*sd(Risk)/sqrt(n()),
      upper_ci = mean(Risk)+1.96*sd(Risk)/sqrt(n()),
      .groups = "drop"
    )
  return(summary_df)
}

get_summary2 <- function(person_df) {
  summary_df <- person_df %>%
    group_by(Age, State) %>%
    summarise(
      mean_risk = mean(Risk),
      lower_ci = quantile(Risk, 0.025),
      upper_ci = quantile(Risk, 0.975),
      #lower_ci = mean(Risk)-1.96*sd(Risk)/sqrt(n()),
      #upper_ci = mean(Risk)+1.96*sd(Risk)/sqrt(n()),
      .groups = "drop"
    )
  return(summary_df)
}

# Plotting function for a given person
plot_risk_for_person <- function(pm, person_idx) {
  person_df <- person_to_df(pm, person_idx)
  summary_df <- get_summary(person_df)
  
  ggplot(summary_df, aes(x = Age, y = mean_risk, color = State, group = State)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = State), alpha = 0.2) +
    labs(title = paste("Predicted Risks for Person", person_idx), y = "Risk", x = "Age", color = "State") +
    theme_classic()
}


plot_risk2_for_person <- function(pm, person_idx) {
  person_df <- person_to_df(pm, person_idx)
  summary_df <- get_summary2(person_df)
  
  ggplot(summary_df, aes(x = Age, y = mean_risk, color = State, group = State)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = State), alpha = 0.2) +
    labs(title = paste("Predicted Risks for Person", person_idx), y = "Risk", x = "Age", color = "State") +
    theme_classic()
}

# Plot for the first person as an example


p1=plot_risk_for_person(pm, "5021941")
# Plot for the first person as an example


p2=plot_risk2_for_person(pm, "5021941")

p3=get_person_summary_plot(pm, s2, "5021941")


o=readRDS("../output/continuosums.rds")

frs=ascvd.30.year.rc[test$identifier%in%"5021941",]

df=data.frame(ages=ages[-41],score=frs/100)

ggplot(df,aes(ages,y=score))+geom_point()
ga=ggarrange(p1,p3[[1]],o,common.legend = T,ncol =3,nrow =1)
ggsave(ga,file="../output/predictedrisksfor5021941.pdf",width = 18,units = "in")



pmidx="5021917"
p1=plot_risk_for_person(pm,pmidx)
# Plot for the first person as an example


p2=plot_risk2_for_person(pm, pmidx)

p3=get_person_summary_plot(pm, s2, pmidx)

ga=ggarrange(p1,p3[[1]],common.legend = T,ncol =2,nrow = 1)
ggsave(ga,file=paste0("../output/predictedrisksfor",pmidx,".pdf"))

