#||||||||||||||||||||||||||||||||
#
# reverse engineer the excel sheet calcualtor
#
#||||||||||||||||||||||||||||||||

# this script shall provide the coeeficients, means and prop from the sample and col2,col3 column G etc from the excel sheet

#::::::::::::::::::::::::::::
# get data ----
#::::::::::::::::::::::::::::
load("~/Library/CloudStorage/Dropbox-Personal///pheno_dir/output/merged_pheno_censor_final_withdrugs_smoke.rds")
dfascvd=readRDS("~/multistate2//output/dfascvd_newbp.rds")

atrisk=merge(dfh,dfascvd[,c("sample_id","sex","age","sbp","hdladj","choladj","Race")],by.x="identifier",by.y="sample_id")
# dat=data.frame("id"=atrisk$identifier,"mysex"=as.factor(atrisk$sex),
#                "myage"=rep(age,nrow(atrisk)),"mysbp"=atrisk$sbp,
#                "mytreat"=ifelse(atrisk$antihtn==1&atrisk$htn_age<age,1,0),
#                "mysmoking"=atrisk$smoke,"mydiabetes"=ifelse(atrisk$Dm_0_Any==2&atrisk$Dm_0_censor_age<age,1,0),
#                "mytotalchol"=atrisk$choladj,"myhdl"=atrisk$hdladj,"Race"=atrisk$Race,"statnow"=ifelse(atrisk$statin==1&atrisk$statin_age<age,1,0)
#                )

dat = data.frame(
  "id" = atrisk$identifier,
  "mysex" = as.factor(atrisk$sex),
  "myage" = atrisk$age,
  "mysbp" = atrisk$sbp,
  "mytreat" = ifelse(atrisk$antihtn == 1 &
                       atrisk$htn_age < atrisk$age, 1, 0),
  "mysmoking" = atrisk$smoke,
  "mydiabetes" = ifelse(atrisk$Dm_0_Any == 2 &
                          atrisk$Dm_0_censor_age < atrisk$age, 1, 0),
  "mytotalchol" = atrisk$choladj,
  "myhdl" = atrisk$hdladj,
  "Race" = atrisk$Race,
  "statnow" = ifelse(atrisk$statin == 1 &
                       atrisk$statin_age < atrisk$age, 1, 0),
  "myevent" = ifelse(atrisk$Cad_0_Any == 2, 1, 0),
  "mytime" = atrisk$Cad_0_censor_age,
  "mycompete" = ifelse(atrisk$Death_censor_Any == 2,1,0))


# this is the original data set
# see .docx for structure
dat$ln_age <- log(dat$myage)
dat$ln_sysbp <- log(dat$mysbp)
dat$ln_totalchol <- log(dat$mytotalchol)
dat$ln_hdl <- log(dat$myhdl)
#dat$ln_bmi <- log(dat$mybmi)

# reverse factors for sex
dat$mysex <- relevel(dat$mysex, ref = "female")



#::::::::::::::::::::::::::::
#  coefficients ----
#::::::::::::::::::::::::::::

mypreds2 <- c("mysex",
              "ln_age",
              "ln_sysbp",
              "mytreat",
              "mysmoking",
              "mydiabetes",   
              "ln_totalchol", 
              "ln_hdl")

#::::::::::::::::::::::::::::
#  coefficients ----
#::::::::::::::::::::::::::::

#...............................
# >>> orig ----
#...............................
# taken form excel sheet

#...............................
# >>> orig ----
#...............................
# taken form excel sheet

# CVD
coeffs_cvd_orig <- c("mysexmale" = 0.34362,
                     "ln_age" = 2.63588,
                     "ln_sysbp" = 1.8803,
                     "mytreatTRUE" = 0.5232,
                     "mysmokingTRUE" = 0.59397, 
                     "mydiabetesTRUE" = 0.68602,
                     "ln_totalchol" = 1.12673,
                     "ln_hdl" = -0.90941)
# competing events
coeffs_compete_orig <- c("mysexmale" = 0.48123,
                         "ln_age" = 3.39222,
                         "ln_sysbp" = 1.39862,
                         "mytreatTRUE" = 0.19035,
                         "mysmokingTRUE" = 0.99858, 
                         "mydiabetesTRUE" = 0.49756,
                         "ln_totalchol" = -0.00439,
                         "ln_hdl" = 0.16081)


## refitted
# 
# create survival object
dat$SurvObj_event <- with(dat, Surv(mytime, myevent))
dat$SurvObj_compete <- with(dat, Surv(mytime, mycompete))

# # fit cox reg
formul_event <- as.formula(paste("SurvObj_event~", paste(mypreds2, collapse = "+")))
formul_compete <- as.formula(paste("SurvObj_compete~", paste(mypreds2, collapse = "+")))
# 
evmod <- coxph(formul_event, data = dat)
compmod <- coxph(formul_compete, data = dat)


#  means and props ----
#::::::::::::::::::::::::::::
source("~/multistate2//code/frs30_URBUT/fun.mean_and_prop.r")
tmp <- fun.mean_and_prop(dat, mypreds2)
mymenas <- unlist(tmp)[c("mysex.male", "ln_age", "ln_sysbp", "mytreat", "mysmoking", "mydiabetes", "ln_totalchol", "ln_hdl" )]


#::::::::::::::::::::::::::::
#  mxbeta and mdxbeta ----
#::::::::::::::::::::::::::::
# taken form excel sheet
mxbeta_orig  <- 21.29326612
mdxbeta_orig  <- 20.12840698

(mxbeta_recali <- sum(mymenas*coeffs_cvd_orig))
# [1]  22.31084 22.74104 22.66713
(mdxbeta_recali <- sum(mymenas*coeffs_compete_orig))
# [1] ] 21.13211 21.58129 21.54294

# 
# 
mybasehaz_ev <- basehaz(evmod, centered = TRUE) # deafult is centered
mybasehaz_comp <- basehaz(compmod, centered = TRUE) # deafult is centered
# 
col2 <- c(1- mybasehaz_ev$hazard)
col3 <- c(1- mybasehaz_comp$hazard)
# # both have dim 1414
# # this is the combined number of events and censoring
# 
write.csv(col2, file = "col2_recali_new.csv", row.names = FALSE,col.names = FALSE)
write.csv(col3, file = "col3_recali_new.csv", row.names = FALSE,col.names = FALSE)

#col2=scan(file = "~/multistate2/code/frs30_URBUT/col2_recali_new.csv")
#col3=scan(file = "~/multistate2/code/frs30_URBUT/col3_recali_new.csv")

# log (i-1)-log(i)
auxi<- c(col2[-1],NA)
columnG <- log(col2)-log(auxi)

columnG_orig <- scan(file = "~/multistate2/code/frs30_URBUT/columnG_orig.csv")
write.csv(columnG, file = "columnG_recali_new.csv", row.names = FALSE)


# these files have been manually created by copying from the excelsheet
col2_orig <- scan(file = "~/multistate2/code/frs30_URBUT/col2_orig.csv")
col3_orig <- scan(file = "~/multistate2/code/frs30_URBUT/col3_orig.csv")


#::::::::::::::::::::::::::::
#  equation A3 == column G ----
#::::::::::::::::::::::::::::


columnG_orig <- scan(file = "~/multistate2/code/frs30_URBUT/columnG_orig.csv")

# cleanup so this canbe sourced
# rm(list=setdiff(ls(),
#                 c("coeffs_compete_orig", "coeffs_compete_refit", "coeffs_cvd_orig",
#   "coeffs_cvd_refit",  "col2", "col3","columnG",
#   "col2_orig", "col3_orig","columnG_orig",
# 
#    "mdxbeta_recali", "mdxbeta_refit", "mxbeta_recali", "mdxbeta_orig", "mxbeta_orig",
#   "mxbeta_refit")))

