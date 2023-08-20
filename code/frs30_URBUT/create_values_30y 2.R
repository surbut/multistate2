#||||||||||||||||||||||||||||||||
#
# reverse engineer the excel sheet calcualtor
#
#||||||||||||||||||||||||||||||||

# this script shall provide the coeeficients, means and prop from the sample and col2,col3 column G etc from the excel sheet
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(survival)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#::::::::::::::::::::::::::::
# get data ----
#::::::::::::::::::::::::::::

# this is the original data set
# see .docx for structure
dat <- readRDS("dat30.RDS")

# prep : logaritmize
dat$ln_age <- log(dat$myage)
dat$ln_sysbp <- log(dat$mysysbp)
dat$ln_totalchol <- log(dat$mytotalchol)
dat$ln_hdl <- log(dat$myhdl)
dat$ln_bmi <- log(dat$mybmi)


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

#...............................
# >>> refitted ----
#...............................

# create survival object
dat$SurvObj_event <- with(dat, Surv(mytime, myevent))
dat$SurvObj_compete <- with(dat, Surv(mytime, mycompete))

# fit cox reg
formul_event <- as.formula(paste("SurvObj_event~", paste(mypreds2, collapse = "+")))
formul_compete <- as.formula(paste("SurvObj_compete~", paste(mypreds2, collapse = "+")))

evmod <- coxph(formul_event, data = dat)
compmod <- coxph(formul_compete, data = dat)

# coefficents
coeffs_cvd_refit <- coef(evmod) 
coeffs_compete_refit <- coef(compmod) 

## >>>>this is for paper tabel----
sumevmod <- summary(evmod)
sumcompmod <- summary(compmod)

coffies <- sumevmod$conf.int[,-2 ]
evcoffies <- as.data.frame(apply(coffies, 2, function(x){as.numeric(sprintf("%.2f", x))}))
rownames(evcoffies) <- rownames(coffies)
evcoffies

coffies <- sumcompmod$conf.int[,-2 ]
compcoffies <- as.data.frame(apply(coffies, 2, function(x){as.numeric(sprintf("%.2f", x))}))
rownames(compcoffies) <- rownames(coffies)
compcoffies

#::::::::::::::::::::::::::::
#  means and props ----
#::::::::::::::::::::::::::::
source("fun.mean_and_prop.r")
tmp <- fun.mean_and_prop(dat, mypreds2)
mymenas <- unlist(tmp)[c("mysex.male", "ln_age", "ln_sysbp", "mytreat.TRUE", "mysmoking.TRUE", "mydiabetes.TRUE", "ln_totalchol", "ln_hdl" )]


#::::::::::::::::::::::::::::
#  mxbeta and mdxbeta ----
#::::::::::::::::::::::::::::
# taken form excel sheet
mxbeta_orig  <- 21.29326612
mdxbeta_orig  <- 20.12840698

mxbeta_recali <- sum(mymenas*coeffs_cvd_orig)
# [1] 21.75659
mdxbeta_recali <- sum(mymenas*coeffs_compete_orig)
# [1] 20.57296


mxbeta_refit <- sum(mymenas*coeffs_cvd_refit)
# [1] 26.22341
mdxbeta_refit <- sum(mymenas*coeffs_compete_refit)
# [1] 19.46958

#::::::::::::::::::::::::::::
#  col2 and col3 ----
#::::::::::::::::::::::::::::

mybasehaz_ev <- basehaz(evmod, centered = TRUE) # deafult is centered
mybasehaz_comp <- basehaz(compmod, centered = TRUE) # deafult is centered

col2 <- c(1- mybasehaz_ev$hazard)
col3 <- c(1- mybasehaz_comp$hazard)
# both have dim 1414
# this is the combined number of events and censoring

write.csv(col2, file = "col2_recali.csv", row.names = FALSE)
write.csv(col3, file = "col3_recali.csv", row.names = FALSE)

# these files have been manually created by copying from the excelsheet
col2_orig <- scan(file = "col2_orig.csv")
col3_orig <- scan(file = "col3_orig.csv")

#::::::::::::::::::::::::::::
#  equation A3 == column G ----
#::::::::::::::::::::::::::::

# log (i-1)-log(i)
auxi<- c(col2[-1],NA)
columnG <- log(col2)-log(auxi)

columnG_orig <- scan(file = "columnG_orig.csv")
write.csv(columnG, file = "columnG_recali.csv", row.names = FALSE)


# cleanup so this canbe sourced
# rm(list=setdiff(ls(),
#                 c("coeffs_compete_orig", "coeffs_compete_refit", "coeffs_cvd_orig",
#   "coeffs_cvd_refit",  "col2", "col3","columnG",
#   "col2_orig", "col3_orig","columnG_orig",
# 
#    "mdxbeta_recali", "mdxbeta_refit", "mxbeta_recali", "mdxbeta_orig", "mxbeta_orig",
#   "mxbeta_refit")))

rm(list=c("auxi", "coffies", "compcoffies", "compmod", "dat", "evcoffies", 
     "evmod", "formul_compete", "formul_event", "fun.mean_and_prop", 
      "mybasehaz_comp", "mybasehaz_ev", 
     "mymenas", "mypreds2", "sumcompmod", "sumevmod", "tmp"))
