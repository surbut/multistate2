## create sample data

testdf=load("~/Library/CloudStorage/Dropbox-Personal/pheno_dir/output/merged_pheno_censor_final_withdrugs_smoke.rds")
testdf=dfh

set.seed(1000);
testdf$Cad_0_censor_age=testdf$Cad_0_censor_age+rnorm(length(testdf$Cad_0_Any),0,sd=1)
testdf$Ht_0_censor_age=testdf$Ht_0_censor_age+rnorm(length(testdf$Cad_0_Any),0,sd=1)
testdf$HyperLip_0_censor_age=testdf$HyperLip_0_censor_age+rnorm(length(testdf$Cad_0_Any),0,sd=1)
testdf$Dm_0_censor_age=testdf$Dm_0_censor_age+rnorm(length(testdf$Cad_0_Any),0,sd=1)
testdf$Death_Censor_Age=testdf$Death_Censor_Age+rnorm(length(testdf$Cad_0_Any),0,sd=1)
testdf$identifier=seq(1:nrow(testdf))
df_full=testdf[,c("identifier","cad.prs","f.31.0.0","Cad_0_Any","Dm_0_Any","Ht_0_Any","Death_censor_Any",
                            "HyperLip_0_Any","Cad_0_censor_age","Dm_0_censor_age","Ht_0_censor_age","Death_Censor_Age",
                            "HyperLip_0_censor_age","statin","antihtn","statin_age","htn_age","smoke")]

saveRDS(df_full,"~/MSGene/data/msgene_sampledf.rds")