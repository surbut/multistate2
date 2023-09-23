library("survMisc")
df=data.table(na.omit(readRDS("~/Library/CloudStorage//Dropbox-Personal//phenotypes/df_ani_new.rds")))

dat=data.frame(df %>%group_by(round(phenos.enrollment,0)) %>%summarise(length(phenos.enrollment)))[c(3:32),]

hazards=matrix(data=NA,nrow=nrow(dat),ncol=6)
#hazards=matrix(data=NA,nrow=length(levels(c)),ncol=3)
df$round_age=round(df$phenos.enrollment,0)
df$ldl.quant=scale(df$ldladj)
for(i in seq(1:nrow(dat))){
  #for(i in 1:length(levels(c))){
  age=dat[i,1]
  #age=levels(c)[i]
  d=df[df$round_age<(age+2)&df$round_age>(age-2)&df$sex=="female",]
  #d=df[df$c==age,]
  sfit=coxph(Surv(phenos.CAD_censor_age,phenos.has_CAD)~prs_quant,data=d)
  hazards[i,1]=exp(coef(sfit)[1])
  reg=glm(phenos.has_CAD~prs_quant,data=d,family="binomial")
  hazards[i,2]= exp(coef(reg)[[2]])
  hazards[i,3]=with(summary(reg), 1 - deviance/null.deviance)
  #hazards[i,4]=rsq(sfit)$mev
  
  d=df[df$round_age<(age+2)&df$round_age>(age-2)&df$sex=="male",]
  sfit=coxph(Surv(phenos.CAD_censor_age,phenos.has_CAD)~prs_quant,data=d)
  hazards[i,4]=exp(coef(sfit)[1])
  reg=glm(phenos.has_CAD~prs_quant,data=d,family="binomial")
  hazards[i,5]= exp(coef(reg))[[2]]
  hazards[i,6]=with(summary(reg), 1 - deviance/null.deviance)
  
  #reg=glm(phenos.has_CAD~ldladj,data=d,family="binomial")
  #hazards[i,8]=rsq(sfit)$mev
  
  
  #hazards[i,3]= exp(coef(coxph(Surv(d$phenos.CAD_censor_age,d$old_cad)~prs_quant,data=d))[1])
  
  #hazards[i,6]= exp(coef(glm(old_cad~prs_quant,data=d,family="binomial"))[[2]])
}


hazards=data.frame(hazards[-1,])
colnames(hazards)=c("prs.HR.female","prs.OR.female","prs.r2.female","prs.HR.male","prs.OR.male","prs.r2.male")

hrmat=hazards[,grep("HR",colnames(hazards))];hrmat$age=dat$round.phenos.enrollment..0.[-1]


ormat=hazards[,grep("OR",colnames(hazards))];ormat$age=dat$round.phenos.enrollment..0.[-1]
r2mat=hazards[,grep("r2",colnames(hazards))];r2mat$age=dat$round.phenos.enrollment..0.[-1]



haz=melt(hrmat,id.vars = "age")
ors=melt(ormat,id.vars = "age")
r2s=melt(r2mat,id.vars = "age")

hr=ggplot(haz, aes(age,value,colour=variable,fill=variable))+
stat_smooth(method ="loess")+labs(y="HR",colour="Coefficient")+theme_classic()+guides(color=guide_legend(override.aes=list(fill=NA)))

or=ggplot(ors, aes(age,value,colour=variable,fill=variable))+ylim(1,2.2)+
stat_smooth(method = "loess")+labs(y="OR",colour="Coefficient")+theme_classic()+guides(color=guide_legend(override.aes=list(fill=NA)))

r2=ggplot(r2s, aes(age,value,colour=variable,fill=variable))+
stat_smooth(method = "loess")+labs(y="R2",colour="Coefficient")+theme_classic()+guides(color=guide_legend(override.aes=list(fill=NA)))

ggarrange(hr,r2,or,nrow=1,common.legend = T)

```

Now let's do it with survsplit:

