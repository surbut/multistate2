## questions for min : diagnosis year distribution

mu=readRDS("~/Dropbox/agedistribution.rds")
m=readRDS("~/Dropbox/agedistribution_mgb.rds")
levels(mu$variable)=levels(m$variable)=c("Diabetes","Hyperlipidemia","CAD","HTN")
ggplot(mu,aes(age,value,fill=variable))+geom_bar(stat="identity")+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
theme_classic()+
  labs(y="Count",x="Age of Diagnosis",fill="Risk Factor",title="Case Counts in UKB")+
  scale_x_discrete(breaks=seq(20, 80, 10))


ggplot(m,aes(age,value,fill=variable))+
geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_classic()+labs(y="Count",x="Age of Diagnosis",fill="Risk Factor",title="Case Counts in MGB")+
  scale_x_discrete(breaks=seq(20, 80, 10))

