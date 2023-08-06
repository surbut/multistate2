## coefficient analysis
c=coefinput(start = "Ht",stop = "Cad",ages = c(40:80),modelfit = abinom)
m=melt(c,id.vars="age")
ggplotly(ggplot(m,aes(age,exp(value),col=variable,fill=variable))+stat_smooth()+
  geom_hline(yintercept = 1,col=black))

coefplotsmoothashr=
function(ages,start,stop,modelfit){
  require(ggplot2)
  agenames=as.character(c(ages))
  s=sapply(agenames,function(x){modelfit$model_list[[x]][[stop]][[start]][,"Estimate"]})
  s=data.frame(s)
  e=sapply(agenames,function(x){modelfit$model_list[[x]][[stop]][[start]][,"Std. Error"]})
  
  colnames(s)=colnames(e)=agenames
  s=t(s)
  e=t(e)
  melt=melt(s)
  g=ggplot(melt,aes(Var1,value,col=Var2))+stat_smooth()+geom_point()
  m=ggplot_build(g)$data[[1]]
  return(list("mat"=m,"plot"=g))}