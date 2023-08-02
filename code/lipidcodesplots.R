## lipid org

cad=readRDS("~/Dropbox/pheno_dir/output/Cad_list_table.rds")
sta_date=readRDS("~/Dropbox/sta_gp.rds")
ldl_gp=readRDS("~/Dropbox/")
events=data.frame(cad$all_event_dt.Include_in_cases)
c=events%>%group_by(identifier)%>%summarise(cad_diag=min(eventdate))
s=sta_date%>%group_by(eid)%>%summarise(sta_start=min(from),age=min(mage))
ldl=readRDS("~/Dropbox/ldl_date.rds")

m=merge(c,s,by.x = "identifier",by.y="eid")

lipidplotfunc(intersect(s$eid,intersect(c$identifier,ldl$eid))[1])
lipidplotfunc=function(id){
  min=min(min(sta_date$from[sta_date$eid==id][1],min(ldl[ldl$eid==id,event_dt])),min(cad$all_event_dt.Include_in_cases[identifier%in%id,eventdate]))-3
  max=max(max(ldl[ldl$eid==id,event_dt],max(sta_date$from[sta_date$eid==id])),min(cad$all_event_dt.Include_in_cases[identifier%in%id,eventdate]))+3
ggplot(ldl[ldl$eid==id,],
         aes(x=event_dt,y=value1,color=value1))+
    geom_point()+geom_line()+
    geom_vline(xintercept=sta_date$from[sta_date$eid==id][1],color="green")+
               xlim(min,max)+
                          geom_vline(xintercept=min(cad$all_event_dt.Include_in_cases[identifier%in%id,eventdate]),
                                      color="red")+
               labs(title="Statin Intro and disease")+labs(y="LDL Cholesterol")+labs(x="Time of Measurement")}


lipidplotfunc(1002769)