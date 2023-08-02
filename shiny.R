library(shiny)
library(ggplot2)
ui <- fluidPage(
  tags$h3("Scatter plot generator"),
  selectInput(inputId = "x", label = "X Axis", choices = names(nstates), selected = "Hm=fread("~/Dropbox (Partners HealthCare)/MGBB (1)/Phenotypes/Manual Curated/MGBB_Phenos_2CODE 2022-06-21.csv")
m=fread("~/Dropbox (Partners HealthCare)/MGBB (1)/Phenotypes/Manual Curated/MGBB_Phenos_2CODE 2022-06-21.csv")
m=fread("~/Dropbox (Partners HealthCare)/MGBB (1)/Phenotypes/Manual Curated/MGBB_Phenos_2CODE 2022-06-21.csv")
m=fread("~/Dropbox (Partners HealthCare)/MGBB (1)/Phenotypes/Manual Curated/MGBB_Phenos_2CODE 2022-06-21.csv")
mu=melt(mu,id.vars="age")
head(mu)
ggplot(mu,aes(age,value,fill=variable))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
),
  selectInput(inputId = "y", label = "Y Axis", choices = names(mtcars), selected = "Age"),
  plotOutput(outputId = "scatterPlot")
)
server <- function(input, output, session) {
  data <- reactive({mtcars})
  
  output$scatterPlot <- renderPlot({
    ggplot(data = data(), aes_string(x = input$x, y = input$y)) + 
      geom_point(aes(size = qsec, color = factor(cyl))) + 
      scale_color_manual(values = c("#3C6E71", "#70AE6E", "#BEEE62")) +
      theme_classic() + 
      theme(legend.position = "none")
    
    m=exp(cc[,"Cad",x]*1+cint[,"Cad",x])/(1+exp(cc[,"Cad",x]+cint[,"Cad",x]))
    m2=exp(cc[,"Cad",x]*2+cint[,"Cad",x])/(1+exp(cc[,"Cad",x]*2+cint[,"Cad",x]))
    m0=exp(cc[,"Cad",x]*0+cint[,"Cad",x])/(1+exp(cc[,"Cad",x]*0+cint[,"Cad",x]))
    mm=exp(cc[,"Cad",x]*-1+cint[,"Cad",x])/(1+exp(cc[,"Cad",x]*-1+cint[,"Cad",x]))
    mt=exp(cc[,"Cad",x]*-2+cint[,"Cad",x])/(1+exp(cc[,"Cad",x]*-2+cint[,"Cad",x]))
    n=event[,"Cad",x]
    
    #mf=data.frame(m,m2,m0,mm,mt)
    
    mf=data.frame(m,m0,mm)
    names(mf)=c("1SD","median","-1SD")
    mf[which(n<10),]="NA"
    mf$age=as.numeric(rownames(mf))
    #mf$m=as.numeric(m)
    mf=na.omit(mf)
    mfr=melt(mf,id.vars="age")
    mfr$value=as.numeric(mfr$value)
    p=ggplot(mfr,aes(age,round((1000*value),0),col=variable))+geom_point(size=0)+geom_smooth(aes(fill=variable),span=1,se=F)+labs(y=paste0("AR",dimnames(event)[[3]][x],"toCAD"))+theme_classic(base_size = 10)
  })
}
shinyApp(ui = ui, server = server)


x=ind[i]
m=exp(cc[,"Cad",x]*1+cint[,"Cad",x])/(1+exp(cc[,"Cad",x]+cint[,"Cad",x]))
m2=exp(cc[,"Cad",x]*2+cint[,"Cad",x])/(1+exp(cc[,"Cad",x]*2+cint[,"Cad",x]))
m0=exp(cc[,"Cad",x]*0+cint[,"Cad",x])/(1+exp(cc[,"Cad",x]*0+cint[,"Cad",x]))
mm=exp(cc[,"Cad",x]*-1+cint[,"Cad",x])/(1+exp(cc[,"Cad",x]*-1+cint[,"Cad",x]))
mt=exp(cc[,"Cad",x]*-2+cint[,"Cad",x])/(1+exp(cc[,"Cad",x]*-2+cint[,"Cad",x]))
n=event[,"Cad",x]

#mf=data.frame(m,m2,m0,mm,mt)

mf=data.frame(m,m0,mm)
names(mf)=c("1SD","median","-1SD")
mf[which(n<10),]="NA"
mf$age=as.numeric(rownames(mf))
#mf$m=as.numeric(m)
mf=na.omit(mf)
mfr=melt(mf,id.vars="age")
mfr$value=as.numeric(mfr$value)
p=ggplot(mfr,aes(age,round((1000*value),0),col=variable))+geom_point(size=0)+geom_smooth(aes(fill=variable),span=1,se=F)+labs(y=paste0("AR",dimnames(event)[[3]][x],"toCAD"))+theme_classic(base_size = 10)