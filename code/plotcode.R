# function to create plots
plot_diamonds <- function(i,plotind) {
  df=int.mat
  df$e <- exp(int.mat[,i])/(1+exp(int.mat[,i]))
  df$Age=ages
  df$Variable=nstates[i]
  df$col=rep(plotind,nrow(df))
  ggplot(df, aes(Age,e,color=col)) +
    #geom_point(color=custom.col[plotind],size=0) +
    stat_smooth(color=custom.col[plotind],method = "loess",fill=custom.col[plotind])+
    theme_classic(base_size = 15)+
    xlab("Age")+
    ylab("Abs to CAD")+ggtitle(nstates[i])
}

### function to create plots
plot_ors <- function(i,mat,plotind) {
  df=int.mat
  df$e <- exp(mat[,i])
  df$Age=ages
  df$Variable=nstates[i]
  df$col=rep(plotind,nrow(df))
  ggplot(df, aes(Age,e,color=col)) +
    #geom_point(color=custom.col[plotind],size=0) +
    theme_classic(base_size = 15)+
    stat_smooth(color=custom.col[plotind],method = "loess",fill=custom.col[plotind])+
    xlab("Age")+
    ylab("OR to CAD")+ggtitle(nstates[i])
}