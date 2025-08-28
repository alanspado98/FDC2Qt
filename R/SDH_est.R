#SDH estimation (Maione et. al, 2003)
library(pals)
library(lubridate)
library(lmomRFA)
library(ppcc)
library(dplyr)
library(png)
library(pracma)

durations<-c(1:48)

#Check if we have hourly time step
file_names<-dir()
Qd_ds<-list()
rd_ds<-list()
d_threshold<-as.numeric(read.csv("Durata_soglia_Downscaling.csv",sep=";",dec=","))

for(fn in file_names)
{
  Q_hh<-read.csv(fn,sep=";",dec=",")
  Q_hh<-Q_hh[-which(is.na(Q_hh$VALORE)),]
  
  fn<-sub("\\..*","",sub("^[^_]*_","",fn))
  if(length(which(c("RenoPracchia","MarecchiaSS16","MontoneCastrocaro","RenoVergato")%in%fn))>0) Q_hh<-Q_hh[which(as.numeric(substring(Q_hh$DATA,first=1,last=4))%in%QMG$Year[which(QMG$Name%in%fn)]),] 
  
  #Trasformazione della serie di portata oraria in serie di portata giornaliera
  Q_gg <- Q_hh %>% group_by(DATA) %>% dplyr::summarize(VALORE = mean(VALORE, na.rm=TRUE)) %>% as.data.frame()

  #Curva di durata (portate medie giornaliere in ordine decrescente+probabilità di superamento Weibull)
  Q=unique(sort(Q_gg$VALORE,decreasing = TRUE))
  porfdc <- cbind.data.frame(p=c((1:length(Q))/(length(Q)+1)),Q=Q)
  
  Q_hh[,1]<-ymd_hms(paste(Q_hh[,1],Q_hh[,2],sep=" "),tz="GMT")
  colnames(Q_hh)[1]<-"DATAORA"
  Q_hh<-Q_hh[,-2]
  
  Q_gg[,1]<-as.Date(Q_gg[,1])
  
  Qgg_d <- Q_gg
  Qgg_d$d<-rep(NA,length(Qgg_d[,1]))
    
  for(q in porfdc$Q) Qgg_d$d[Qgg_d$VALORE%in%q] <- porfdc$p[porfdc$Q==q] #daily streamflow durations
  
  index<-c()
  for(i in 2:(length(Qgg_d$DATA)-1))
  {
    if(Qgg_d$DATA[i]-Qgg_d$DATA[i-1]==1 
       && Qgg_d$DATA[i+1]-Qgg_d$DATA[i]==1)
    {
     if(Qgg_d$VALORE[i]>Qgg_d$VALORE[i-1] && 
        Qgg_d$VALORE[i]>Qgg_d$VALORE[i+1] && 
        Qgg_d$d[i]<=d_threshold) index<-append(index,i)
    }
  }
  index<-unique(index)
  
  Qgg_d<-Qgg_d[index,]

  Q_d <- data.frame(matrix(NA,nrow=length(index),ncol=1+length(durations),dimnames = list(x<-NULL,y<-c("date",paste("value",as.character(durations),sep="_"))))) 
  rd_d <- data.frame(matrix(NA,nrow=length(index),ncol=1+length(durations),dimnames = list(x<-NULL,y<-c("date",paste("value",as.character(durations),sep="_"))))) 
  Q_d$date <- Qgg_d$DATA
  rd_d$date <- Qgg_d$DATA
  
  for(j in durations)
  {
    Media_Mobile <- function(x, n = j){stats::filter(x,rep(1/n,n),sides = 1)}
    Res_Media_Mobile <- cbind.data.frame(DATAORA=Q_hh$DATAORA,VALORE=Media_Mobile(Q_hh$VALORE))
    Res_Media_Mobile[1:(length(Res_Media_Mobile[,1])-j+1),2] <- Res_Media_Mobile[j:length(Res_Media_Mobile[,1]),2]
    Res_Media_Mobile<-Res_Media_Mobile[1:(length(Res_Media_Mobile[,1])-j+1),]
    
    for(d in 1:length(Qgg_d$DATA))
    { 
      #t_start <- Res_Media_Mobile$DATAORA[year(Res_Media_Mobile$DATAORA)%in%y][which.max(Res_Media_Mobile$VALORE[year(Res_Media_Mobile$DATAORA)%in%y])-j+1][1]
      #t_peak  <- Q_hh$DATAORA[year(Q_hh$DATAORA)%in%y][which.max(Q_hh$VALORE[year(Q_hh$DATAORA)%in%y])]

      Index_MM <- which(substring(Res_Media_Mobile$DATAORA,1,10)%in%c(Qgg_d$DATA[d]-1,Qgg_d$DATA[d],Qgg_d$DATA[d]+1))
      Index_HH <- which(substring(Q_hh$DATAORA,1,10)%in%c(Qgg_d$DATA[d]-1,Qgg_d$DATA[d],Qgg_d$DATA[d]+1))
      t_start <- Res_Media_Mobile$DATAORA[Index_MM][which.max(Res_Media_Mobile$VALORE[Index_MM])]
      t_peak  <- Q_hh$DATAORA[Index_HH][which.max(Q_hh$VALORE[Index_HH])]
    
      #Salvataggio dei risultati
      #rd_d[rd_d$year%in%y,colnames(rd_d)%in%paste("value",as.character(j),sep="_")] <- as.numeric(ymd_hms(t_peak)-ymd_hms(t_start))
      #Q_d[rd_d$year%in%y,colnames(Q_d)%in%paste("value",as.character(j),sep="_")] <- max(Res_Media_Mobile$VALORE[year(Res_Media_Mobile$DATAORA)%in%y],na.rm=T)
      rd_d[rd_d$date%in%Qgg_d$DATA[d],colnames(rd_d)%in%paste("value",as.character(j),sep="_")] <- as.numeric(t_peak-t_start)/j
      Q_d[Q_d$date%in%Qgg_d$DATA[d],colnames(Q_d)%in%paste("value",as.character(j),sep="_")] <- max(Res_Media_Mobile$VALORE[Index_MM],na.rm=T)
    } 
  }
  Qd_ds<-append(Qd_ds,list(Q_d))
  rd_ds<-append(rd_ds,list(rd_d))
}

names(Qd_ds)<-sub("\\..*","",sub("^[^_]*_","",file_names))
names(rd_ds)<-sub("\\..*","",sub("^[^_]*_","",file_names))

Qd_ds<-Q_ds
rd_ds<-r_ds
load("Event analysis/Q_d_events.RData")

load("Event analysis/r_d_events.RData")
#save(Qd_ds,file="Event analysis/Q_d_events.RData")
for(n in names(rd_ds))
{
  r_d<-rd_ds[[n]]
  for(c in 1:24) r_d[,colnames(r_d)%in%paste("value",as.character(c),sep="_")]<- r_d[,colnames(r_d)%in%paste("value",as.character(c),sep="_")]/c
  rd_ds[[n]]<-r_d
}
# save(rd_ds,file="r_d_events.RData")
# save(Qd_ds,file="Q_d_events.RData")

# pal.bands(kelly()[c(19,18,16,13,12,10,9,7,5,4,3,2)],show.names=FALSE)
colori <- kelly()
colori <- colori[c(3:10,12,16,19,20,21)] #Si rimuovono colori che possono confondersi con quelli scelti, fino ad averne 14
colori[3]<-"red"
par(mar=c(4.2,4.2,1,10),xpd=F)
plot(NULL,xlab="duration (-)",ylab=expression("CV("*Q[D]*")"),xlim=c(1,24),ylim=c(0,2),xaxt="n")
grid(nx=NA,ny=NULL,lwd=1)
axis(1,at=c(1,seq(3,24,3)),labels=c(1,seq(3,24,3)))
abline(v=c(1,seq(3,24,3)),col="lightgray",lty="dotted",lwd=1)

for(n in names(Qd_ds))
{
  Q_d<-Qd_ds[[n]]
  CV<-apply(Q_d[,-1],2,sd,na.rm=T)/apply(Q_d[,-1],2,mean,na.rm=T)
  lines(x=as.numeric(sub("^[^_]*_","",names(CV))),y=as.numeric(CV),col=colors[names(Qd_ds)%in%n],lwd=2)
}
par(xpd=T)
legend(x=25,y=2,legend=names(Qd_ds),col=colors,lwd=2,lty=1,cex=.8)
#Lamone StCasale -> Massimo annuale 2014

## Reduction ratio curve 
set_cex=1.5
par(mar=c(5.1, 6.1, 2.1, 2.1),mgp=c(4,1,0))
plot(NULL,xlab="Duration (h)",ylab=expression(paste(r[D]," (-)"),""),
     xlim=c(0,24),ylim=c(0,.6),
     xaxt="n",
     yaxt="n",
     cex.lab=set_cex)
axis(1,at=seq(0,24,3),labels=seq(0,24,3),cex.axis=set_cex)
axis(2,at=c(0,seq(0,0.6,0.1)),labels=c(0,seq(0,0.6,0.1)),cex.axis=set_cex,las=2)
abline(h=seq(0,0.6,0.1),col="gray",lty="solid",lwd=1)#griglia orizzontale
abline(v=seq(0,24,3),col="gray",lty="solid",lwd=1)#griglia vericale


durations<-c(1:48)
rd_d_mean<-data.frame(matrix(NA,1,length(durations)+1,dimnames=list(x=NULL,y=c("name",paste("value",as.character(durations),sep="_")))))
for(n in names(rd_ds))
{
  rd_d<-rd_ds[[n]]
  y<-c()
  for(c in 2:length(rd_d[1,]))
  {
    rd_d_c<-rd_d[,c]
    rd_d_c<-rd_d_c[which(rd_d_c>0&rd_d_c<1)]
    y<-append(y,mean(rd_d_c))
  }
  y<-c(n,y)
  names(y)<-c("name",paste("value",as.character(durations),sep="_"))
  rd_d_mean<-rbind(rd_d_mean,y)
}
rd_d_mean<-rd_d_mean[-1,]
rd_d_mean[,-1]<-apply(rd_d_mean[,-1],2,as.numeric)

rd_d_mean<-rd_d_mean[,1:25]
durations<-durations[1:24]

rd_d_mean<-rd_d_mean[order(factor(rd_d_mean[,1],levels=Tc$Name)),]
#durations<-durations[-which(durations%in%2)]

durations<-1:46
durations<-c(.001,durations)
for(n in names(rd_ds)) lines(x=durations,y=as.numeric(rd_d_mean[names(rd_ds)%in%n,-c(1,2)]),col=df_colori$values[which(df_colori$names%in%n)],lwd=3)

#par(xpd=T)
#legend(x=25,y=.6,legend=names(rd_ds),col=colors,lwd=2,lty=1,cex=.8)

#Plot for each station
img <- readPNG("Event analysis/CV_d.png")
for(n in names(rd_ds))
{
  n<-names(rd_ds)[1]
  rd_d<-rd_ds[[n]]
  #png(paste0("Event analysis/",n,".png"),height=dim(img)[1],width=dim(img)[2])

  par(mar=c(4.2,4.2,3,11),xpd=F)
  plot(NULL,xlab="duration (-)",ylab=expression(r[D]),xlim=c(1,24),ylim=c(0,1),xaxt="n",
       main=n)
  
  axis(1,at=c(1,seq(3,24,3)),labels=c(1,seq(3,24,3)))
  abline(v=c(1,seq(3,24,3)),col="lightgray",lty="dotted",lwd=1)
  index<-c()
  for(c in 2:length(rd_d[1,]))
  {
    if(length(which(rd_d[,c]<0))>0) index<-append(index,which(rd_d[,c]<0))
    if(length(which(rd_d[,c]>1))>0) index<-append(index,which(rd_d[,c]>1))
  }
  index<-unique(index)
  rd_d<-rd_d[-index,]

  for(i in 1:length(rd_d[,1])) lines(x=c(1,3:48),y=as.numeric(rd_d[i,-1:-2]),col=t_col("blue",70),lwd=2,lty=2,cex.axis=1.5,cex.lab=1.5,cex.main=2)
  lines(x=c(1,3:48),y=apply(rd_d[,-1:-2],2,mean),col="red",lwd=3)
  par(xpd=T)
  legend(x=25,y=1,legend=c("events","mean"),col=c(t_col("blue",70),"red"),lwd=2,lty=c(2,1),cex=1.5)
  dev.off()
}

#Reduction ratio model: rather than eq.16 Maione et al. 2003, Tasker and Stedinger, 1989
exponential_model_2<-function(pars)
{
  lambda1<-pars[1]
  lambda2<-pars[2]
  est<-.5*exp(-lambda1*norm_durations/(1+lambda2*norm_durations))
  sumsq<-sum((est-obs)^2)
  return(sumsq)
}

linear_model<-function(y) stats::lm(y ~ x)
#quadratic_model<-function(y) stats::lm(y ~ poly(x,3,raw = TRUE))

#Regression for lambda1 and lambda2 parameter
setwd("..")
setwd("..")
Tc<-read.csv("Marecchia/04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/abc_sezioni_chiusura.csv",sep=";",dec=",")
Tc<-Tc[,c(2,6)]
rd_d_mean<-rd_d_mean[order(factor(rd_d_mean[,1],levels=Tc$Name)),]
#rd_d_mean<-rd_d_mean[,1:25]
#durations<-durations[1:23]

rd_d_pars<-data.frame(matrix(NA,12,3,dimnames=list(x=NULL,y=c("name","lambda1","lambda2"))))
rd_d_pars$name<-names(rd_ds)

for(n in names(rd_ds))
{
  norm_durations<-durations/Tc$Tc[Tc[,1]%in%n]
  obs<-as.numeric(rd_d_mean[rd_d_mean[,1]%in%n,-c(1,2)])
  rd_d_model<-optim(c(1,1),exponential_model_2)
  rd_d_pars[rd_d_pars[,1]%in%n,-1]<-rd_d_model$par
}

x<-Tc$Tc#[-which(Tc$Name%in%c("SamoggiaCalcara","BorelloBorello"))]
#y<-rd_d_pars$lambda1[-which(rd_d_pars$name%in%c("SamoggiaCalcara","BorelloBorello"))]
#Initialization of the exponential model parameters dataframe for fitting lambda1=f(Tc), lambda2=f(Tc) relationships
lambda_Tc_pars <- data.frame(matrix(NA,nrow=2,ncol=3,dimnames = list(x=NULL,y=c("parameter","q","m"))))
lambda_Tc_pars$parameter <- c("lambda1","lambda2")
for(i in 1:2) lambda_Tc_pars[i,2:3]<-linear_model(rd_d_pars[,i+1])$coefficients

library(openxlsx)
ab_values<-read.xlsx("Tc_a_b_values.xlsx")


#Create a 2x1 layout: 2x2 plots + 1 column for legend
layout(matrix(c(1,2,3), nrow = 1, byrow = TRUE),
       widths = c(1, 1, 0.4))
par(mar = c(6, 4, 4, 1),xpd=F)
x_poly<-seq(min(x),max(x),length.out=1000)
y_poly<-lambda_Tc_pars$m*x_poly+lambda_Tc_pars$q

plot(ab_values$Tc_from_b,ab_values$b,xlab="time of concentration (h)",ylab="a (-)",col=colors,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5)
grid()

plot(Tc$Tc,rd_d_pars$lambda1,xlab="time of concentration (h)",ylab="a (-)",col=colors,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5)
grid()
lines(x_poly,y_poly,col="red",lwd=2)
x_poly<-seq(min(df$x),max(df$x),length.out=1000)
y_poly<-predict(b_model,newdata = data.frame(x = x_poly))
plot(Tc$Tc,rd_d_pars$lambda2,xlab="time of concentration (h)",ylab="b (-)",col=colors,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5)
grid()
lines(x_poly,y_poly,col="red",lwd=2)
# Empty plot for legend
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",legend=names(rd_ds),col=colors,pch=16,cex=1.2)
mtext(expression(r[D] == 1/2*exp(- (a * D/Tc) / (1 + b * D/Tc))),
      outer = TRUE, cex = 1.5, font = 2, line = -2.7)

##Reduction ratio curve
for(n in names(rd_ds))
{
  n<-names(Qd_ds)[1]
  Q_d<-Qd_ds[[n]]
  Q_d[,-1]<-Q_d[,-1]/Q_d[,2]
  #png(paste0("Event analysis/εd_d_empirical/",n,".png"),height=dim(img)[1],width=dim(img)[2])
  
  par(mar=c(4.2,4.2,3,11),xpd=F)
  plot(NULL,xlab="duration (-)",ylab=expression(Q[D]),xlim=c(1,24),ylim=c(0,1),xaxt="n",
       main=n)
  grid(nx=NA,ny=NULL,lwd=1)
  axis(1,at=c(1,seq(3,24,3)),labels=c(1,seq(3,24,3)))
  abline(v=c(1,seq(3,24,3)),col="lightgray",lty="dotted",lwd=1)
  
  index<-c()
  for(c in 2:length(Q_d[1,]))
  {
    if(length(which(Q_d[,c]<0))>0) index<-append(index,which(Q_d[,c]<0))
    if(length(which(Q_d[,c]>1))>0) index<-append(index,which(Q_d[,c]>1))
  }
  index<-unique(index)
  Q_d<-Q_d[-index,]
  
  for(i in 1:length(Q_d[,1])) lines(x=2:49,y=as.numeric(Q_d[i,-1]/Q_d[i,2]),col="blue",lwd=2,lty=2,cex.axis=1.5,cex.lab=1.5,cex.main=2)
  lines(x=2:49,y=apply(Q_d[,-1],2,mean),col="red",lwd=3)
  par(xpd=T)
  legend(x=25,y=1,legend=c("events","mean"),col=c(t_col("blue",70),"red"),lwd=2,lty=c(2,1),cex=1.5)  
  dev.off()
}

set_cex=1.5
par(mar=c(5.1, 6.1, 2.1, 2.1),mgp=c(4,1,0))
plot(NULL,xlab="Duration (h)",ylab=expression(paste(epsilon[D]," (-)"),""),
     xlim=c(0,24),ylim=c(0,1),
     xaxt="n",
     yaxt="n",
     cex.lab=set_cex)
axis(1,at=seq(0,24,3),labels=seq(0,24,3),cex.axis=set_cex)
axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2),cex.axis=set_cex,las=2)
abline(h=seq(0,1,0.2),col="gray",lty="solid",lwd=1)#griglia orizzontale
abline(v=seq(0,24,3),col="gray",lty="solid",lwd=1)#griglia vericale

durations<-as.numeric(sub("^[^_]*_","",names(Q_d)))
durations[1]<-0.001
durations<-c(durations,47)
for(n in names(Qd_ds))
{
  Q_d<-Qd_ds[[n]][,-1]
  Q_0<-mean(Q_d[,1])
  lines(x=durations,y=as.numeric(apply(Q_d,2,mean)/Q_0),col=df_colori$values[df_colori$names%in%n],lwd=3)
}
par(xpd=T)
legend(x=25,y=1,legend=names(Qd_ds),col=colors,lwd=2,lty=1,cex=.8)  


#Peak-Duration curve model: Fiorentino's Formula, eq.13 Maione et al.,2003
exponential_model_1<-function(pars)
{
  tau<-pars[1]
  est<-tau/norm_durations*(1-exp(-norm_durations/tau))
  sumsq<-sum((est-obs)^2)
  return(sumsq)
}


#Regression for tau parameter
epsd_d_mean<-data.frame(matrix(NA,12,24,dimnames=list(x<-NULL,y<-c("Name",paste("value",2:24,sep="_")))))
epsd_d_mean$Name<-names(Qd_ds)

for(n in names(Qd_ds))
{
  Q_d<-Qd_ds[[n]][,-1]
  Q_0<-mean(Q_d[,1])
  epsd_d_mean[epsd_d_mean$Name%in%n,-1]<-as.numeric(apply(Q_d,2,mean)/Q_0)
}

epsd_d_pars<-data.frame(matrix(NA,12,2,dimnames=list(x=NULL,y=c("name","tau"))))
epsd_d_pars$name<-names(rd_ds)
durations<-c(.001,durations[-1])
for(n in names(Qd_ds))
{
  norm_durations<-durations/Tc$Tc[Tc$Name%in%n]
  obs<-as.numeric(epsd_d_mean[epsd_d_mean[,1]%in%n,-1])
  epsd_d_model<-optim(1,exponential_model_1,method = "Brent",lower = 0.001,upper = 100)
  epsd_d_pars[epsd_d_pars[,1]%in%n,-1]<-epsd_d_model$par
}
write.table(lambda_Tc_pars,"Parametri_lambdaTc.csv",row.names=FALSE,quote=T,sep=";",dec=",")
write.table(tau_Tc_pars,"Parametri_tauTc.csv",row.names=FALSE,quote=T,sep=";",dec=",")

#Regression for a and b parameter
x<-Tc$Tc[-which(Tc$Name%in%"SamoggiaCalcara")]
y<-epsd_d_pars$tau
y<-epsd_d_pars$tau[-which(Tc$Name%in%"SamoggiaCalcara")]
tau_Tc_pars <- data.frame(matrix(NA,nrow=1,ncol=3,dimnames = list(x=NULL,y=c("parameter","q","m"))))
tau_Tc_pars$parameter <- "tau"
tau_Tc_pars[,2:3]<-linear_model(epsd_d_pars[,2])$coefficients

df_colori<-data.frame(name=Tc$Name,value=colori)
df_colori<-read.csv("colors.csv",sep=";")

tau<-c(6,4.85,4.15,7.3,1.85,2,2.6,4,2.4,2.8,2.35,2.2)
plot(sort(Tc$Tc),tau)
x<-seq(min(Tc$Tc),max(Tc$Tc),.1)
Tc<-Tc %>% arrange(Tc)
df_colori<-df_colori[order(factor(df_colori[,1],levels=Tc$Name)),]
set_cex=1.5
par(mar=c(5.1, 6.1, 2.1, 2.1),mgp=c(4,1,0))
plot(NULL,
     xlim = c(3,max(Tc$Tc)),ylim=c(min(tau),max(tau)),
     xaxt="n",yaxt="n",
     #main=paste0("b = ",as.character(round(Esponenziale_b$par[1],2)),"+",as.character(round(Esponenziale_b$par[2],2)),"*exp(-",as.character(round(Esponenziale_b$par[3],2)),"*Tc)"),col.main="red",
     xlab="time of concentration (h)",ylab=expression(paste(tau," (-)")),cex.axis=set_cex,cex.lab=set_cex)
axis(side=1, at=seq(3,round(max(Tc$Tc))),cex.axis=set_cex)
axis(side=2, at=seq(1,7,1),cex.axis=set_cex,las=2)
abline(h=seq(1,7,1),col="gray",lty="solid",lwd=1)#griglia orizzontale
abline(v=seq(3,10,1),col="gray",lty="solid",lwd=1)#griglia vericale
points(Tc$Tc,tau,pch=16,col=df_colori$values,cex=1.8)#nuvola di punti Qcolmo/Qgg=f(Tc) (uno per ciscuna sezione)
x<-seq(from = 3,
       to   = max(Tc$Tc),
       by   = ((max(Tc$Tc) - 3)/(1000 - 1)))
lines(x,-0.598234389088741*x+7.99805402791602,col="cyan",lwd=3)
lines(x,.54+39.72*exp(-x/1.8),col="magenta",lwd=3)#Esponenziale_b$par[1]+Esponenziale_b$par[2]*exp(-Esponenziale_b$par[3]*x)#Modello esponenziale
text(Par_Tc$Tc+(max(Par_Tc$b)-min(Par_Tc$b))/40,Par_Tc$b+(max(Par_Tc$b)-min(Par_Tc$b))/40,Par_Tc$Code,cex=0.75)#Nomi delle sezioni


# Create a 2x1 layout: 2x2 plots + 1 column for legend
par(mar = c(6, 4, 4, 1),xpd=F)
x_poly<-seq(min(x),max(x),length.out=1000)
y_poly<-predict(linear_model(y),newdata = data.frame(x = x_poly))
plot(Tc$Tc,epsd_d_pars$tau,xlab="time of concentration (h)",ylab="a (-)",col=colors,pch=16,cex=1.2,cex.lab=1.2,cex.axis=1.2)
grid()
lines(x_poly,y_poly,col="red",lwd=2)
title(main = expression(epsilon[D] == (a/(D/Tc)) * (1 - exp(-(D/Tc)/a))))
legend("center",legend=names(rd_ds),col=colors,pch=16,cex=1.2)
mtext(expression(r[D] == (1/2)*exp(- (a * D/Tc) / (1 + b * D/Tc))),
      outer = TRUE, cex = 1.5, font = 2, line = -2.7)

#Ex. 
Tc=8 #8 ore

# Create a 2x1 layout: 2x2 plots + 1 column for legend
layout(matrix(c(1,2,3,3), nrow = 2, byrow = TRUE), 
       widths = c(1, 1, 0.4))
par(mar = c(6, 4, 4, 1),xpd=F)


#rd_d
a_test<-as.numeric(predict(a_model,newdata = data.frame(x = Tc)))
b_test<-as.numeric(predict(b_model,newdata = data.frame(x = Tc)))
durations<-seq(0,24,by=.2)
norm_durations<-durations/Tc
plot(norm_durations,.5*exp(-a_test*norm_durations/(1+b_test*norm_durations)),
     col="black",lwd=2,
     xlab="d/Tc (-)",
     ylab="rd (-)",
     ylim=c(0,.6),
     type="l")
grid()

#epsilond_d
tau_test<-as.numeric(predict(tau_model,newdata = data.frame(x = 8)))
plot(norm_durations,intercept*tau_test/norm_durations*(1-exp(-norm_durations/tau_test)),
     col="black",lwd=2,
     xlab="d/Tc (-)",
     ylab="rd (-)",
     ylim=c(0.4,1),
     type="l")

#Rising limb
#Q(t) = rd*D*Qd(T)
mtext(expression(r[D] == (1/2)*exp(- (a * D/Tc) / (1 + b * D/Tc))),
      outer = TRUE, cex = 1.5, font = 2, line = -2.7)

#Falling limb
#(1-rd)*D*Q(T)
Qmax=120 #m3/s
plot(norm_durations,Qmax*intercept*tau_test/norm_durations*(1-exp(-norm_durations/tau_test)),
     col="black",lwd=2,
     xlab="d (h)",
     ylab="Qd (m3/s)",
     ylim=c(40,Qmax),
     type="l")
grid()

#d(rd)/dd
drd_dd<--1/2*((a_test/Tc)/(1+b_test*norm_durations)^2)*exp((-a_test*norm_durations)/(1+b_test*norm_durations))
plot(norm_durations,drd_dd,
     col="black",lwd=2,
     xlab="duration (h)",
     ylab="rd (-)",
     ylim=c(-.07,0),
     type="l")

#d(εd)/dd
dQd_dd<-Qmax*tau_test*(1/norm_durations)*(-1/(Tc*norm_durations)*(1-exp(-norm_durations/tau_test))+(1/(tau_test*Tc)*exp(-norm_durations/tau_test)))
plot(norm_durations,dQd_dd,
     col="black",lwd=2,
     xlab="duration (h)",
     ylab="εd (-)",
     ylim=c(-10,1),
     type="l")

durations<-seq(8.001,.001,-.2)
durations<-seq(.001,16,.2)
norm_durations<-durations/Tc
drd_dd<--1/2*((a_test/Tc)/(1+b_test*norm_durations)^2)*exp((-a_test*norm_durations)/(1+b_test*norm_durations))
dQd_dd<-Qmax*tau_test*(1/norm_durations)*(-1/(Tc*norm_durations)*(1-exp(-norm_durations/tau_test))+(1/(tau_test*Tc)*exp(-norm_durations/tau_test)))
rd<-.5*exp(-a_test*norm_durations/(1+b_test*norm_durations))
Qd<-Qmax*tau_test/norm_durations*(1-exp(-norm_durations/tau_test))

grid()

par(mar = c(5, 4, 2,1),xpd=F)
plot(NULL,
     type="l",
     xlab="d (h)",
     ylab="streamflows (m3/s)",
     lwd=2,
     xlim=c(0,24),
     ylim=c(50,120))
grid()
lines(rev(durations),(drd_dd*durations*Qd+rd*Qd+rd*durations*dQd_dd)/(drd_dd*durations+rd),type="l",lwd=2)
lines(durations+Tc,(-drd_dd*durations*Qd+(1-rd)*Qd+(1-rd)*durations*dQd_dd)/(-drd_dd*durations+1-rd),type="l",lwd=2)
mtext(expression(paste("Ex: ","Tc = 8 ore, ",D[f]," = 24 ore, ",r[D[f]]," = 1/3, "),sep=""),font=2,cex=1.3,line=21)

