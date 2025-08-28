#' Estimation of total duration of flood events for a given river site
#'
#' @return
#' @export
#'
#' @examples

#Df e Rdf estimation
library(dplyr)
library(lubridate)
library(pracma)
library(pals)
  
setwd("..")
setwd("..")
#Check if we have hourly time step
d_threshold<-as.numeric(read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/Durata_soglia_Downscaling.csv",sep=";",dec=","))
file_names<-dir("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/CSV_HH")
#Caricamento del dataset delle portate medie giornaliere definitivo (Capitolo 2) per limare le serie di portata oraria delle sezioni che hanno una nuvola di punti Qcolomo/Qgg=f(Durata) sospetta
QMG <- read.csv("02_Riempimento&Revisione_QMG/02_2_Revisione_QMG/Dataset_QMG_definitivo.csv",sep=";",dec=",",check.names = F)

Df <- list()
rDf<- list()
#rising_limb<-list()
for(fn in file_names)
{
  Q_hh<-read.csv(paste0("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/CSV_HH/",fn),sep=";",dec=",")
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
  
  Qs_d <- Q_gg
  Qs_d$d<-rep(NA,length(Qs_d[,1]))
  
  for(q in porfdc$Q) Qs_d$d[Qs_d$VALORE%in%q] <- porfdc$p[porfdc$Q==q] #daily streamflow durations
  d_threshold<-porfdc$p[which.min(abs(porfdc$Q/mean(porfdc$Q)-1))]/10
  
  index<-c()
  for(i in 2:(length(Qs_d$DATA)-1))
  {
    if(Qs_d$DATA[i]-Qs_d$DATA[i-1]==1 
       && Qs_d$DATA[i+1]-Qs_d$DATA[i]==1)
    {
      if(Qs_d$VALORE[i]>Qs_d$VALORE[i-1] && 
         Qs_d$VALORE[i]>Qs_d$VALORE[i+1] && 
         Qs_d$d[i]<=d_threshold) index<-append(index,i)
    }
  }
  index<-unique(index)
  Qs_d<-Qs_d[index,]
  
  Df_site <- data.frame(matrix(NA,nrow=length(index),ncol=3,dimnames = list(x<-NULL,y<-c("date","value","t_start"))))
  rDf_site<- data.frame(matrix(NA,nrow=length(index),ncol=2,dimnames = list(x<-NULL,y<-c("date","value"))))
  Df_site$date <- Qs_d$DATA
  rDf_site$date<- Qs_d$DATA
  
  #rising_limb_site <- data.frame(matrix(NA,nrow=length(index),ncol=2,dimnames = list(x<-NULL,y<-c("date","value"))))  
  #rising_limb_site$date<-Qs_d$DATA
  
  for(j in 1:length(Qs_d$DATA))
  {
    i_start<- which(Q_hh$DATAORA%in%paste0(Qs_d$DATA[j]-1," 12:00:00"))
    i_end  <- which(Q_hh$DATAORA%in%paste0(Qs_d$DATA[j]+1," 12:00:00"))
    if(length(i_start)==0 || length(i_end)==0) next
    Q_hh_j<-Q_hh[i_start:i_end,]
    peaks<-findpeaks(Q_hh_j$VALORE,threshold=0,zero="-",sortstr=FALSE)#Minimum difference between a peak and its neighbors
    peaks<-data.frame(peaks)
    colnames(peaks)<-c("Qc","i_Qc","i_start","i_end")
    peaks<-cbind.data.frame(peaks,pos=1:length(peaks[,1]))
    peaks<-peaks[order(peaks[,1],decreasing = T),]
    
    #Piu' colmi
    if(dim(peaks)[1]>1)
    {
      peaks_event <- data.frame(matrix(,nrow=0,ncol=5,dimnames=list(x=NULL,y=colnames(peaks))))
      k<-0
      flag=T
      while(flag)
      {
        k<-k+1
        Qc<-as.numeric(peaks[k,1])
        if(abs(Qc/max(peaks[,1])-1)<=.4 &&
           length(which(abs(Q_hh_j$VALORE[min(c(peaks[k,2],peaks[which.max(peaks[,1]),2])):max(c(peaks[k,2],peaks[which.max(peaks[,1]),2]))]/peaks[which.max(peaks[,1]),1]-1)>.5))==0)
          peaks_event<-rbind.data.frame(peaks_event,peaks[k,])

        if(k==dim(peaks)[1]) 
        {
          flag=F
          #Index_max<-which(peaks_event[,1]==max(peaks[,1]))
          #peaks_event<-peaks_event[-Index_max[-1],]
        }  
      }
    }else peaks_event<-peaks #Un solo colmo

    peaks_event<-peaks_event[order(peaks_event[,5]),]
    
    # plot(Q_hh_j,type="l",xaxt="n",xlab="ore",ylab="portata (m3/s)")
    #  axis(1,at=Q_hh_j$DATAORA,label=Q_hh_j$DATAORA,las=2)
    #  points(Q_hh_j$DATAORA[peaks_event[which.max(peaks_event[,1]),2]],max(peaks_event[,1]),col="red",pch=16)
    # points(Q_hh_j$DATAORA[c(peaks_event[1,3],peaks_event[length(peaks_event[,1]),4])],
    #        Q_hh_j$VALORE[c(peaks_event[1,3],peaks_event[length(peaks_event[,1]),4])],col="green",pch=16)
     #abline(h=Q_threshold,col="blue")
    Df_site$value[j]<-peaks_event[length(peaks_event[,1]),4]-peaks_event[1,3]
    Df_site$t_start[j]<-as.character(Q_hh_j$DATAORA[peaks_event[1,3]])
    i_Qc <- round(peaks_event[,2]%*%peaks_event[,1]/sum(peaks_event[,1]))
    #abline(v=Q_hh_j$DATAORA[i_Qc],col="blue")
    rDf_site$value[j]<-(i_Qc-peaks_event[1,3])/(peaks_event[length(peaks_event[,1]),4]-peaks_event[1,3])
  } 
  Df <- append(Df,list(Df_site))
  rDf<- append(rDf,list(rDf_site))
}
names(Df) <-sub("\\..*","",sub("^[^_]*_","",file_names))
names(rDf)<-sub("\\..*","",sub("^[^_]*_","",file_names))
Df<-Df[1:12]
mean(rDf_site$value,na.rm=T)
#Rising limb
rising_limb_mean<-rep(NA,length(Df))
names(rising_limb_mean)<-names(Df)
for(i in 1:length(Df)) rising_limb_mean[i]<-mean(Df[[i]]$value*rDf[[i]]$value,na.rm = T)

save(Df,file="Df_event.RData")
save(rDf,file="rDf_event.RData")

#Time of concentration
Tc<-read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/abc_sezioni_chiusura.csv",sep=";",dec=",")
Tc<-Tc[,c(2,6)]
Tc<-Tc[order(factor(Tc[,1],levels=names(Df))),]

layout(matrix(c(1,2), nrow = 1, byrow = TRUE), 
       widths = c(1, 1, 0.4))

Df_mean<-rep(NA,length(Df))
names(Df_mean)<-names(Df)
for(n in names(Df)) Df_mean[names(Df_mean)%in%n]<-mean(Df[[n]][,2],na.rm=T)

df<-cbind.data.frame(x=Tc$Tc,y=rising_limb_mean)
model <- lm(y ~ x,data=df)#poly(x, 2, raw = TRUE)
summary(model)
par(mar = c(5,5,2,12),xpd=F)
colors<-kelly()[c(19,18,16,13,12,10,9,7,5,4,3,2)]
plot(Tc$Tc,Df_mean,xlab="time of concentration (h)",ylab=expression(paste(D[f]," (h)",sep="")),col=colors,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5)

par(xpd=T)
legend(x=10.9,y=max(rising_limb_mean),legend=names(Df),col=colors,pch=16,cex=1.1)

rDf_mean<-rep(NA,length(rDf))
names(rDf_mean)<-names(rDf)
for(n in names(Df)) rDf_mean[names(rDf_mean)%in%n]<-mean(rDf[[n]][,2],na.rm=T)

par(mar = c(5,5,2,12),xpd=F)
colors<-kelly()[c(19,18,16,13,12,10,9,7,5,4,3,2)]
plot(basin_descriptors$`DrainageDensity_1/km`[-which(rownames(basin_descriptors)%in%"SanternoBTossignano")],rDf_mean[-which(names(rDf_mean)%in%"SanternoBTossignano")],xlab="time of concentration (h)",ylab=expression(paste(r[D[f]]," (-)",sep="")),col=colors,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5)
par(xpd=T)
legend(x=10,y=max(rDf_mean),legend=names(Df),col=colors,pch=16,cex=1.1)
plot(Tc$Tc,rDf_mean,xlab="area",ylab=expression(paste(D[f]," (h)",sep="")),col=colors,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5)

#Morfological descriptors
load("FDC2Qt_datasets/basin_descriptors.RData")
basin_descriptors$Name<-gsub("@","",basin_descriptors$Name)
basin_descriptors$Name[3]<-"SavioSCarlo"
basin_descriptors$Name[4]<-"LamoneStCasale"
basin_descriptors$Name[5]<-"SanternoBTossignano"
basin_descriptors$Name[6]<-"MarecchiaSS16"

basin_descriptors<-basin_descriptors[basin_descriptors$Name%in%names(Df),]
basin_descriptors<-basin_descriptors[order(factor(basin_descriptors[,1],levels=names(Df))),]
rownames(basin_descriptors)<-basin_descriptors$Name
basin_descriptors<-basin_descriptors[,-1]
#Trasformazione in scala log di portate e descrittori (poiche' si ipotizza un legame di potenza tra le due grandezze)
Index_c<-which(sapply(basin_descriptors,function(column) any(column<=0)))#Negative values
Index_c<-append(Index_c,which(sapply(basin_descriptors,function(column) any(is.na(column)))))#NA values
if(length(Index_c)>0) basin_descriptors <- basin_descriptors[,-Index_c] 

x<-basin_descriptors[-which(rownames(basin_descriptors)%in%"SanternoBTossignano"),]
y<-Df_mean[-which(names(Df_mean)%in%"SanternoBTossignano")]


# x<-x[,-which(colnames(x)%in%"1stOrderStreamSlopes_%")]
# x<-x[,-which(colnames(x)%in%"1stOrderStreamLengths_km")]
#x<-x[,-which(colnames(x)%in%"2ndOrderStreamLengths_km")]
# x<-x[,-which(colnames(x)%in%"CentroidNorth_m")]
# x<-x[,-which(colnames(x)%in%"CentroidEst_m")]
x<-x[,-which(colnames(x)%in%"DrainageAreaRatio")]
# x<-x[,-which(colnames(x)%in%"BifurcationRatio")]
x<-x[,-which(colnames(x)%in%"StreamSlopeRatio")]
x<-x[,-which(colnames(x)%in%"ElongationRatio")]
# x<-x[,-which(colnames(x)%in%"1stOrderDrainageAreas_km2")]
# x<-x[,-which(colnames(x)%in%"2ndOrderDrainageAreas_km2")]
# x<-x[,-which(colnames(x)%in%"3rdOrderDrainageAreas_km2")]
# x<-x[,-which(colnames(x)%in%"2ndOrderStreamSlopes_%")]
# x<-x[,-which(colnames(x)%in%"3ndOrderStreamSlopes_%")]
# x<-x[,-which(colnames(x)%in%"2ndOrderStreamLengths_km")]
# x<-x[,-which(colnames(x)%in%"3rdOrderStreamLengths_km")]
# x<-x[,-which(colnames(x)%in%"WidthFunction5%_m")]
# x<-x[,-which(colnames(x)%in%"WidthFunction10%_m")]
# x<-x[,-which(colnames(x)%in%"WidthFunction15%_m")]
# x<-x[,-which(colnames(x)%in%"WidthFunction40%_m")]
# x<-x[,-which(colnames(x)%in%"WidthFunction50%_m")]
x<-x[,-which(colnames(x)%in%"StDevMAEELast20_mm")]
x<-x[,-which(colnames(x)%in%"StDevMAPELast20_mm")]
x<-x[,-which(colnames(x)%in%"MeanMAPELast20_mm")]
# x<-x[,-which(colnames(x)%in%"MeanMAAELast20_mm")]
# x<-x[,-which(colnames(x)%in%"HypsographicCurve90%_masl")]
# x<-x[,-which(colnames(x)%in%"HypsographicCurve95%_masl")]
# x<-x[,-which(colnames(x)%in%"HypsographicCurve75%_masl")]

#R2 adjusted metrics temporary variable 
R2_temp_1  <- 0  #Highest R2 adjusted of the k-th while cycle
R2_temp_2  <- -1 #Highest R2 adjusted of the (k-1)-th while cycle
#lm summary temporaty variables
Model_best_sum_1 <- list() #Summary of the highest R2 adjusted of the k-th while cycle
Model_best_sum_2 <- list() #Summaries of the highest R2 adjusted of the (k-1)-th while cycle
#Model descriptors
Best_descriptors <- matrix(,nrow=length(x[,1]),ncol=0,dimnames=list(x=rownames(x),y=NULL))
N<-2
while(R2_temp_1 > R2_temp_2 && length(Best_descriptors[1,])<N)#repeat while adding another descriptor increase R2 adjusted metrics 
{
  R2_temp_2    <- R2_temp_1 #Highest adjusted R2 (best model) of the previous cycle
  R2_temp_1    <- 0 #Temporary variable update
  
  #Searching for a new descriptor between the residual ones
  for (i in 1:length(x[1,]))
  {
    if(R2_temp_2 == 0)#Enter in the for cycle for the first time
    {
      Model <- lm(y ~ x[,i])#Linear model 
    }else if(R2_temp_2 != 0)#Enter in the for cycle from the second time
    {
      Model <- lm(y ~ Best_descriptors + x[,i])#Multiple linear regression
    }
    Model_sum <- summary(Model) #Model summary
    
    if(is.na(Model_sum$adj.r.squared)) break #Esc while cycle in case R2 adjusted is NA 
    if(Model_sum$adj.r.squared > R2_temp_1)
    {
      R2_temp_1 <- Model_sum$adj.r.squared #Update R2 adjusted
      Model_best_sum_1 <- Model_sum #Update summary 
      Pos_best_descriptor <- i #index i-th descriptor (in order to save it in "Best_descriptors" and remove it from "x")
    }
  }
  
  #Adding the best descriptor found in the previous iteration 
  Best_descriptors <- cbind(Best_descriptors,x[,Pos_best_descriptor])
  colnames(Best_descriptors)[which(is.null(colnames(Best_descriptors))==TRUE)]<- colnames(x)[Pos_best_descriptor]#First column 
  colnames(Best_descriptors)[which(colnames(Best_descriptors)=="")] <- colnames(x)[Pos_best_descriptor]#Column consecutive to the first one: "" 
  x<-x[,-Pos_best_descriptor]#Removal of the best descriptors from the dataset
  Model_best_sum_2 <- append(Model_best_sum_2,Model_best_sum_1)#Appending best model summary of the previous while cycle
}


df<-cbind.data.frame(x1=Best_descriptors[,1],x2=Best_descriptors[,2],y=y)
model<-lm(y ~x1+x2,data = df)
summary(model)
y_est<-predict(model,newdata = data.frame(x1=df$x1,x2=df$x2))
par(mar = c(5,5,2,12),xpd=F)
plot(df$y,y_est,xlim=c(min(c(df$y,y_est)),max(c(df$y,y_est))),ylim=c(min(c(df$y,y_est)),max(c(df$y,y_est))),
     col=colors,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5,
     xlab="",
     ylab="")
colnames(Best_descriptors)
title(main = bquote(bold(D[f]) == .(round(model$coefficients[1],1)) * .(round(model$coefficients[2],1)) * bold(ShapeFactor) + .(round(model$coefficients[3],2)) * bold(MeanSlope2)))
abline(a=0,b=1,col="red",lwd=2)
par(xpd=T)
legend(x=max(c(min(c(df$y,y_est)),max(c(df$y,y_est))))*1.02,y=max(c(min(c(df$y,y_est)),max(c(df$y,y_est)))),legend=names(rDf)[-8],col=colors[-8],pch=16,cex=1)

#Quando iniziano gli eventidi piena?
Df$MarzenoRivalta
