#Adimensional regional period of record flow duration curve
target_section<-"Marecchia@RiminiSS16"
load("data/basin_descriptors.rda")
reg_adim_porfdc<-function(basin_descriptors,target_section)
{
  streamflows<-use_internal_data(list_name = "streamflows",dataset_name = "GG")
  
  basin_descriptors<-basin_descriptors[basin_descriptors[,1]%in%unique(streamflows[,1]),]
  rownames(basin_descriptors)<-basin_descriptors[,1]
  basin_descriptors<-basin_descriptors[,-1]
  
  #Principal Component Analysis
  #Removing descriptors with at least one negative or NA value
  Index_c<-which(sapply(basin_descriptors,function(column) any(is.na(column))))
  if(length(Index_c)>0) basin_descriptors <- basin_descriptors[,-Index_c]
  
  PCA <- stats::prcomp(basin_descriptors,scale = TRUE) #scale=TRUE: descriptors normalization (variance=1)
  PCA_sum<-summary(PCA)
  #Number of principal components; they have to "explain" at least 75% of cumulative proportion of variance
  N_pc<-min(which(PCA_sum$importance[3,]>.8))
  
  #Values of principal components truncated to N_pc
  PCA_x_N<-PCA$x[,1:N_pc]
  par(xpd=F)
  par(mar=c(5,5,3,1))
  barplot(PCA_sum$importance[3,1:10],
          xlim=c(0,12.5),ylim=c(0,1),las=2,col="blue",
          xlab="Principal component (-)",ylab = "Cumulative proportion of variance (-)",cex.axis=1,xaxs="i",yaxt="n")
  legend("topleft",legend="0.8",lty=2,col="red",lwd=2,box.col ="transparent") 
  axis(2,at=seq(0,1,.1),labels=seq(0,1,.1),cex.axis=1,las=2)
  abline(h=0,col="black")
  abline(h=.8,col="red",lty=2,lwd=2)
  
  PCs_target <- PCA_x_N[rownames(PCA_x_N)%in%target_section,]#Target section principal components
  PCs_all    <- PCA_x_N[-which(rownames(PCA_x_N)%in%target_section),]#All the other principal components
  
  E_distance <- 0 #Euclidean distance of principal components in N_pc-dimensional space
  for (i in 1:N_pc) E_distance <- E_distance + (PCs_target[i] - PCs_all[,i])^2
  E_distance <- sqrt(E_distance)
  
  #Calcolo dei pesi da assegnare a ciscun bacino
  exponent <- 3 #Choice of the exponent
  weights <- 1/E_distance^exponent/sum(1/E_distance^exponent)

  #Regional adimensional PORFDC: weighting avarage of asimensional FDC of each gauged section
  #Weights must be applied to PORFDC with the same length: linear interpolation to the minimum length
  min_length<-min(rle(streamflows[,1])$lengths)
  PORFDCs <- matrix(NA,dim(basin_descriptors)[1],min_length,dimnames=list(x<-rownames(basin_descriptors),y=NULL))

  #Weibull exceedance probability (final p-vector in the interpolation) 
  p_0<-1:(min_length)/(min_length+1)
  
  #Threshold duration: duration of the lowest streamflow value | Q/Qmean>=1   
  d_threshold<-1
  
  for(r in rownames(PORFDCs))
  {
    str_r <- streamflows[streamflows[,1]%in%r,]#Selezione della stazione di interesse 
    Q <- sort(str_r[,3],decreasing=T) #Ordinamento descrescente

    # #Salvataggio della PORFDC dimensionale per il sito target
    if(r==target_section)
    {
      PORFDC_target <- Q/mean(Q)
      PORFDC_target <- cbind(c(1:length(Q)/(length(Q)+1)),PORFDC_target)
      colnames(PORFDC_target) <- c("p","Q")
    }
    
    #Vettore probabilita' di superamento (Weibull)
    p_r<-c(1:length(Q)/(length(Q)+1))
    
    adim_Q<-Q/mean(Q)
    
    if(d_threshold>p_r[which.min(adim_Q>=1)]) d_threshold <-p_r[which.min(adim_Q>=1)]
    
    #Ricampionamento della PORFDC adimensionale, in modo che tutte le serie abbiamo lunghezza pari a quella 
    #di lunghezza minima
    PORFDCs[rownames(PORFDCs)%in%r,]<-approx(x=p_r,y=adim_Q,xout=p_0,method="linear")$y
  }
  ##Approccio RoI per la pesatura delle curve di durata adimensionali (PORFDC regionali, da osservazioni empiriche) al fine di ottenere
  #la curva adimensionale del bacino target stimata su base regionale
  PORFDCs <- PORFDCs[-which(rownames(PORFDCs)==target_section),]#Non si considera MarecchiaSS16 nella pesatura nel caso in cui sia il target(il peso sarebbe 0...)
  reg_PORFDC <- as.numeric(t(PORFDCs)%*%weights)
  
  reg_PORFDC<-unique(reg_PORFDC) #remove streamflows redundancies 
  reg_PORFDC<-data.frame(p=1:length(reg_PORFDC)/(length(reg_PORFDC)+1),Q=reg_PORFDC)

  #Extension of regional PORFDC in order to be applicable for every durations in continuous daily streamflow generation
  reg_PORFDC<- rbind(rep(NA,2),reg_PORFDC,rep(NA,2))
  
  #First data: maximum adimensional streamflow overall the gauged sections, associated to 0 exceedence probabiliy
  reg_PORFDC$p[1] <- 0
  reg_PORFDC$Q[1] <- max(PORFDCs[,1])
  #Last data: minimum adimensional streamflow of the target section, associated to a 1 exceedence probability
  reg_PORFDC$p[length(reg_PORFDC$p)] <- 1
  reg_PORFDC$Q[length(reg_PORFDC$Q)] <- reg_PORFDC$Q[length(reg_PORFDC$Q)-1]
  
  ##Figura - PORFDC adimensionali per i siti "donatori" e PORFDC adimensionale per il sito di interesse stimata su base regionale 
  
  par(mar=c(5,6,2,1))
  par(mgp=c(3.5,1,0))
  options(scipen=999)#No scientific axis notation
  
  #Adimensional flow duration cruve of gauged sections
  plot(NULL,type="l",lty=1,lwd=3,col="black",
       log="y",
       yaxt="n",
       xlim=c(0,1),
       ylim=c(0.001,max(PORFDCs[,1])),
       xaxt="n",
       main="Adimensional PORFDCs",
       xlab="duration (-)",ylab=expression(paste("streamflow (",m^3,"/s)",sep="")),cex.lab=1)
  axis(1,at=seq(0,1,0.1),cex.axis=1)#Inserimento tick marks asse x
  axis(2,at=c(0.01,0.1,1,10),cex.axis=1,las=2)#Inserimento tick marks asse y
  abline(v=seq(0,1,0.1),col="lightgray",lty="dotted",lwd=0.1)
  abline(h=c(0.01,0.1,1,10),col="lightgray",lty="dotted",lwd=0.1)
  
  for(r in rownames(PORFDCs)) lines(p_0,PORFDCs[rownames(PORFDCs)%in%r,],lty=1,lwd=1,col="grey")
  lines(reg_PORFDC,lty=1,lwd=2,col="black")
  lines(PORFDC_target,lty="dotted",lwd=2,col="black")                          
  
  legend(.7,max(PORFDCs[,1])-.001,
         title="target site",legend=c("regional PORFDC","empirical PORFDC"),
         lwd=c(3,3,2),col=c("black","black"),
         lty=c("solid","dotted"),title.font = 2, title.adj =.1,cex=.8)
  
  legend(.7,(max(PORFDCs[,1])-.001)*.12,
         legend=expression(bold("gauged sites")),
         lwd=2,col=rgb(.6,.6,.6),
         lty="solid",cex=.8)
  
  # Store all dataframes inside a list
  .internal_env$regPORFDC[["weights"]] <- weights
  .internal_env$regPORFDC[["PORFDC"]] <-  reg_PORFDC
  
  #Store
  .internal_env$regstreamflows[["d_threshold"]] <-  d_threshold
}