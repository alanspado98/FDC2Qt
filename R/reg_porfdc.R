#Dimensionless regional period of record flow duration curve

reg_adim_porfdc<-function(streamflows_GG,basin_descriptors,target_section,weighing_exponent)
{
  
  basin_descriptors<-basin_descriptors[basin_descriptors[,1]%in%unique(streamflows_GG[,1]),]
  rownames(basin_descriptors)<-basin_descriptors[,1]
  basin_descriptors<-basin_descriptors[,-1]
  
  #Principal Component Analysis
  #Removing descriptors with at least one negative or NA value
  Index_c<-which(sapply(basin_descriptors,function(column) any(is.na(column))))
  if(length(Index_c)>0) basin_descriptors <- basin_descriptors[,-Index_c]
  
  PCA <- stats::prcomp(basin_descriptors,scale = TRUE) #scale=TRUE: descriptors normalization (variance=1)
  PCA_sum<-summary(PCA)
  #Number of principal components; they have to explain at least 80% of the cumulative proportion of variance
  N_pc<-min(which(PCA_sum$importance[3,]>.8))
  PCA_x_N<-PCA$x[,1:N_pc]
  
  PCs_target <- PCA_x_N[rownames(PCA_x_N)%in%target_section,]#Target section principal components
  PCs_all    <- PCA_x_N[-which(rownames(PCA_x_N)%in%target_section),]#All the other principal components
  
  E_distance <- 0 #Euclidean distance of principal components in N_pc-dimensional space
  for (i in 1:N_pc) E_distance <- E_distance + (PCs_target[i] - PCs_all[,i])^2
  E_distance <- sqrt(E_distance)
  
  #weights to be assigned for each catchment
  exponent <- weighting_exponent
  weights <- 1/E_distance^exponent/sum(1/E_distance^exponent)

  #Regional adimensional PORFDC: weighting avarage of asimensional FDC of each gauged section
  #Weights must be applied to PORFDC with the same length: linear interpolation to the minimum length
  min_length<-min(rle(streamflows_GG[,1])$lengths)
  PORFDCs <- matrix(NA,dim(basin_descriptors)[1],min_length,dimnames=list(x<-rownames(basin_descriptors),y=NULL))

  #Weibull exceedance probability (final p-vector in the interpolation) 
  p_0<-1:(min_length)/(min_length+1)
  
  #Threshold duration: duration of the lowest streamflow value | Q/MAF>=1   
  d_threshold<-1
  
  for(r in rownames(PORFDCs))
  {
    str_r <- streamflows_GG[streamflows_GG[,1]%in%r,] 
    Q <- sort(str_r[,3],decreasing=T)
    
    #exceded probability vector (Weibull)
    p_r<-c(1:length(Q)/(length(Q)+1))
    
    adim_Q<-Q/mean(Q)
    
    if(d_threshold>p_r[which.min(adim_Q>=1)]) d_threshold <-p_r[which.min(adim_Q>=1)]
    
    #Resampling of the dimensionless PORFDC, so that all series have the same length as the shortest one
    PORFDCs[rownames(PORFDCs)%in%r,]<-approx(x=p_r,y=adim_Q,xout=p_0,method="linear")$y
  }
  ##RoI approach for weighting the dimensionless duration curves (regional PORFDCs, from empirical observations) 
  #in order to derive the dimensionless curve of the target catchment on a regional basis
  PORFDCs <- PORFDCs[-which(rownames(PORFDCs)==target_section),]
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

  # Store regional porfdc and d_threshold
  reg_FDC_results<-list(PORFDC=reg_PORFDC,d_threshold=d_threshold)
  return(reg_PORFDC)
}