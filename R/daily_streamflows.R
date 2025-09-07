#Daily synthetic streamflow series

#' Title
#'
#' @param basin_descriptors 
#' @param target_section 
#'
#' @return
#' @export
#'
#' @example

daily_streamflows<-function(reg_PORFDC,streamflows_GG,reg_MAF,target_section,hydro_consistency)
{
  #Donator site, choice depens on the aim:
  #1) hydrologically consistent long streamflow series -> longest series within the region
  #2) reconstruct historical streamflow series (target site must be gauged!) -> R2 
  hydro_consistency<-T
  
  if(hydro_consistency) donor_section<-names(which.max(table(streamflows_GG[,1])))

  if(!hydro_consistency && !is.na(match(target_section,unique(streamflows_GG[,1]))))
  {
    str_no_target<- streamflows_GG[-which(streamflows_GG[,1]%in%target_section),]
    str_target   <- streamflows_GG[streamflows_GG[,1]%in%target_section,]
    
    str_no_target<-str_no_target[str_no_target[,2]%in%str_target[,2],]
    
    R2<-rep(NA,length(unique(str_no_target[,1])));names(R2)<-unique(str_no_target[,1])
    for(r in unique(str_no_target[,1]))
    {
      str_no_target_r<-str_no_target[str_no_target[,1]%in%r,]
      str_target_r<-str_target[str_target[,2]%in%str_no_target_r[,2],]
      R2[names(R2)%in%r]<-cor(str_target_r[,3],str_no_target_r[,3])^2
    }
   donor_section<-names(which.max(R2))
  }
  
  #PORFDCdonor site
  str_donor<-streamflows_GG[streamflows_GG[,1]%in%donor_section,]
  Q <- sort(unique(str_donor[,3]),decreasing=T)
  donor_PORFDC <- data.frame(p=c(1:length(Q)/(length(Q)+1)),Q)

  #Resampling regional PORFDC according to PORFDC of donor site
  reg_PORFDC[,2]<-reg_PORFDC[,2]*reg_MAF
  reg_PORFDC <- data.frame(p=donor_PORFDC$p,Q=approx(x=reg_PORFDC$p,y=reg_PORFDC$Q,xout=donor_PORFDC$p,method="linear")$y)
  
  #Synthetic hydrograph for the target site  
  #Dataframe intialization
  reg_daily_str <- data.frame(matrix(NA,ncol = 2,nrow = length(str_donor[,1]),dimnames=list(x=NULL,y=c("Date","Value"))))
  reg_daily_str[,1] <- str_donor[,2]
  
  for(i in 1:length(reg_daily_str[,1]))
  { 
    index<-which(donor_PORFDC[,2]==str_donor[i,3])#duration in the donor site PORFDC corresponding to i-th daily streamflow
    reg_daily_str[i,2]<-reg_PORFDC[index,2]#streamflow in regional PORFDC corrisponding to the duration of the i-th daily streamflow
  }

  # regional daily streamflows dataframe
  return(reg_daily_str)
}