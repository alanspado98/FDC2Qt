#Daily synthetic streamflow series

#' Title
#'
#' @param basin_descriptors 
#' @param target_section 
#'
#' @return
#' @export
#'
#' @examples
daily_streamflows<-function(basin_descriptors,target_section)
{
  reg_PORFDC<-use_internal_data(list_name = "regPORFDC",dataset_name = "PORFDC")#Regional PORFDC
  reg_LT_str<-use_internal_data(list_name = "regPORFDC",dataset_name = "LT")#Regional streamflow
  weights<-use_internal_data(list_name = "regPORFDC",dataset_name = "weights")#Weights
  LT_streamflow<-use_internal_data(list_name = "regPORFDC",dataset_name = "LT")#Daily streamflows
  
  #Donator site
  #1) if gauged -> R2
  #2) if not gauged -> descriptors
  
  #1) Gauged target section
  if(!is.na(match(target_section,unique(streamflows[,1]))))
  {
    str_no_target<- streamflows[-which(streamflows[,1]%in%target_section),]
    str_target   <- streamflows[streamflows[,1]%in%target_section,]
    
    str_no_target<-str_no_target[str_no_target[,2]%in%str_target[,2],]
    
    R2<-rep(NA,length(unique(str_no_target[,1])));names(R2)<-unique(str_no_target[,1])
    for(r in unique(str_no_target[,1]))
    {
      str_no_target_r<-str_no_target[str_no_target[,1]%in%r,]
      str_target_r<-str_target[str_target[,2]%in%str_no_target_r[,2],]
      R2[names(R2)%in%r]<-cor(str_target_r[,3],str_no_target_r[,3])^2
    }
    donator_section<-names(which.max(R2))
  }else#2) if not gauged -> descriptors
  {
    donator_section<-names(which.max(weights))
  }
  
  #PORFDC donator site
  str_donator<-streamflows[streamflows[,1]%in%donator_section,]
  Q <- sort(unique(str_donator[,3]),decreasing=T)
  donator_PORFDC <- data.frame(p=c(1:length(Q)/(length(Q)+1)),Q)
  
  #Resampling regional PORFDC according to PORFDC of donator site
  reg_PORFDC[,2]<-reg_PORFDC[,2]*LT_streamflow
  reg_PORFDC <- data.frame(p=donator_PORFDC$p,Q=approx(x=reg_PORFDC$p,y=reg_PORFDC$Q,xout=donator_PORFDC$p,method="linear")$y)
  
  #Synthetic hydrograph for the target site  
  #Dataframe intialization
  reg_daily_str <- data.frame(matrix(NA,ncol = 2,nrow = length(str_donator[,1]),dimnames=list(x=NULL,y=c("Date","Value"))))
  reg_daily_str[,1] <- str_donator[,2]
  
  for(i in 1:length(reg_daily_str[,1]))
  { 
    index<-which(donator_PORFDC[,2]==str_donator[i,3])#index in PORFDC of the donator site corresponding to i-th daily streamflow
    reg_daily_str[i,2]<-reg_PORFDC[index,2]#corrisponding streamflow in regional PORFDC
  }
  
  # Store regional daily streamflows inside a list
  .internal_env$reg_streamflows[["GG"]] <- reg_daily_str
}