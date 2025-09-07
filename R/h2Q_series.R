#' Streamflow estimates from observed stream stages and rating curves
#'
#' @param stream_stages a data frame with three columns: station name
#'   (character), time istant of observation (POSIXct), stream stage
#'   observation (numeric)
#'   
#' @param rating_curves a data frame with six columns: station name
#'   (character), start date of applicability (character), end date of
#'   applicability (character), minimum applicable stage (numeric), maximum
#'   applicable stage (numeric), rating curve equation (character)
#'   
#' @param streamflows a data frame with three columns: station name
#'   (stations in common with stream_stages must match), time istant of
#'   observation (character), streamflow observation (numeric)
#'   
#' @param availability_matrices a list with `stream_stages` and `obs_streamflow` availability matrix,
#'   `rating_curves` binary matrix
#'
#' @return
#' @export
#'
#' @examples

streamflows_eval<-function(stream_stages,rating_curves,streamflows,hl_DA_matrix,rc_DA_matrix,str_DA_matrix) 
{
  #Date-Time Strings standardization
  stream_stages[,2] <- std_datetime(stream_stages[,2])
  rating_curves[,2] <- std_datetime(rating_curves[,2])
  rating_curves[,3] <- std_datetime(rating_curves[,3])
  streamflows[,2]   <- std_datetime(streamflows[,2])
  
  #Hourly series generation (section that has at least one stream stages and a rating curve)
  `%>%` <- magrittr::`%>%`
  streamflows_HH<-data.frame(matrix(NA,0,3,dimnames=list(x=NULL,y=colnames(streamflows))))
  streamflows_GG<-data.frame(matrix(NA,0,3,dimnames=list(x=NULL,y=colnames(streamflows))))
  
  #Section with hydrometric levels and rating curves
  stations<-unique(stream_stages[,1])[unique(stream_stages[,1])%in%rownames(rc_DA_matrix)[which(apply(rc_DA_matrix,1,sum)>0)]]

  for(r in stations) 
  {
    years_rc<-as.numeric(names(which(rc_DA_matrix[rownames(rc_DA_matrix)%in%r,]==1)))
    hl_r<-stream_stages[stream_stages[,1]%in%r & lubridate::year(stream_stages[,2])%in%years_rc,-1]
    rc_r<-rating_curves[grep(gsub("([a-z])([A-Z])","\\1 \\2",strsplit(r,"@")[[1]][2]),rating_curves[,1],ignore.case = T),]
    str_eval_r<-rating_curve_eval(hl_r,rc_r)
  
    #Hourly streamflows
    str_eval_r<- str_eval_r %>%
      dplyr::group_by(Date = as.Date(Datetime), Hour = lubridate::hour(Datetime)) %>%
      dplyr::summarise(Value = mean(Value, na.rm = TRUE),.groups="drop") %>%
      as.data.frame()
   
    if(length(which(is.na(str_eval_r[,2])))>0) str_eval_r<-str_eval_r[-which(is.na(str_eval_r[,2])),]
    str_r_HH<-str_eval_r
    str_r_HH[,2]<-as.POSIXct(paste0(str_r_HH[,1]," ",str_r_HH[,2],":00"),format="%Y-%m-%d %H:%M",tz="GMT")
    colnames(str_r_HH)[2]<-"Datetime"
    str_r_HH[,1]<-r
    colnames(str_r_HH)[1]<-"Name"
    
    streamflows_HH<-rbind.data.frame(streamflows_HH,str_r_HH)

    #Daily streamflows
    str_eval_r<- str_eval_r %>%
      dplyr::group_by(Date) %>%
      dplyr::summarise(Value= mean(Value, na.rm = TRUE),.groups = "drop")%>%
      as.data.frame

    str_r_GG<-streamflows[streamflows[,1]%in%r,-1]
    colnames(str_r_GG)[1]<-"Date"
    
    #Filling missing values
    str_r_GG[is.na(str_r_GG[,2]),2][str_r_GG[is.na(str_r_GG[,2]),1]%in%str_eval_r[,1]] <- str_eval_r[str_eval_r[,1]%in%str_r_GG[is.na(str_r_GG[,2]),1],2]

    #Adding values outside the temporal window of streamflows
    str_r_GG<-rbind.data.frame(str_r_GG,str_eval_r[!str_eval_r[,1]%in%str_r_GG[,1],])
    str_r_GG<-str_r_GG %>% dplyr::arrange(Date)

    #Linear interpolation if less than 7 consecutive NA daily streamflows
    NA_streamflows_r <- split(which(is.na(str_r_GG[,2])),cumsum(c(1,diff(which(is.na(str_r_GG[,2])))!=1)))
    for(n in names(NA_streamflows_r))
    {
      i_start <- NA_streamflows_r[[n]][1]-1 #First day
      i_end <- NA_streamflows_r[[n]][length(NA_streamflows_r[[n]])]+1 #Last day
      
      if(i_end-i_start-1<15 && !is.na(str_r_GG[i_start,2]) && !is.na(str_r_GG[i_end,2]))
      {
        str_r_GG[i_start:i_end,2] <- approx(x=c(i_start,i_end),y=c(str_r_GG[i_start,2],str_r_GG[i_end,2]),method="linear",n=i_end-i_start+1,na.rm=TRUE)$y 
      }
    }
    year_to_remove <- unique(lubridate::year(str_r_GG[is.na(str_r_GG[,2]),1]))
    if(length(year_to_remove)>0) str_r_GG<-str_r_GG[-which(lubridate::year(str_r_GG[,1])%in%year_to_remove),]
    str_r_GG<-cbind.data.frame(Name=r,str_r_GG)
    colnames(str_r_GG)[2]<-"Datetime"
    streamflows_GG<-rbind.data.frame(streamflows_GG,str_r_GG)
  } 
  
  streamflows_LT <-streamflows[is.na(streamflows[,2]),-2]
  rownames(streamflows_LT)<-NULL
  
  # Store all dataframes inside a list
  streamflows_h2Q <- list(HH=streamflows_HH,GG=streamflows_GG,LT=streamflows_LT)
  
  return(streamflows_h2Q)
}
