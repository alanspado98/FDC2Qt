
#' Title
#'
#' @return
#' @export
#'
#' @examples
streamflows_eval<-function() 
{
  hydrometric_levels<-use_internal_data(list_name = "input_data",dataset_name = "stream_stages")
  streamflows<-use_internal_data(list_name = "input_data",dataset_name = "obs_streamflows")
  rating_curves<-use_internal_data(list_name = "input_data",dataset_name = "rating_curves")

  #Hourly series generation (section that has at least one hydrometric level and a rating curve)
  hl_DA_matrix<-use_internal_data(list_name = "DA_matrices",dataset_name = "stream_stages")
  rc_DA_matrix<-use_internal_data(list_name = "DA_matrices",dataset_name = "rating_curves")
  str_DA_matrix<-use_internal_data(list_name = "DA_matrices",dataset_name = "obs_streamflows")
  
  `%>%` <- magrittr::`%>%`
  streamflows_HH<-data.frame(matrix(NA,0,3,dimnames=list(x=NULL,y=colnames(streamflows))))
  streamflows_GG<-data.frame(matrix(NA,0,3,dimnames=list(x=NULL,y=colnames(streamflows))))
  
  #Section with hydrometric levels and rating curves
  stations<-unique(hydrometric_levels[,1])[unique(hydrometric_levels[,1])%in%rownames(rc_DA_matrix)[which(apply(rc_DA_matrix,1,sum)>0)]]

  for(r in stations) 
  {
    years_rc<-as.numeric(names(which(rc_DA_matrix[rownames(rc_DA_matrix)%in%r,]==1)))
    hl_r<-hydrometric_levels[hydrometric_levels[,1]%in%r & lubridate::year(hydrometric_levels[,2])%in%years_rc,-1]
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
    str_r_HH[,1]<-r
      
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

    #Adding values outside the upper and lower limits
    str_r_GG<-rbind.data.frame(str_r_GG,str_eval_r[!str_eval_r[,1]%in%str_r_GG[,1],])
    str_r_GG<-str_r_GG %>% dplyr::arrange(Date)

    #Linear interpolation if less than 7 consecutive NA daily streamflows
    NA_streamflows_r <- split(which(is.na(str_r_GG[,2])),cumsum(c(1,diff(which(is.na(str_r_GG[,2])))!=1)))
    for(n in names(NA_streamflows_r))
    {
      i_start <- NA_streamflows_r[[n]][1]-1 #Fisrt day
      i_end <- NA_streamflows_r[[n]][length(NA_streamflows_r[[n]])]+1 #Last day
      
      if(i_end-i_start-1<7 && !is.na(str_r_GG[i_start,2]) && !is.na(str_r_GG[i_end,2]))
      {
        str_r_GG[i_start:i_end,2] <- approx(x=c(i_start,i_end),y=c(str_r_GG[i_start,2],str_r_GG[i_end,2]),method="linear",n=i_end-i_start+1,na.rm=TRUE)$y 
      }
    }
    if(length(which(is.na(str_r_GG[,2])))>0) str_r_GG<-str_r_GG[-which(is.na(str_r_GG[,2])),]
    str_r_GG<-cbind.data.frame(Name=r,str_r_GG)
    streamflows_GG<-rbind.data.frame(streamflows_GG,str_r_GG)
  } 

  streaflows_LT <-streamflows[is.na(streamflows[,2]),-2]
  rownames(streamflows)<-NULL
  
  # Store all dataframes inside a list
  .internal_env$streamflows <- list(HH=streamflows_HH,GG=streamflows_GG,LT=streaflows_LT)
  
  #How does streamflows matrix change 
  input_data(hydrometric_levels,rating_curves,use_internal_data("streamflows","GG"))
}

