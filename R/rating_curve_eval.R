#' Rating Curve Evaluation
#' 
#' @param hydro_levels 
#'
#' @param rc_name 
#'
#' @export

rating_curve_eval<-function(hydro_levels,rc_name)
{
  streamflows<-hydro_levels
  streamflows[,2]<-NA
  
  hydro_levels<-hydro_levels[!is.na(hydro_levels[,2]),]#Remove NA streamflows
  
  datetime_indexes<-c(1,which(diff.POSIXt(rc_name$StartDate)!=0)+1)
  
  for(datetime_index in datetime_indexes)
  {
    hl_t<-hydro_levels[hydro_levels[,1]>=rc_name[datetime_index,2] & hydro_levels[,1]<rc_name[datetime_index,3],]
    rc_name_datetime<-rc_name[rc_name[,2]%in%rc_name[datetime_index,2],]
    
    if(length(hl_t[,1])>0)
    {
      for(value_index in 1:length(rc_name_datetime[,1])) #Hydrometric level interval
      {
        hl_t_l<-hl_t[hl_t[,2]>=rc_name_datetime[value_index,4] & hl_t[,2]<rc_name_datetime[value_index,5],]
        if(length(hl_t_l[,1])>0)
        {
          X<-hl_t_l[,2]
          streamflows[streamflows[,1]%in%hl_t_l[,1],2]<-eval(parse(text=rc_name_datetime[value_index,6]))
        }
      }
    }
  }
  return(streamflows)
}

  