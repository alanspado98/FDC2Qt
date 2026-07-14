#' Rating Curve Evaluation
#'
#' Generates at-site streamflow estimates by applying the stage–discharge
#' relationship to observed stream stages.
#'
#' @param stream_stages A zoo time series. Time index is of class `Date` or
#'    `POSIXct`. The list names represent site names
#'
#' @param rating_curves A six column-data frame: station name (`character`),
#'   start date of applicability and end date of applicability (`Date`or`POSIXct`),
#'   minimum and maximum applicable stage (`numeric`), rating curve equation
#'  (`character`)
#'
#' @return A zoo time series of the computed streamflows
#'
#' @export

rating_curve_eval<-function(stream_stages,rating_curves)
{
  #Streamflow series inizialization
  streamflows <- zoo::zoo(rep(NA,length(stream_stages)), order.by = time(stream_stages))

  stream_stages<-stream_stages[!is.na(stream_stages[,2])]#Remove NA stream stages

  indeces<-c(1,which(diff.POSIXt(rating_curves[,2])!=0)+1)

  for(i in indeces)
  {
    stages_i<-stream_stages[time(stream_stages)>=rating_curves[i,2] &
                          time(stream_stages)<rating_curves[i,3],]
    if(length(stages_i)>0)
    {
      rating_i<-rating_curves[rating_curves[,2]%in%rating_curves[i,2],]
      for(j in 1:length(rating_i[,1])) #stages intervals of the rating curves
      {
        stages_ij<-stages_i[stages_i>=rating_i[j,4] & stages_i<rating_i[j,5],]
        if(length(stages_ij)>0)
        {
          X<-stages_ij
          streamflows[time(streamflows)%in%time(stages_ij)]<-eval(parse(text=rating_i[j,6]))
        }
      }
    }
  }
  return(streamflows)
}
