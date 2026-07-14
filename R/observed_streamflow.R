#' Observed streamflow estimation
#'
#' Computes observed hourly, daily, and long-term streamflow series.
#'
#' @param stream_stages A named list of zoo time series (one per site). Each
#'    element must be a zoo object whose time index is of class `Date` or
#'    `POSIXct`. The list names represent site names
#'
#' @param rating_curves A six column-data frame: station name (`character`),
#'   start date of applicability and end date of applicability (`Date`or`POSIXct`),
#'   minimum and maximum applicable stage (`numeric`), rating curve equation
#'  (`character`)
#'
#' @param daily_streamflows A named list of zoo time series (one per site). Each
#'    element must be a zoo object whose time index is of class `Date`.
#'    The list names represent site identifiers
#'
#' @return A list of two data frames, one for each time step: hourly (`hourly`),
#'   daily (`daily`)
#'
#' @examples
#' # Immagine you have stream_stages, rating_curves and streamflows dataframes loaded
#' obs_streamflow<-observed_streamflow(stream_stages,rating_curves,streamflows)
#'
#' obs_streamflows$hourly
#' obs_streamflows$daily
#'
#' @export
observed_streamflow<-function(stream_stages,rating_curves,daily_streamflows)
{
  #Streamflow series generation (sites that has at least one stream stage and a rating curve)
  hourly_streamflows<- list()

  #Site with stream stages and rating curves
  stations<-names(stream_stages)[which(lengths(sapply(
      names(stream_stages),function(x) grep(x,rating_curves[,1],ignore.case = T)))>0)]

  for(s in stations)
  {
    years<-names(which(rating_avail_matrix[rownames(rating_avail_matrix)%in%s,]==1))
    stages_s<-stream_stages[[s]]
    stages_s<-stages_s[lubridate::year(time(stages_s))%in%years]

    rating_s<-rating_curves[grep(s,rating_curves[,1],ignore.case = T),]
    flows_eval_s<-rating_curve_eval(stages_s,rating_s)
    flows_eval_s<-flows_eval_s[!is.na(flows_eval_s)]#remove NAs

    #Append hourly streamflows for site s
    hourly_streamflows<-append(hourly_streamflows,list(flows_eval_s))
    names(hourly_streamflows)[which(stations%in%s)]<-s

    flows_eval_s <- aggregate(flows_eval_s,as.Date,mean)
    if(s%in%names(daily_streamflows))#sites with daily streamflows already
    {
      #Filling gaps
      daily_streamflows[[s]][is.na(daily_streamflows[[s]]) & time(daily_streamflows[[s]])%in%time(flows_eval_s)]<-
        flows_eval_s[time(flows_eval_s)%in%time(daily_streamflows[[s]][is.na(daily_streamflows[[s]])])]

      #Adding values outside the streamflows time-span
      daily_streamflows[[s]]<-rbind(daily_streamflows[[s]],flows_eval_s[!time(flows_eval_s)%in%
                                                                   time(daily_streamflows[[s]])])
    }else#Sites with no daily streamflows
    {
      daily_streamflows<-append(daily_streamflows,list(flows_eval_s))
      names(daily_streamflows)[length(daily_streamflows)]<-s
    }

    #Linear interpolation if less than 15 consecutive NA daily streamflows
    na_flows_s <- split(which(is.na(daily_streamflows[[s]])),
                              cumsum(c(1,diff(which(is.na(daily_streamflows[[s]])))!=1)))
    for(n in names(na_flows_s))
    {
      if(length(na_flows_s[[n]])>0)
      {

        i_start <- na_flows_s[[n]][1]-1 #First day
        i_end <- na_flows_s[[n]][length(na_flows_s[[n]])]+1 #Last day

        if(i_end-i_start-1<15 && !is.na(daily_streamflows[[s]][i_start]) && !is.na(daily_streamflows[[s]][i_end]))
        {
          daily_streamflows[[s]][i_start:i_end] <- approx(x=c(i_start,i_end),y=c(daily_streamflows[[s]][i_start],daily_streamflows[[s]][i_end]),
                                                          method="linear",n=i_end-i_start+1,na.rm=TRUE)$y
        }
      }
    }
    year_to_remove <- unique(lubridate::year(time(daily_streamflows[[s]][is.na(daily_streamflows[[s]])])))#years with at least one NA
    if(length(year_to_remove)>0) daily_streamflows[[s]]<-daily_streamflows[[s]][-which(lubridate::year(time(daily_streamflows[[s]]))%in%year_to_remove)]
  }

  #Sort in alphabetical order
  daily_streamflows<-daily_streamflows[order(names(daily_streamflows))]
  hourly_streamflows<-hourly_streamflows[order(names(hourly_streamflows))]

  # Store all lists inside a list
  output <- list(hourly=hourly_streamflows,daily=daily_streamflows)

  return(output)
}
