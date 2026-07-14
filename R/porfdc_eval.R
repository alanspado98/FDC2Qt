#' Flow duration curves construction
#'
#' Compute the duration-discharge relationship for the entire streamflow series
#' (period-of-record), according to Weibull plotting position
#'
#' @param streamflows A zoo time series
#'
#' @return A zoo object representing the period of record flow duration curve:
#'       sorted discharges in descending order and associated duration as time
#'       index (`numeric`)
#' @export
porfdc_eval<-function(streamflows)
{
  streamflows_core<-zoo::coredata(streamflows)
  streamflows_sort<-streamflows_core[order(streamflows,decreasing=T,na.last=NA)]
  porfdc<-zoo::zoo(streamflows_sort)
  zoo::index(porfdc)<-1:length(streamflows_sort)/(length(streamflows_sort)+1)
  return(porfdc)
}
