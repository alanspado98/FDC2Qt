#' Precipitation elasticity evaluation
#'
#' Compute climate elasticity of streamflows, given mean annual flow and total
#' annual precipitation series
#'
#' @param annual_flows A named list of zoo time series. Each element must be a
#'     zoo object whose time index is of class `Date` or `POSIXct`. The list names
#'     represent site names
#' @param annual_precipitations A named list of zoo time series. Each element must be
#'    a zoo object whose time index is of class `numeric`. The list names represent
#'     site names
#' @param lower,upper values of type `numeric` representing typical bounds of the
#'    climate elasticity in the study area
#'
#' @return A list of three lists containing computed zoo time series according to lower
#'    and upper bounds of climate elasticity: annual flows, annual precipitations,
#'    elasticites
#'
#' @export
prec_elasticity_eval<-function(annual_flows,annual_precipitations,lower,upper)
{
  mean_annual_flows<-unlist(lapply(annual_flows,mean))#mean annual flows computation
  mean_annual_precipitations<-unlist(lapply(annual_precipitations,mean))#mean annual precipitation computation

  #Deviations computation
  deviations_eval<-function(values,mean_values) Map(function(x,y) x/y-1,values,mean_values)

  annual_flows_dev<-deviations_eval(annual_flows,mean_annual_flows)
  annual_precipitations_dev<-deviations_eval(annual_precipitations,mean_annual_precipitations)

  #Precipitation deviations >= 10% in order to not have higher elasticity values
  criterion_mask<-lapply(annual_precipitations_dev,function(x) abs(x)>=.1)
  values_kept<-function(x,mask) Map(function(x, mask) x[mask],x,mask)
  annual_precipitations_dev <- values_kept(annual_precipitations_dev,criterion_mask)

  elasticities<-Map(function(flows,precipitations) flows/precipitations, annual_flows_dev,annual_precipitations_dev)
  bounds_mask<-lapply(elasticities,function(x) x>=lower & x<=upper)

  annual_flows<-values_kept(annual_flows,bounds_mask)
  annual_precipitations<-values_kept(annual_precipitations,bounds_mask)
  elasticities<-values_kept(elasticities,bounds_mask)

  output<-list(flows=annual_flows,
               precipitations=annual_precipitations,
               elasticities=elasticities)

  return(output)
}
