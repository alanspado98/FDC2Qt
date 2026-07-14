#' Synthetic daily streamflow series of the target site
#'
#' Detects donor site according to user purpose and computes daily streamflow
#' series in the target site through a nonlinear interpolation of PORFDCs
#'
#' @param daily_streamflows A named list of zoo time series. Each element must be
#'   a zoo object whose time index is of class `Date` or `POSIXct`. The list names
#'   represent site names
#' @param porfdc A data frame of the regional period-of-record flow duration curve
#'   of the traget site: duration in the first column (`numeric`) and discharge
#'   in the second (`numeric`).
#' @param maf numeric. Regional mean annual flow of the target site
#' @param target Target site (`character`)
#' @param generation logical. Defines the strategy for selecting the donor sites.
#'   If TRUE the longest record in the region is selected (generation purpose),
#'   if FALSE the most-time correlated record with the target sites is selected
#'   (prediction purpose)
#'
#' @return A zoo time series of synthetic daily streamflows
#'
#' @example
#' # Assuming you have daily_streamflows, porfdc, maf and target loaded:
#' synthetic_daily_streamflow(daily_streamflows,porfdc,maf,target,generation=TRUE)
#'
#' # This computes a synthetic daily streamflow series at the target site using
#' # the longest observed daily discharge series in the region
#' (`basin_descriptors`).

#Days in common
days_in_common<-function(x,y)
{
  y_time<-time(y)
  lapply(x,function(z) z[time(z)%in%y_time])
}


#' @export
synthetic_daily_streamflow<-function(daily_streamflows,porfdc,maf,target,generation=TRUE)
{
  #Donator site, choice depens on the aim:
  #1) hydrologically consistent long streamflow series -> longest series within the region
  #2) reconstruct historical streamflow series (target site must be, at least partially, gauged!)
  # -> Perason correlation coefficient of daily streamflow series

  #1)
  if(generation) donor_section<-names(which.max(lengths(daily_streamflows)))

  #2)
  if(!generation && !is.na(match(target,names(daily_streamflows))))
  {
    flows_notarget <- daily_streamflows[-which(names(daily_streamflows)%in%target)]
    flows_target <- daily_streamflows[names(daily_streamflows)%in%target][[1]]

    flows_notarget<-days_in_common(flows_notarget,flows_target)#days in common
    index_s <- which(lengths(flows_notarget)==0)
    if(length(index_s)>0) flows_notarget<-flows_notarget[-index_s]

    R2<-rep(NA,length(flows_notarget));names(R2)<-names(flows_notarget)
    for(s in names(R2))
    {
      flows_notarget_s<-flows_notarget[names(flows_notarget)%in%s][[1]]
      flows_target_s<-flows_target[time(flows_target)%in%time(flows_notarget_s)]
      R2[names(R2)%in%s]<-cor(as.numeric(flows_notarget_s),
                              as.numeric(flows_target_s))^2
    }
   donor_section<-names(which.max(R2))
  }

  #PORFDCdonor site
  obs_daily_streamflows<-daily_streamflows[names(daily_streamflows)%in%donor_section][[1]]
  Q <- sort(unique(as.numeric(obs_daily_streamflows)),decreasing=T)
  donor_porfdc <- data.frame(p=c(1:length(Q)/(length(Q)+1)),Q)

  #Resampling regional PORFDC according to PORFDC of donor site
  reg_porfdc<-porfdc
  reg_porfdc[,2]<-reg_porfdc[,2]*maf
  reg_porfdc <- data.frame(p=donor_porfdc$p,Q=approx(x=reg_porfdc$p,y=reg_porfdc$Q,xout=donor_porfdc$p,method="linear")$y)

  #Synthetic hydrograph for the target site
  #Dataframe intialization
  syn_daily_streamflows <- zoo::zoo(vector("numeric",length(obs_daily_streamflows)))
  zoo::index(syn_daily_streamflows) <- zoo::index(obs_daily_streamflows)

  for(i in 1:length(syn_daily_streamflows))
  {
    index<-which(donor_porfdc[,2]==as.numeric(obs_daily_streamflows[i]))#duration in the donor site PORFDC corresponding to i-th daily streamflow
    syn_daily_streamflows[i]<-reg_porfdc[index,2]#streamflow in regional PORFDC corrisponding to the duration of the i-th daily streamflow
  }

  # regional daily streamflows dataframe
  return(syn_daily_streamflows)
}
