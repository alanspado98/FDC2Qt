#' Synthetic daily streamflow series of the target site
#'
#' Detects the donor site according to user purpose and computes daily streamflow
#' series in the target site through a nonlinear interpolation of PORFDCs
#'
#' @param daily_streamflows A named list of zoo time series. Each element must be
#'   a zoo object whose time index is of class `Date` or `POSIXct`. The list names
#'   represent site names
#' @param porfdc A data frame of the regional period-of-record flow duration curve
#'   of the traget site: duration in the first column (`numeric`) and discharge
#'   in the second (`numeric`).
#' @param target Target site (`character`)
#' @param n Number of donor sites (`numeric`). Defaults to `1`.
#' @param generation logical. Defines the strategy for selecting the donor sites.
#'   If TRUE the longest record in the region is selected (generation purpose),
#'   if FALSE the most-time correlated record with the target sites is selected
#'   (prediction purpose)
#'
#' @return A zoo time series of synthetic daily streamflows and a dataframe of
#'   the regional prediction of the period-of-record flow duration curve at the
#'   target site
#'
#' @examples
#' # Assuming you have daily_streamflows, porfdc, maf and target loaded:
#' synthetic_daily_streamflow(daily_streamflows,porfdc,target,n=1,generation=TRUE)
#'
#' # This computes a synthetic daily streamflow series and the regional PORFDC at
#' # the target site using the longest observed daily discharge series in the region

#Days in common
days_in_common<-function(x,y)
{
  y_time<-time(y)
  lapply(x,function(z) z[time(z)%in%y_time])
}


#' @export
synthetic_daily_streamflow<-function(daily_streamflows,porfdc,target="Marecchia@RiminiSS16",n=1,generation=TRUE)
{
  #Donator site, choice depens on the aim:
  #1) hydrologically consistent long streamflow series -> longest series within the region
  #2) reconstruct historical streamflow series (target site must be, at least partially, gauged!)
  # -> Spearman correlation coefficient of daily streamflow series

  #1)
  if(generation)
  {
    l<-lengths(daily_streamflows)
    potential_donors<-names(l)[order(l,decreasing=T)]
  }

  #2)
  if(!generation && !is.na(match(target,names(daily_streamflows))))
  {
    #target streamflows
    flows_target <- daily_streamflows[names(daily_streamflows)%in%target][[1]]
    #all the other daily streamflows
    flows_notarget <- daily_streamflows[-which(names(daily_streamflows)%in%target)]

    flows_notarget<-days_in_common(flows_notarget,flows_target)#days in common
    index_s <- which(lengths(flows_notarget)==0)#remove sites with no days in common with the target
    if(length(index_s)>0) flows_notarget<-flows_notarget[-index_s]

    rho<-rep(NA,length(flows_notarget));names(rho)<-names(flows_notarget)#Spearman coefficient
    for(s in names(rho))
    {
      flows_notarget_s<-flows_notarget[names(flows_notarget)%in%s][[1]]
      flows_target_s<-flows_target[time(flows_target)%in%time(flows_notarget_s)]
      rho[names(rho)%in%s]<-cor(as.numeric(flows_notarget_s),
                                as.numeric(flows_target_s),
                                method = "spearman",
                                use="complete.obs")
    }
    potential_donors<-names(rho)[order(rho,decreasing=T)]
  }

  #Dimensional period-of-record flow duration curve at the target site
  reg_porfdc<-porfdc
  reg_porfdc[,2]<-reg_porfdc[,2]*maf

  #PORFDCdonor site
  donors<-potential_donors[1:n]
  syn_flows<-list()
  for(s in donors)#cycle on the number of sites chosen as donors
  {
    flows_s  <- daily_streamflows[names(daily_streamflows)%in%s][[1]]
    porfdc_s <- porfdc_eval(flows_s)

    #Resampling regional PORFDC according to PORFDC of donor site
    reg_porfdc_s <- data.frame(p=zoo::index(porfdc_s),
                               Q=approx(x=zoo::index(reg_porfdc),y=zoo::coredata(reg_porfdc),
                                        xout=zoo::index(porfdc_s),
                                        method="linear",rule=2)$y)

    #Synthetic hydrograph for the target site
    #zoo series intialization
    syn_flows_s <- zoo::zoo(vector("numeric",length(flows_s)))
    zoo::index(syn_flows_s) <- zoo::index(flows_s)

    #Non-linear interpolation based on PORFDCs
    for(q in unique(flows_s))
    {
      index<-zoo::index(porfdc_s[porfdc_s%in%q])#duration in the donor site PORFDC corresponding to i-th daily streamflow
      syn_flows_s[flows_s%in%q]<-mean(reg_porfdc_s[reg_porfdc_s$p%in%index,2])#streamflow in regional PORFDC corrisponding to the duration of the i-th daily streamflow
    }
    syn_flows<-append(syn_flows,list(syn_flows_s))#append synthetic daily streamflow for site s
  }
  syn_flows <- do.call(merge,syn_flows)
  syn_daily_streamflows <- zoo::zoo(rowMeans(syn_flows, na.rm = TRUE),index(syn_flows))
  output<-list(syn_daily=syn_flows_s,
               porfdc=reg_porfdc)
  # regional daily streamflow zoo object
  return(output)
}
