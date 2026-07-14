#' Streamflow Quality Check
#'
#' Computes climate elasticity of streamflows and period-of-record flow duration
#' cyrve (POR-FDC) in order to remove years or sites that exhibit unphysical behavior
#'
#' @param daily_streamflows A named list of zoo time series. Each element must be a
#'     zoo object whose time index is of class `Date` or `POSIXct`. The list names
#'     represent site names
#'
#' @param annual_precipitations A named list of zoo time series. Each element must be
#'    a zoo object whose time index is of class `numeric`. The list names represent
#'     site names
#'
#' @param lower,upper values of type `numeric` representing typical bounds of the
#'    climate elasticity in the study area
#'
#' @return A list of two zoo time series: hourly and daily streamflows.
#'
#' @examples Assuming you have loaded daily streamflows and annual precipitation
#' lists of catchments located in a warm temperate climate
#' streamflows_check<-check_streamflow_quality(daily_streamflows,hourly_streamflows,
#'                                             annual_precipitations,annual_precipitation,
#'                                             lower=0,upper=5)
#'
#' Removal based on climate elasticity is performed automatically based on lower and
#' upper precipitation elasticity bounds
#'
#' Sites removal based on POR-FDC trend is performed manually by the user, based on
#' the computed number of inflection points of log(Q)=f(duration)
#'
#' #Access streamflow series
#' streamflows_check$hourly
#' streamflows_check$daily
#'

#Years in common
years_in_common<-function(x,y)
{
  if(inherits(time(x[[1]]),c("Date","POSIXct","POSIXt")) &&
     !inherits(time(y[[1]]),c("Date","POSIXct","POSIXt")))
    x_in_y <- Map(function(x, y) x %in% y, lapply(lapply(x,time),lubridate::year),
                  lapply(y,time))
  else
    x_in_y <- Map(function(x, y) x %in% y, lapply(x,time), lapply(y,time))

  x_kept <- Map(function(x, mask) x[mask], x, x_in_y)
  return(x_kept)
}


#Prompts for user removal descriptors
ask_which_descriptor <- function(prompt = "Which site/s?:")
{
  repeat
  {
    response <- trimws(readline(prompt))
    response <- strsplit(response,",")[[1]]
    if(length(which(response %in% names(flow_duration_curves)))>0) return(response)
    else cat("Invalid input. Please enter the proper descriptor/s name/s.\n")
  }
}
ask_yes_no <- function(prompt = "Do you want to delete one of these sites?[Y/N]:")
{
  repeat
  {
    response <- toupper(trimws(readline(prompt)))
    if(response %in% c("Y", "N")) return(response)
    else cat("Invalid input. Please enter 'Y' or 'N'.\n")
  }
}

#' @export
check_streamflow_quality<-function(daily_streamflows,hourly_streamflows,
                                   annual_precipitations,lower,upper)
{
  #Annual flows evaluation
  #removal of sites with less then 5 years
  count_years<-function(x) length(unique(lubridate::year(time(x))))
  daily_streamflows <- daily_streamflows[lapply(daily_streamflows,count_years)>=5]
  annual_flows_eval<-function(x) aggregate(x,lubridate::year,mean,na.rm=T)
  annual_flows<-lapply(daily_streamflows,annual_flows_eval)

  #Remove NA values
  annual_flows <- lapply(annual_flows,function(x) x[!is.na(x)])
  annual_precipitations <- lapply(annual_precipitations,function(x) x[!is.na(x)])

  #Sites in common
  annual_flows<-annual_flows[names(annual_precipitations)]
  annual_precipitations<-annual_precipitations[names(annual_flows)]
  annual_flows<-years_in_common(annual_flows,annual_precipitations)
  annual_precipitations<-years_in_common(annual_precipitations,annual_flows)


  #Check #1: Precipitation elasticity of streamflows
  loop<-T;lengths_sites<-0#initial conditions
  while(loop)
  {
    attempt<-prec_elasticity_eval(annual_flows,
                                  annual_precipitations,
                                  lower,upper)
    if(!identical(lengths_sites,lengths(attempt$flows)))
    {
      annual_flows<-attempt$flows;annual_flows<-annual_flows[lengths(annual_flows)>=5]

      annual_precipitations<-attempt$precipitations;annual_precipitations<-
        annual_precipitations[lengths(annual_precipitations)>=5]

      lengths_sites<-lengths(annual_flows)
    }else loop<-F
  }
  elasticities<-attempt$elasticities
  cat("Precipitation elasticities of streamflows:\n")
  print(lapply(elasticities,function(x) round(x,digits = 2)))

  #Mask hourly and daily streamflow series according to precipitation elasticity check
  #Sites in common
  daily_streamflows<-daily_streamflows[which(names(daily_streamflows)%in%
                                         names(elasticities))]
  hourly_streamflows<-hourly_streamflows[which(names(hourly_streamflows)%in%
                                          names(elasticities))]
  #Years in common
  daily_streamflows<-years_in_common(daily_streamflows,elasticities)
  hourly_streamflows<-years_in_common(hourly_streamflows,elasticities)


  #Check #2: number of inflection points in the logarithmic POR-FDC
  daily_streamflows_log<-lapply(daily_streamflows,log)
  daily_streamflows_log<-Map(function(x, mask) x[mask],
                        daily_streamflows_log,lapply(daily_streamflows_log,is.finite))

  #Uniforming PORFDCs: same vector of durations
  n<-max(lengths(daily_streamflows_log))#number of points
  p_inf<-1/(1+min(lengths(daily_streamflows_log)))#lower bound of durations
  p_sup<-min(lengths(daily_streamflows_log))/
          (1+min(lengths(daily_streamflows_log)))#upper bound for durations
  p_0  <-1:n/(n+1)#lower bound for durations
  p_0<-p_0[p_0>=p_inf&p_0<=p_sup]

  #Evaluation of period of record flow duration curve
  flow_duration_curves<-lapply(daily_streamflows_log,porfdc_eval)
  #porfdc computed at same durations
  n_inflection_points<-c()
  for(s in names(flow_duration_curves))
  {
    #change of sign in the second derivative
    porfdc_unif<-approx(x = time(flow_duration_curves[[s]]),#unifromed
                       y = flow_duration_curves[[s]],
                       xout = p_0 ,method="linear")$y
    dy_dx <- diff(porfdc_unif) / diff(p_0)# Compute first derivative
    d2y_d2x <- diff(dy_dx) / diff(p_0[-1])# Compute second derivative
    n_inflection_points<-append(n_inflection_points,
                                length(which(diff(sign(d2y_d2x))!=0)))#Number of inflection points
  }
  names(n_inflection_points)<-names(flow_duration_curves)

  cat("Number of inflection points: \n")
  print(n_inflection_points)

  for(s in names(flow_duration_curves))#plot
  {
    plot(flow_duration_curves[[s]],main=s,lwd=2,
         xlab="duration (-)",ylab=expression(paste("streamflow (",m^3,"/s)",sep="")))
    grid()
  }
  user_input_1<-ask_yes_no() #to user: sites removal?
  if(user_input=="Y")
  {
    user_input_2<-ask_which_descriptor()#to user: which site/s?
    daily_streamflows<- daily_streamflows[!names(daily_streamflows)%in%user_input_2]
    hourly_streamflows<- hourly_streamflows[!names(hourly_streamflows)%in%user_input_2]
  }

  # Store all dataframes inside a list
  output <- list(hourly=hourly_streamflows,
                 daily=daily_streamflows)
  return(output)
}
