#' Synthetic hourly streamflow series of the target site
#'
#' Downscales synthetic daily streamflow series to hourly time-step based on linear
#'  (normal and low flows) and polynomial (high flows) interpolation by means of
#'  regionalized parameters of the hydrograph
#'
#' @param syn_daily_streamflows A zoo time series of the regional daily streamflow series
#' @param porfdc A zoo object representing the regional prediction of the period
#'   of record flow duration curve
#' @param target A `character` with the target site name (by default set as )
#' @param peak_pars A data frame with one row for each regionalized parameter
#'   of the at-site hourly to daily peak ratio relationship. Columns give the
#'   parameter name (`character`), the formula used for regionalization (`character`),
#'   and the formula's coefficients (`numeric`).
#' @param shape_pars A data frame with one row for each regionalized parameter
#'   of the at-site shape ratio relationships. Columns give the parameter name
#'   (`character`), the formula used for regionalization (`character`), and
#'   the formula's coefficients (`numeric`).
#' @param D_target A `numeric` value giving the duration (in hours) of
#'   high-flow events at the target section. Defaults to `25`.
#' @param d_threshold A `numeric` value giving the duration threshold for
#'   high-flow events, as used in [regional_porfdc()]. Defaults to `0.2`.
#' @param mae_threshold A `numeric` representing the mean absolute error thresholds
#'   of high-flow hydrograph and high-flow volumes. Defaults to  `0.2`.
#' @param max_iter A `numeric` indicating the number of maximum iterations in polynomial
#'   interpolation. Defaults to  `20`.
#'
#' @return A zoo time series of synthetic hourly streamflows at the target site
#'
#' @examples
#' #Assuming you have syn_daily_streamflows, porfdc, target, peak_pars and shape_pars
#' #loaded
#' synthetic_hourly_streamflow(syn_daily_streamflows,porfdc,target, peak_pars,shape_pars,
#'                             D_target,d_threshold)
#'
#' #This computes a synthetic hourly streamflow series at the target site using
#' #by downscaling synthetic daily series according to peak-ratio and shape-ratio
#' # relationships in high-flow conditions.

# Snthetic Design Hydrograph (from Maione et al, Int. J. River Basin Management, 2003)
synthetic_design_hydrograph<-function(j)
{
  if(j==1)#rising limb
  {
    t_sdh<--rD*Durations#sdh time-steps
    q_sdh<-(drD_dD*Durations*QD+rD*QD+rD*Durations*dQD_dD)/(drD_dD*Durations+rD)#sdh hourly streamflows
  }

  if(j==3)#falling limb
  {
    t_sdh<-Durations*(1-rD)#sdh time-steps
    q_sdh<-(-drD_dD*Durations*QD+(1-rD)*QD+(1-rD)*Durations*dQD_dD)/(-drD_dD*Durations+1-rD)#sdh hourly streamflows
  }
  t_sdh<-t_sdh+max_pos#time-steps
  output<-data.frame(x=t_sdh,y=q_sdh)
  return(output)
}

# Constrained polynomial y = a*x^2 + b*x + c for representing rising and falling limbs, by imposing:
#   - first and last values equal to daily streamflow of previous and following days
#   - passage through the peak
#   - minimize the mean absolute error with respect to the synthetic design hydrograph
#   - minimize the absolute error with respect to daily flood volume
constrained_polynomial <- function(c)
{
  A <- matrix(c(x[j]^2,x[j],#Coefficient matrix
                x[which.max(y)]^2,x[which.max(y)]),
              2, 2, byrow=TRUE)

  colnames(A) <- c("a","b")#parameters

  b <- c(y[j]-c,y[which.max(y)]-c)#Vector of known terms

  pars <- c(solve(A, b),c) #parameters vector

  if (j==1)#rising limb
  {
    t <- seq(from=1,to=max_pos,by=1)#time-steps
    q <- pars[1]*t^2 + pars[2]*t + pars[3]#hourly streamflows of the high-flow conditions
  }

  if (j==3)#falling limb
  {
    t <- seq(from=max_pos,to=49,by=1)#time-steps
    q <- pars[1]*t^2 + pars[2]*t + pars[3]#hourly streamflows of the high-flow conditions
  }

  if(result=="calibration")#parameters calibration: returns a weighted average between differences with daily flood volumes and SDH
  {
    sdh<-synthetic_design_hydrograph(j)#computation of the synthetic design hydrograph
    hourly_volume <- pracma::trapz(t,q)#hourly volumes of calibrated hydrograph at the current step

    deltaQ<-sum(abs((pars[1]*sdh$x^2 + pars[2]*sdh$x + pars[3])-sdh$y)/sdh$y)/length(sdh$y)#streamflow mean absolute error
    deltaV <- abs(hourly_volume/daily_volume-1)#volume absolute error

    delta<-(deltaQ+deltaV)*0.5#stremflow and volume MAEs equally weighted
    return(delta)
  }
  if(result=="generation") return(q) #generation of hourly streamflow series of the event with calib parameters
}

# Convert a date to a datetime set at 12:00:00
noon <- function(x) lubridate::ymd_hms(paste0(x," 12:00:00"))

synthetic_hourly_streamflow<-function(syn_daily_streamflows,porfdc,target="Marecchia@RiminiSS16",
                                      peak_pars,shape_pars,D_target=25,d_threshold=0.2,mae_threshold=0.2,max_iter=20)
{
  # Grouping consecutive years of the regional daily streamflow series
  group_cons_years<-split(unique(lubridate::year(syn_daily_streamflows)),
                          cumsum(c(1,diff(unique(lubridate::year(syn_daily_streamflows)))!=1)))

  # Time of concentration for target site
  descriptors_target<-basin_descriptors[names(tc)%in%target,]
  tc_target<-time_of_concentration(descriptors_target)

  # Peak-ratio and shape hydrograph parameters for the target site
  x<-as.numeric(tc_target) #where to evaluate the formulas
  peak_pars_target <- c()#Hourly-to-daily peak ratio parameters
  for(i in 1:length(peak_pars[,1]))
  {
    if(is.na(peak_pars[i,5])) {P1<-peak_pars[i,3];P2<-peak_pars[i,4]} else {P1<-peak_pars[i,3];P2<-peak_pars[i,4];P3<-peak_pars[i,5]}
    peak_pars_target <- append(peak_pars_target,eval(parse(text=substring(peak_pars[i,2],first=3,last=nchar(peak_pars[i,2])))))
  }
  shape_pars_target<-c() #Peak-duration curve parameters
  for(i in 1:length(shape_pars[-3,1]))
  {
    q<-shape_pars[i,3];m<-shape_pars[i,4]
    shape_pars_target <- append(shape_pars_target,eval(parse(text=substring(shape_pars[i,2],first=3,last=nchar(shape_pars[i,2])))))
  }
  shape_pars_target <- append(shape_pars_target,shape_pars[3,3]+shape_pars[3,4]*x)# Reduction ratio parameter

  Durations<-c(0.001,1:D_target)#high-flows event durations
  norm_Durations<-Durations/tc_target#normalized Durations (by time of concentration)

  # Compute the variables that define the synthetic design hydrograph (Maione et al.,2003)
  rD<-.5*exp(-shape_pars_target[1]*norm_Durations/(1+shape_pars_target[2]*norm_Durations))#Peak-duration curve
  epsD<-shape_pars_target[3]/norm_Durations*(1-exp(-norm_Durations/shape_pars_target[3]))#Reduction ratio
  drD_dD<--((shape_pars_target[1]/tc_target)/(1+shape_pars_target[2]*norm_Durations)^2)*rD#Derivative of peak-duration (with respect to event duration)
  depsD_dD<-1/norm_Durations*(-shape_pars_target[3]/(tc_target*norm_Durations)*(1-exp(-norm_Durations/shape_pars_target[3]))+
                                (1/tc_target*exp(-norm_Durations/shape_pars_target[3])))#Derivative of the reduction ratio (with respect to event duration)
  for(k in 1:length(group_cons_years))#cycle on each group of consecutive years
  {
    dailyflows_k<-syn_daily_streamflows[lubridate::year(syn_daily_streamflows)%in%group_cons_years[[k]],] #daily streamflow in group k

    flows_d <- data.frame(date=zoo::index(dailyflows_k),
                          value=zoo::coredata(dailyflows_k),
                          p=NA)

    #Daily streamflow durations
    for(q in unique(reg_porfdc)) flows_d[flows_d[,2]%in%q,3] <- mean(zoo::index(reg_porfdc[reg_porfdc%in%q]))

    #Downscaling to hourly streamflows: linear interpolation, daily streamflows at noon
    z    <- zoo(coredata(dailyflows_k), noon(zoo::index(dailyflows_k)))
    datetime <- seq(from = start(z), to = end(z), by = "hour")
    hourlyflows_k   <- na.approx(z, xout = datetime)

    #High-flow conditions (d<=d_threshold and daily streamflow higher than contiguous days): hourly peak fixed at 12:00:00
    #and two 2nd order polynomials describe the rising and the falling limb (y=a*t^2+b*t+c); its three parameters
    #are estimated by constraining the passage through the hourly peak, the passage through the daily streamflows of the previous and
    # the following day, conservation of the flow volumes and of the shape of the hydrograph

    #High flows detection (daily streamflow higher than contiguous days and with duration lower or equal to d_threshold)
    index_flood<-c()
    for(i in 2:(length(flows_d[,1])-1))
    {
      if(flows_d[i,1]-flows_d[i-1,1]==1 && flows_d[i+1,1]-flows_d[i,1]==1)
      {
        if(flows_d[i,2]>=flows_d[i-1,2] &&
           flows_d[i,2]>=flows_d[i+1,2] &&
           flows_d[i,3]<=d_threshold) index_flood<-append(index_flood,i)
      }
    }

    for(i in index_flood)#cycle on high-flow day
    {
      Qmax<-(peak_pars_target[1]+peak_pars_target[2]*exp(-peak_pars_target[3]*flows_d[i,3]))*flows_d[i,2]#hourly peak
      QD<-epsD*Qmax #maximum average discharge at varying flood durations
      dQD_dD<-depsD_dD*Qmax #first derivative of QD

      max_pos<-25#Position of the flood peak (at 12:00:00 of the high-flow day)
      x<-c(1,max_pos,49)#hours in which hourly streamflow is known
      y<-c(flows_d[i-1,2],Qmax,flows_d[i+1,2])#known hourly streamflows

      #Polynomial paramters' estimates, for the rising (1st>25th hour) and for the falling (25th>49th hour) of the flood event
      q_complete<-c()#Initializing the variable of the calibrated high-flow hydrograph
      lower_bound <- -10#lower end point of the interval where to search
      upper_bound <- 10 #upper end point of the interval where to search

      for(j in c(1,3))#cycle on the hydrograph limbs  (1: rising, 3: falling, referring to the known streamflow in vector y)
      {
        if(j==1) daily_volume<-12*flows_d[i-1,2]+(max_pos-13)*flows_d[i,2]#Volumi idrogramma giornaliero
        if(j==3) daily_volume<-12*flows_d[i+1,2]+(37-max_pos)*flows_d[i,2]#Volumi idrogramma giornaliero

        result<-"calibration";mae_above_threshold <- TRUE;cont <- 0#hydrograph calibration
        while(mae_above_threshold && cont < max_iter)
        {
          cont<-cont+1
          delta <- optimize(constrained_polynomial,lower=lower_bound,upper=upper_bound)#

          if (delta$objective < mae_threshold) {
            mae_above_threshold <- FALSE
          } else # lower and upper bounds expanded adaptively
          {
            par <- delta$minimum

            if (abs(par - upper_bound) <= abs(par - lower_bound)) {
              lower_bound <- par
              upper_bound <- upper_bound + abs(upper_bound) * 0.5
            } else {
              upper_bound <- par
              lower_bound <- lower_bound - abs(lower_bound) * 0.5
            }
          }
        }
        result<-"generation";q_half<-constrained_polynomial(par)

        #Ensure monotonicity and overwrite linearly interpolated daily streamflows
        if(j==1)#rising limb
        {
          index_monotonic<- which(!q_half < cummax(q_half))
          window(hourlyflows_k,start=noon(flows_d[i-1,1]),end=noon(flows_d[i,1]))[index_monotonic]<-q_half
        }
        if(j==3)#falling limb
        {
          index_monotonic<- which(!q_half > cummin(q_half))
          window(hourlyflows_k,start=noon(flows_d[i,1]),end=noon(flows_d[i+1,1]))[index_monotonic]<-q_half
        }
      }
    }
    if(k==1) syn_hourly_streamflows<-hourlyflows_k#First group of consecutive years
    if(k>1) syn_hourly_streamflows<-c(syn_hourly_streamflows,hourlyflows_k)#Second groups of consecutive years
  }

  #regional hourly streamflow object
  return(syn_hourly_streamflows)
}
