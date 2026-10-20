#' Hourly to daily peak ratio parameters
#'
#' Computes ratios of daily streamflow under maximum hourly streamflow as a function
#' of the associated daily duration, fit an exponential model and regionalizes its
#'  parameters according to time of concentration
#'
#' @param hourly_streamflows A named list of zoo time series (one per site). Each
#'    element must be a zoo object whose time index is of class `POSIXct`.
#'    The list names represent site identifiers.
#' @param basin_descriptors A data frame whose rows represents sites, columns
#'   multiple columns: station name in the first (`character`) and morpho‑climatic
#'   descriptors in the remaining columns (`numeric`).
#' @param d_threshold A numeric of duration threshold for high-flow events
#'   (evaluated in [(regional_porfdc)]).
#'
#' @return A data frame of with the exponential-model parameters as a function of time of concentration

#Exponential model function: y=a+b*exp(-c*x) for HH/GG = f(duration) fit
exponential_model_1 <- function(c)
{
  b <- (y[1]-y[length(x)])/(exp(-c*x[1])-exp(-c*x[length(x)]))#a,b model parameters
  a <- y[1] - b*exp(-c*x[1])
  est <- a + b*exp(-c*x)#Model estimates
  sumsq<-sum((y - est)^2) #Sum of squared residuals
  return(sumsq)
}

#Exponential models function (no constraint point to satisfy)
#1) y=Par[1]*exp(-Par[2]*Tc) with two free parameters
#2) y=Par[1]+Par[2]*exp(-Par[3]*Tc) with three free parameters
exponential_model_2o3 <- function(pars)
{
  if(length(pars)==2)#2 parameters
  {
    est <- pars[1]*exp(-pars[2]*x) #Model estimates
    sumsq<-sum((y - est)^2)#Sum of squared residuals
    return(sumsq)
  }
  if(length(pars)==3)#3 parameters
  {
    est <- pars[1]+pars[2]*exp(-pars[3]*x)#Model estimates
    sumsq<-sum((y - est)^2)#Sum of squared residuals
    return(sumsq)
  }
}

#' @export
peak_ratio_parameters<-function(hourly_streamflows,basin_descriptors,d_threshold=0.2)
{
  #Hourly streamflows with basin descriptors
  hourly_streamflows <- hourly_streamflows[names(hourly_streamflows)%in%basin_descriptors[,1]]
  #Basin descriptors with hourly streanflow series
  basin_descriptors<-basin_descriptors[basin_descriptors[,1]%in%names(hourly_streamflows),]

  hourly_streamflows<-lapply(hourly_streamflows,na.omit)#NA removal

  #Time of concentration estimates
  tc<-time_of_concentration(basin_descriptors)

  #Initialization of the exponential model parameter dataframe for fitting Qhh/Qd=a+b*exp(-c*duration)
  qq_d_pars<- data.frame(matrix(NA,nrow=length(basin_descriptors[,1]),ncol=4,dimnames=list(x=NULL,y=c("name","a","b","c"))))
  qq_d_pars[,1]<-basin_descriptors[,1]

  #Estimation of exponential model parameters
  for (s in names(hourly_streamflows))
  {
    #Hourly and daily streamflows for s-th site
    hourlyflows_s <- hourly_streamflows[[s]]

    #Corresponding daily streamflows
    dailyflows_s <- aggregate(hourlyflows_s,as.Date,mean)

    streamflows<-dailyflows_s
    #s-th site PORFDC
    porfdc_s<-porfdc_eval(dailyflows_s)
    porfdc_s<-unique(porfdc_s) #remove streamflows redundancies
    porfdc_s<-data.frame(p=1:length(porfdc_s)/(length(porfdc_s)+1),
                         Q=porfdc_s,row.names = NULL)

    #Initialization of the dataframe used by the fitting model:
    # x: daily streamflow durations
    # y: maximum hourly streamflow/daily streamflow
    data <- data.frame(matrix(NA,nrow=length(dailyflows_s),ncol=2,dimnames = list(x=NULL,y=c("x","y"))),check.names = FALSE)
    for(Q in unique(dailyflows_s)) data[dailyflows_s%in%Q,1] <- porfdc_s[porfdc_s[,2]%in%Q,1] #daily streamflow durations

    #Maximum hourly streamflow
    max_hourlyflows_s<- aggregate(hourlyflows_s,as.Date,max)

    data[,2]<-zoo::coredata(max_hourlyflows_s/dailyflows_s) #maximum hourly streamflow/daily streamflow

    #Only durations below the threshold, corresponding to ratios between daily and long term streamflows above 1
    data <- data[data[,1]<=d_threshold,]
    data<- data[order(data[,1]),] #Dataframe sorting according to durations

    #Moving Median: 30-day step
    smooth_median <- caTools::runquantile(data[,2],30,probs = 0.5,type=7,endrule="NA")
    smooth_median <- data.frame(x=data[,1],y=smooth_median)
    smooth_median <- smooth_median[-which(is.na(smooth_median[,2])),]

    x<-smooth_median[,1]
    y<-smooth_median[,2]
    #Exponential model fitting: y = a+b*exp(-c*x)
    exp_model <- optim(1,exponential_model_1,method="Brent",upper=2000,lower=-100)

    #parameters writing for s-th site
    qq_d_pars$c[qq_d_pars$Name%in%s] <- exp_model$par
    qq_d_pars$b[qq_d_pars$Name%in%s] <- (smooth_median$y[1] - smooth_median$y[length(smooth_median[,1])])/ (exp(-exp_model$par*smooth_median$x[1])-
                                                                                                              exp(-exp_model$par*smooth_median$x[length(smooth_median$x)]))
    qq_d_pars$a[qq_d_pars$Name%in%s] <- smooth_median$y[1] - qq_d_pars$b[qq_d_pars$Name%in%s]*exp(-smooth_median$x[1]*exp_model$par)
  }

  #Regionalization of exponential models' parameters
  #Initialization of the exponential model parameters dataframe for a=f(Tc), b=f(Tc), c=f(Tc) relationships
  peak_pars <- data.frame(matrix(NA,nrow=3,ncol=5,dimnames = list(x=NULL,y=c("parameter","formula","p1","p2","'p3"))))
  peak_pars$parameter <- c("a","b","c")

  #a=f(Tc)
  y <- qq_d_pars$a
  x <- as.numeric(tc)
  exp_model<-optim(c(1,1),exponential_model_2o3)#2 parameters, sufficient to characterize the trend
  if(length(exp_model$par)==2)
  {
    peak_pars$formula[1]<-"y=p1*exp(-p2*x)"
    peak_pars[1,3:4]<-exp_model$par

  }else
  {
    peak_pars$formula[1]<-"y=p1+p2*exp(-p3*x)"
    peak_pars[1,3:5]<-exp_model$par
  }

  #b=f(Tc)
  y <- qq_d_pars$b
  x <- as.numeric(tc)
  exp_model<-optim(c(1,1,1),exponential_model_2o3)#3 parameters
  if(length(exp_model$par)==2)
  {
    peak_pars$formula[2]<-"y=p1*exp(-p2*x)"
    peak_pars[2,3:4]<-exp_model$par

  }else
  {
    peak_pars$formula[2]<-"y=p1+p2*exp(-p3*x)"
    peak_pars[2,3:5]<-exp_model$par
  }

  #c=f(Tc)
  y <- log(qq_d_pars$c)#logarithmic transformation (more variance of c for the gauged section)
  x <- as.numeric(tc)
  exp_model<-optim(c(1,1),exponential_model_2o3)#2 parameters, sufficient to characterize the trend
  if(length(exp_model$par)==2)
  {
    peak_pars$formula[3]<-"y=p1*exp(-p2*x)"
    peak_pars[3,3:4]<-exp_model$par

  }else
  {
    peak_pars$formula[3]<-"y=p1+p2*exp(-p3*x)"
    peak_pars[3,3:5]<-exp_model$par
  }

  # Save a data frame with functional relationships of a, b, c and the time of concetration
  return(peak_pars)
}
