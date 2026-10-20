#' Peak-Duration curve and Reduction Ratio's parameters
#'
#' Computes peak-duration and reduction ratio parameters as a function
#' of associated event durations, fit an exponential model and regionalizes its
#'  parameters according to time of concentration
#'
#' @param hourly_streamflows A named list of zoo time series (one per site). Each
#'    element must be a zoo object whose time index is of class `POSIXct`.
#'    The list names represent site identifiers.
#' @param basin_descriptors A data frame whose rows represents sites, columns
#'   multiple columns: station name in the first (`character`) and morpho‑climatic
#'   descriptors in the remaining columns (`numeric`).
#' @param D_target A `numeric` of high-flow events duration at the target section
#'   in hours (by default is set to 25).
#' @param d_threshold A `numeric` of duration threshold for high-flow events,
#'   evaluated in [(regional_porfdc)] (by default is set to 0.2).
#'
#'#' @return A data frame of with the exponential-model parameters as a function of time of concentration

# 2-parameter exponential model (y = 0.5*exp(-lambda1*x/(1+lambda2*x)), for peak-duration curve fitting)
exponential_model_2<-function(pars)
{
  lambda1<-pars[1]
  lambda2<-pars[2]
  est<-.5*exp(-lambda1*norm_Durations/(1+lambda2*norm_Durations))
  sumsq<-sum((est-obs)^2)
  return(sumsq)
}

# 1-parameter exponential model (y = tau/x*(1-exp(-x/tau)), for reduction ratio fitting)
exponential_model_1<-function(pars)
{
  tau<-pars[1]
  est<-tau/norm_Durations*(1-exp(-norm_Durations/tau))
  sumsq<-sum((est-obs)^2)
  return(sumsq)
}

# 2-Parameter linear model (y = q+m*x, for reduction ratio and peak-duration curve parameters fitting)
linear_model<-function(y) stats::lm(y ~ x)

#' @export
shape_ratio_parameters<-function(hourly_streamflows,basin_descriptors,D_target=25,d_threshold=0.2)
{
  #Hourly streamflows with basin descriptors
  hourly_streamflows <- hourly_streamflows[names(hourly_streamflows)%in%basin_descriptors[,1]]
  #Basin descriptors with hourly streanflow series
  basin_descriptors<-basin_descriptors[basin_descriptors[,1]%in%names(hourly_streamflows),]

  #Time of concentration (used later for regionalization)
  tc<-time_of_concentration(basin_descriptors)

  hourly_streamflows<-lapply(hourly_streamflows,na.omit)#NA removal
  daily_streamflows<-lapply(hourly_streamflows,function(x) aggregate(x,as.Date,mean,na.rm=T))

  stations<-names(hourly_streamflows)
  Durations<-c(1:(D_target+2))#event durations (to shift durations vector and keep D_target)

  #Lists inizialization
  QD_D<-list()#maximum average discharge in given durations for each high-flow condition
  rD_D<-list()#time-to-peak associated to the maximum average discharge in given durations for each high-flow condition

  for(s in stations)
  {
    hourlyflows_s<-hourly_flows[[s]]

    porfdc_s <- porfdc_eval(daily_streamflows[[s]])#porfdc
    flows_d <- data.frame(date=zoo::index(daily_streamflows[[s]]),
                          value=zoo::coredata(daily_streamflows[[s]]),
                          p=NA)

    #Daily streamflow durations
    for(q in unique(reg_porfdc)) flows_d[flows_d[,2]%in%q,3] <- mean(zoo::index(reg_porfdc[reg_porfdc%in%q]))

    #High-flows detection: daily streamflows higher than contiguous days and
    #with duration lower or equal to d_threshold)
    hf_index<-c()
    for(i in 2:(length(flows_d[,1])-1))
    {
      if(flows_d[i,1]-flows_d[i-1,1]==1 && flows_d[i+1,1]-flows_d[i,1]==1)#check on consecutive days
      {
        if(flows_d[i,2]>=flows_d[i-1,2] &&
           flows_d[i,2]>=flows_d[i+1,2] &&
           flows_d[i,3]<=d_threshold) hf_index<-append(hf_index,i)
      }
    }
    highflows_d <- flows_d[hf_index,]

    #Dataframe initialization
    QD_D_s <- data.frame(matrix(NA,nrow=length(hf_index),ncol=1+length(Durations),#maximum average discharge at site s
                                dimnames = list(x<-NULL,y<-c("date",paste0("D_",as.character(Durations))))))
    rD_D_s <- data.frame(matrix(NA,nrow=length(hf_index),ncol=1+length(Durations),#time-to-peak at site s
                                dimnames = list(x<-NULL,y<-c("date",paste0("D_",as.character(Durations))))))
    QD_D_s[,1] <- highflows_d[,1]
    rD_D_s[,1] <- highflows_d[,1]

    for(D in Durations)
    {
      hourly_means_s<-rollmean(hourlyflows_s,D,align="left")

      for(i in 1:length(hf_index))
      {
        before  <- as.POSIXct(format(highflows_d[i,1] - 1))#day before
        after  <- as.POSIXct(format(highflows_d[i,1] + 2))-1#day after

        means_s_date<-window(hourly_means_s, start = before, end = after)#volumes subset
        flows_s_date<-window(hourlyflows_s, start = before, end = after)#streamflows subset

        t_start <- zoo::index(means_s_date[which.max(means_s_date)])#time of maximum volume
        t_peak  <- zoo::index(flows_s_date[which.max(flows_s_date)])#time of maximum streamflow

        #Saving the results
        rD_D_s[rD_D_s[,1]%in%highflows_d[i,1],colnames(rD_D_s)%in%paste("D",as.character(D),sep="_")] <- (as.numeric(t_peak-t_start)+1)/D
        QD_D_s[QD_D_s[,1]%in%highflows_d[i,1],colnames(QD_D_s)%in%paste("D",as.character(D),sep="_")] <- means_s_date[which.max(means_s_date)]
      }
    }
    QD_D<-append(QD_D,list(QD_D_s))
    rD_D<-append(rD_D,list(rD_D_s))
  }
  names(QD_D)<-stations; names(rD_D)<-stations#station names

  #Shift -2 hours (in accordance to Maione et al., 2003)
  Durations<-c(0.001,Durations[-c(length(Durations)-1,length(Durations))])

  # Peak-duration curve
  # At-site computation
  rD_D_mean<-data.frame(matrix(NA,0,length(rD_D[[s]][1,-1]),dimnames=list(x=NULL,y=colnames(rD_D[[s]])[-1])))
  for(s in stations)
  {
    rD_D_s<-rD_D[[s]][,-1]
    rD_D_s_mean<-colMeans(replace(rD_D_s, rD_D_s<=0 | rD_D_s>=1, NA), na.rm = TRUE)
    rD_D_mean<-rbind.data.frame(rD_D_mean,data.frame(t(rD_D_s_mean)))
  }
  rownames(rD_D_mean)<-stations
  rD_D_mean<-rD_D_mean[,-1]#all NAs for 1-h Duration since t_peak-t_start=1
  colnames(rD_D_mean)<-paste("D",as.character(Durations),sep="_")

  rD_D_mean<-rD_D_mean[order(factor(rownames(rD_D_mean),levels=descriptors[,1])),]

  #Peak-duration curve parameters
  rD_D_pars<-data.frame(matrix(NA,length(rD_D_mean[,1]),3,dimnames=list(x=NULL,y=c("name","lambda1","lambda2"))))
  rD_D_pars[,1]<-rownames(rD_D_mean)
  for(s in stations)
  {
    norm_Durations<-Durations/tc[names(tc)%in%s]
    obs<-as.numeric(rD_D_mean[rownames(rD_D_mean)%in%s,])
    rD_D_model<-optim(c(1,1),exponential_model_2)
    rD_D_pars[rD_D_pars[,1]%in%s,-1]<-rD_D_model$par
  }

  # Reduction ratio
  # At-site computation
  epsD_D_mean<-data.frame(matrix(NA,0,length(QD_D[[s]][1,-1]),dimnames=list(x=NULL,y=colnames(QD_D[[s]])[-1])))
  for(s in stations)
  {
    QD_D_s<-QD_D[[s]]
    epsD_D_s<-QD_D_s[,-1]/QD_D_s[,2]
    epsD_D_s_mean<-colMeans(replace(epsD_D_s,epsD_D_s<=0 | epsD_D_s>=1, NA), na.rm = TRUE)
    epsD_D_mean<-rbind.data.frame(epsD_D_mean,data.frame(t(epsD_D_s_mean)))
  }
  rownames(epsD_D_mean)<-stations
  epsD_D_mean<-epsD_D_mean[,-1]#all NAs for 1h-Duration since t_peak-t_start=1
  colnames(epsD_D_mean)<-paste("D",as.character(Durations),sep="_")

  #Reduction ratio parameters
  epsD_D_par<-data.frame(matrix(NA,length(epsD_D_mean[,1]),2,dimnames=list(x=NULL,y=c("name","tau"))))
  epsD_D_par[,1]<-rownames(epsD_D_mean)

  for(s in stations)
  {
    norm_Durations<-Durations/tc[names(tc)%in%s]
    obs<-as.numeric(epsD_D_mean[rownames(epsD_D_mean)%in%s,])
    epsD_D_model<-optim(1,exponential_model_1,method = "Brent",lower = 0.001,upper = 100)
    epsD_D_par[epsD_D_par[,1]%in%s,2]<-epsD_D_model$par
  }

  #Regionalization of shape parametrs
  #Peak-duration curve
  shape_pars <- data.frame(matrix(NA,nrow=3,ncol=4,dimnames = list(x=NULL,y=c("parameter","formula","q","m"))))
  shape_pars[-3,1] <- colnames(rD_D_pars)[-1]
  x<-as.numeric(tc)
  shape_pars[,2]<-"y=q+m*x"
  for(i in 1:2) shape_pars[i,3:4]<-linear_model(rD_D_pars[,i+1])$coefficients
  #Reduction ratio
  shape_pars[3,1] <- colnames(epsD_D_par)[2]
  shape_pars[3,3:4]<-linear_model(epsD_D_par[,2])$coefficients

  # Store a data frame with functional relationships of lambda1, lambda2, tau and the time of concetration
  return(shape_pars)
}
