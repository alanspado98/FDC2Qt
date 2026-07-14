#Estimates of Peak-Duration curve and Reduction Ratio's parameters as a function of time of concentration


# Exponential model function (Tasker and Stedinger, 1989): y = 0.5*exp(-lambda1*x/(1+lambda2*x)) for Flood-volume reduction ratio
exponential_model_2<-function(pars)
{
  lambda1<-pars[1]
  lambda2<-pars[2]
  est<-.5*exp(-lambda1*norm_durations/(1+lambda2*norm_durations))
  sumsq<-sum((est-obs)^2)
  return(sumsq)
}

# Exponential model function (Fiorentino, 1985): y = tau/x*(1-exp(-x/tau)) for Peak-duration curve
exponential_model_1<-function(pars)
{
  tau<-pars[1]
  est<-tau/norm_durations*(1-exp(-norm_durations/tau))
  sumsq<-sum((est-obs)^2)
  return(sumsq)
}

#Linear model function: y = q+m*x with two free parameter
linear_model<-function(y) stats::lm(y ~ x)

#' @export
Df.target<-3
rD_eD_D_pars<-function(hourly_streamflows,basin_descriptors,Df.target,d.threshold)
{
  hourly_streamflows<-lapply(hourly_streamflows,na.omit)#NA removal
  daily_streamflows<-lapply(hourly_streamflows,function(x) aggregate(x,lubridate::ymd,mean,na.rm=T))
  as.Date(zoo::time(hourly_streamflows[[2]]))

  stations<-names(hourly_streamflows)
  durations<-c(1:(Df.target+1))#event durations

  #Lists inizialization
  QD_D<-list()#maximum average discharge in a given duration for each flood event
  rD_D<-list()#time-to-peak for a given duration for each flood event

  `%>%` <- magrittr::`%>%`

  for(s in stations)
  {
    hourly_streamflows_s<-hourly_streamflows[which(hourly_streamflows[,1]%in%r),-1]
    daily_streamflows_s <-

    #PORFDC
    Q=unique(sort(str_GG_r[,2],decreasing = TRUE))
    PORFDC <- cbind.data.frame(p=c((1:length(Q))/(length(Q)+1)),Q=Q)

    str_GG_d <- str_GG_r
    str_GG_d[,3]<-rep(NA,length(str_GG_d[,1]))
    colnames(str_GG_d)[3]<-"d"

    for(q in PORFDC[,2]) str_GG_d[str_GG_d[,2]%in%q,3] <- PORFDC[PORFDC[,2]==q,1] #daily streamflow durations

    #Flood event detection (daily streamflow higher than contiguous days and with duration lower or equal to d_threshold)
    index_flood<-c()
    for(i in 2:(length(str_GG_d[,1])-1))
    {
      if(str_GG_d[i,1]-str_GG_d[i-1,1]==1 && str_GG_d[i+1,1]-str_GG_d[i,1]==1)
      {
        if(str_GG_d[i,2]>=str_GG_d[i-1,2] &&
           str_GG_d[i,2]>=str_GG_d[i+1,2] &&
           str_GG_d[i,3]<=d_threshold) index_flood<-append(index_flood,i)
      }
    }
    index_flood<-unique(index_flood)
    str_GG_d<-str_GG_d[index_flood,]

    QD_D_r <- data.frame(matrix(NA,nrow=length(index_flood),ncol=1+length(durations),dimnames = list(x<-NULL,y<-c("Date",paste0("D_",as.character(durations))))))
    rD_D_r <- data.frame(matrix(NA,nrow=length(index_flood),ncol=1+length(durations),dimnames = list(x<-NULL,y<-c("Date",paste0("D_",as.character(durations))))))
    QD_D_r[,1] <- str_GG_d[,1]
    rD_D_r[,1] <- str_GG_d[,1]
    moving_average <- function(x, n = j){stats::filter(x,rep(1/n,n),sides = 1)}

    for(j in durations)
    {
      vol_r_HH <- cbind.data.frame(Datetime=str_HH_r[,1],Value=moving_average(str_HH_r[,2]))
      vol_r_HH[1:(length(vol_r_HH[,1])-j+1),2] <- vol_r_HH[j:length(vol_r_HH[,1]),2]
      vol_r_HH<-vol_r_HH[1:(length(vol_r_HH[,1])-j+1),]

      for(e in 1:length(str_GG_d[,1]))
      {
        index_vol <- which(substring(vol_r_HH[,1],1,10)%in%c(str_GG_d[e,1]-1,str_GG_d[e,1],str_GG_d[e,1]+1))
        index_str <- which(substring(str_HH_r[,1],1,10)%in%c(str_GG_d[e,1]-1,str_GG_d[e,1],str_GG_d[e,1]+1))

        t_start <- vol_r_HH[index_vol,1][which.max(vol_r_HH[index_vol,2])]
        t_peak  <- str_HH_r[index_str,1][which.max(str_HH_r[index_str,2])]

        #Saving the results
        rD_D_r[rD_D_r[,1]%in%str_GG_d[e,1],colnames(rD_D_r)%in%paste("D",as.character(j),sep="_")] <- as.numeric(t_peak-t_start)/j
        QD_D_r[QD_D_r[,1]%in%str_GG_d[e,1],colnames(QD_D_r)%in%paste("D",as.character(j),sep="_")] <- max(vol_r_HH[index_vol,2],na.rm=T)
      }
    }
    QD_D<-append(QD_D,list(QD_D_r))
    rD_D<-append(rD_D,list(rD_D_r))
  }
  names(QD_D) <- stations#station names
  names(rD_D) <- stations#station names

  #Shift -1 hour (in accordance to Maione et al., 2003)
  durations<-c(0.001,durations[-length(durations)])

  ## Reduction ratio (time to flood peak)
  rD_D_mean<-data.frame(matrix(NA,1,length(durations)))
  for(r in stations)
  {
    rD_D_r<-rD_D[[r]][,-1]
    rD_D_r_mean<-c()
    for(c in 1:length(durations))
    {
      rD_D_r_c<-rD_D_r[,c]
      rD_D_r_c<-rD_D_r_c[which(rD_D_r_c>0 & rD_D_r_c<1)]
      rD_D_r_mean<-append(rD_D_r_mean,mean(rD_D_r_c))
    }
    rD_D_mean<-rbind.data.frame(rD_D_mean,rD_D_r_mean)
  }
  rD_D_mean<-rD_D_mean[-1,]#First row: NAs
  rD_D_mean <- rD_D_mean[,-1]#First column: NAs
  rownames(rD_D_mean)<-stations
  colnames(rD_D_mean)<-paste("D",as.character(durations[-length(durations)]),sep="_")

  basin_descriptors <- basin_descriptors[basin_descriptors[,1]%in%rownames(rD_D_mean),]
  rD_D_mean<-rD_D_mean[order(factor(rownames(rD_D_mean),levels=basin_descriptors[,1])),]

  #Time of concentration estimates
  descriptors_tc<-basin_descriptors[,1:5]
  ToC<-Time_of_Concentration(descriptors_tc)

  #Initialization lambda1 and lambda2 dataframe
  rD_D_pars<-data.frame(matrix(NA,length(rD_D_mean[,1]),3,dimnames=list(x=NULL,y=c("Name","lambda1","lambda2"))))
  rD_D_pars[,1]<-rownames(rD_D_mean)

  for(r in stations)
  {
    norm_durations<-durations[-length(durations)]/ToC[names(ToC)%in%r]
    obs<-as.numeric(rD_D_mean[rownames(rD_D_mean)%in%r,])
    rD_D_model<-optim(c(1,1),exponential_model_2)
    rD_D_pars[rD_D_pars[,1]%in%r,-1]<-rD_D_model$par
  }

  #Initialization of the exponential model parameters of the functional relationships: lambda1=f(Tc), lambda2=f(Tc)
  λs_Tc_pars <- data.frame(matrix(NA,nrow=2,ncol=4,dimnames = list(x=NULL,y=c("Parameter","Formula","q","m"))))
  λs_Tc_pars[,1] <- colnames(rD_D_pars)[-1]
  x<-as.numeric(ToC)

  for(i in 1:2) λs_Tc_pars[i,3:4]<-linear_model(rD_D_pars[,i+1])$coefficients
  λs_Tc_pars[,2]<-"y=q+m*x"

  ## Peak-duration curve (maximum flood volume in D/maximum flood volume in 1h)
  εD_D_mean<-data.frame(matrix(NA,1,length(durations)))

  for(r in stations)
  {
    QD_D_r<-QD_D[[r]]
    εD_D_r<-QD_D_r[,-1]/QD_D_r[,2]

    εD_D_r_mean<-c()
    for(c in 1:length(durations))
    {
      εD_D_r_c<-εD_D_r[,c]
      εD_D_r_c<-εD_D_r_c[which(εD_D_r_c>=0 & εD_D_r_c<=1)]
      εD_D_r_mean<-append(εD_D_r_mean,mean(εD_D_r_c))
    }
    εD_D_mean<-rbind.data.frame(εD_D_mean,εD_D_r_mean)
  }
  εD_D_mean<-εD_D_mean[-1,]#First row: NAs
  rownames(εD_D_mean)<-stations
  colnames(εD_D_mean)<-paste("D",as.character(durations),sep="_")

  # Initialization tau dataframe
  εD_D_par<-data.frame(matrix(NA,12,2,dimnames=list(x=NULL,y=c("Name","tau"))))
  εD_D_par[,1]<-rownames(εD_D_mean)

  for(r in stations)
  {
    norm_durations<-durations/ToC[names(ToC)%in%r]
    obs<-as.numeric(εD_D_mean[rownames(εD_D_mean)%in%r,])
    εD_D_model<-optim(1,exponential_model_1,method = "Brent",lower = 0.001,upper = 100)
    εD_D_par[εD_D_par[,1]%in%r,-1]<-εD_D_model$par
  }

  #Initialization of the exponential model parameters of the functional relationship: tau=f(Tc)
  tau_Tc_pars <- data.frame(matrix(NA,nrow=1,ncol=4,dimnames = list(x=NULL,y=c("Parameter","Formula","q","m"))))
  tau_Tc_pars[,1] <- colnames(εD_D_par)[-1]
  x<-as.numeric(ToC)

  tau_Tc_pars[,3:4]<-linear_model(εD_D_par[,2])$coefficients
  tau_Tc_pars[,2]<-"y=q+m*x"

  #Store functional relationships of reduction ratio and peak-duration curve's parameters inside a list
  rD_εD_D_pars <- list(λs_Tc_pars=λs_Tc_pars,tau_Tc_pars=tau_Tc_pars)
  return(rD_εD_D_pars)
}

