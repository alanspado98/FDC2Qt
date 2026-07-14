#' Synthetic hourly streamflow series of the target site
#' @param reg_daily_str
#'
#' @param target_section
#' @param abc_Tc_pars
#' @param τ_Tc_pars
#' @param λs_Tc_pars
#' @param Df_target
#' @param d_threshold
#'
#' @export
#'
#'
#' @param target_section a character with the cross-section name where to estimate streamflows
#' @return

#Hourly hydrograph in case of flood events: daily streamflows with duration < d_threshold and daily streamflow value higher than previous and next day
#1) Synthetic Design Hydrograph (Maione et al, Int. J. River Basin Management, 2003)
synthetic_design_hydrograph<-function(k)
{
  if(k==1)
  {
    t<--rD*durations
    Q_hh<-(drD_dD*durations*QD+rD*QD+rD*durations*dQD_dD)/(drD_dD*durations+rD)
  }

  if(k==2 || k==3)
  {
    t<-durations*(1-rD)
    Q_hh<-(-drD_dD*durations*QD+(1-rD)*QD+(1-rD)*durations*dQD_dD)/(-drD_dD*durations+1-rD)
  }
  t<-t+pos_max
  sdh<-data.frame(x=t,y=Q_hh)

  return(sdh)
}

#2) constrained polynomial y = a*x^2 + b*x + c for representing rising and falling limbs, by imposing:
#   - first and last value equal to daily streamflow of previous and next days
#   - passage through the peak (given by h2d_d_pars.R)
#   - minimize the difference with respect to Synthetic Design Hydrograph
#   - minimize the difference with respect to daily flood volume
constrained_polynomial <- function(c)
{
  A <- matrix(c(x[k]^2,x[k],#Coefficient matrix
                x[which.max(y)]^2,x[which.max(y)]),
              2, 2, byrow=TRUE)

  colnames(A) <- c("a","b")#parameters

  b <- c(y[k]-c,y[which.max(y)]-c)#Vector of known terms

  Par <- solve(A, b)

  P <- c(Par,c) #parameters vector

  if (k==1)#rising limb
  {
    t_hh <- seq(from=1,to=pos_max,by=1)#time vector of the event
    Q_24h <- P[1]*t_hh^2 + P[2]*t_hh + P[3]#hourly streamflows of the flood event
  }

  if (k==3)#falling limb
  {
    t_hh <- seq(from=pos_max,to=49,by=1)#time vector of the event
    Q_24h <- P[1]*t_hh^2 + P[2]*t_hh + P[3]#hourly streamflows of the flood event
  }

  if(result=="Calibration")#parameters calibration: returns a weighted average between differences with daily flood volumes and SDH
  {
    Q_sdh<-synthetic_design_hydrograph(k)
    Volume_hh <- pracma::trapz(t_hh,Q_24h)

    Delta_Q<-sum(abs((P[1]*Q_sdh$x^2 + P[2]*Q_sdh$x + P[3])-Q_sdh$y)/Q_sdh$y)/length(Q_sdh$y)
    Delta_V <- abs(Volume_hh/Vol_GG-1)

    Delta<-(Delta_Q+Delta_V)*0.5
    return(Delta)
  }

  if(result=="Generation")#generation of hourly streamflow series of the event with calib parameters
  {
    return(Q_24h)
  }
}

synthetic_hourly_streamflow<-function(,target,abc_Tc_pars,tau_Tc_pars,lambda_Tc_pars,
                                      Df.target,d.threshold,mae.threshold=0.25,max.iter=20)
{
  (s)

  # Regional period-of-record flow duration curve
  reg_PORFDC<-unique(reg_daily_str[,2]) #remove streamflows redundancies
  reg_PORFDC<-data.frame(p=1:length(reg_PORFDC)/(length(reg_PORFDC)+1),Q=sort(reg_PORFDC,decreasing = TRUE))

  # Grouping consecutive years of the regional daily streamflow series
  Group_cons_years<-split(unique(lubridate::year(reg_daily_str[,1])),cumsum(c(1,diff(as.numeric(unique(lubridate::year(reg_daily_str[,1]))))!= 1)))

  # Time of concentration for target site
  descriptors_target_tc<-basin_descriptors[names(ToC)%in%target_section,1:5]
  ToC_target<-Time_of_Concentration(descriptors_target_tc)

  # Hourly to daily peak ratio parameters estimation
  abc_target <- c()#Parameters' vector initialization
  x<-as.numeric(ToC_target) #where to evaluate the formulas
  for(i in 1:length(abc_Tc_pars[,1]))
  {
    if(is.na(abc_Tc_pars[i,5])) {P1<-abc_Tc_pars[i,3];P2<-abc_Tc_pars[i,4]} else {P1<-abc_Tc_pars[i,3];P2<-abc_Tc_pars[i,4];P3<-abc_Tc_pars[i,5]}
    abc_target <- append(abc_target,eval(parse(text=substring(abc_Tc_pars[i,2],first=3,last=nchar(abc_Tc_pars[i,2])))))
  }

  # Reduction ratio parameters
  λs_target<-c()
  for(i in 1:length(λs_Tc_pars[,1]))
  {
    q<-λs_Tc_pars[i,3];m<-λs_Tc_pars[i,4]
    λs_target <- append(λs_target,eval(parse(text=substring(λs_Tc_pars[i,2],first=3,last=nchar(λs_Tc_pars[i,2])))))
  }

  # Peak-duration curve parameter
  τ_target <- τ_Tc_pars[,3]+τ_Tc_pars[,4]*x

  durations<-c(0.001,1:Df_target)#durations
  norm_durations<-durations/ToC_target#normalized durations (by Time of Concentration)

  #Synthetic design hydrograph's contructing variables
  rD<-.5*exp(-λs_target[1]*norm_durations/(1+λs_target[2]*norm_durations))#Reduction ratio
  εD<-τ_target/norm_durations*(1-exp(-norm_durations/τ_target))#Peak-duration curve
  drD_dD<--((λs_target[1]/ToC_target)/(1+λs_target[2]*norm_durations)^2)*rD#first derivative of the reduction ratio (with respect to event duration)
  dεD_dD<-1/norm_durations*(-τ_target/(ToC_target*norm_durations)*(1-exp(-norm_durations/τ_target))+(1/ToC_target*exp(-norm_durations/τ_target)))#first derivative of the peak-duration curve (with respect to event duration)

  for(l in 1:length(Group_cons_years))#cycle on each group of consecutive years
  {
    str_GG_l<-reg_daily_str[lubridate::year(reg_daily_str[,1])%in%Group_cons_years[[l]],] #daily streamflow in group l
    rownames(str_GG_l)<-NULL
    str_GG_l <- cbind.data.frame(str_GG_l,p=rep(NA,length(str_GG_l[,1])))
    for(q in reg_PORFDC[,2]) str_GG_l[str_GG_l[,2]%in%q,3] <- reg_PORFDC[reg_PORFDC[,2]==q,1] #daily streamflow durations

    str_HH_l<-rep(NA,length(str_GG_l[,1])*24-12-11)#Initilization of the hourly series

    #Linear interpolation of daily streamflow series
    for(i in 1:(length(str_GG_l[,1])-1)) str_HH_l[((i-1)*24+1):(i*24+1)] <- approx(x=i:(i+1),y=str_GG_l[c(i,i+1),2],xout=seq(i,i+1,1/24),method="linear")$y

    #temporal axis to hourly streamflow series
    str_HH_l <- data.frame(Datetime=rep(NA,length(str_HH_l)),Value=str_HH_l)
    str_HH_l[,1]<-seq.POSIXt(from=lubridate::ymd_hms(paste0(str_GG_l[1,1]," 12:00:00")),
                             to=lubridate::ymd_hms(paste0(str_GG_l[length(str_GG_l[,1]),1]," 12:00:00")),
                             by="1 hour",tz="GMT")


    #Flood events (d<=d_threshold and daily streamflow higher than contiguous days): hourly flood peak fixed at 12:00
    #and 2 polynomials are calibrated: y=a*t^2+b*t+c for rapresenting the rising and falling limbs; the three parameters
    #are estimated by constraining the passage through the hourly flood peak, the passage through the daily streamflows of the previous and
    # the following day, the respect of the daily flood volumes and of the shape of the hydrograph

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


    for(i in index_flood)#cycle on day
    {
      Qmax<-(abc_target[1]+abc_target[2]*exp(-abc_target[3]*str_GG_l[i,3]))*str_GG_l[i,2]#hourly flood peak
      QD<-εD*Qmax #maximum flood volumes at varying flood durations
      dQD_dD<-dεD_dD*Qmax #first derivative of QD

      pos_max<-25#Position of the flood peak
      x<-c(1,pos_max,49)#hours in which hourly streamflow is known
      y<-c(str_GG_l[i-1,2],Qmax,str_GG_l[i+1,2])#known hourly streamflows

      #Polynomial paramters' estimates, for the rising (1st>25th hour) and for the falling (25th>49th hour) of the flood event
      result<-"Calibration"
      c_opt<-c()#Initializing the variable of the optima c values of the polynomials
      lower_bound <- -10#lower end point of the interval where to search
      upper_bound <- 10 #upper end point of the interval where to search

      for(k in c(1,3))#Indice di x e y relativo alle portate orarie note oltre al colmo
      {
        if(k==1) Vol_GG<-12*str_GG_l[i-1,2]+(pos_max-13)*str_GG_l[i,2]#Volumi idrogramma giornaliero
        if(k==3) Vol_GG<-12*str_GG_l[i+1,2]+(37-pos_max)*str_GG_l[i,2]

        dV_above_threshold <- TRUE
        cont <- 0
        while(dV_above_threshold && cont < max_iter)
        {
          cont<-cont+1
          Delta <- optimize(constrained_polynomial,lower=lower_bound,upper=upper_bound)

          if (Delta$objective < Delta_threshold) {
            dV_above_threshold <- FALSE
          } else
          {
            par <- Delta$minimum
            # Bounds expanded adaptively
            if (abs(par - upper_bound) <= abs(par - lower_bound)) {
              lower_bound <- par
              upper_bound <- upper_bound + abs(upper_bound) * 0.5
            } else {
              upper_bound <- par
              lower_bound <- lower_bound - abs(lower_bound) * 0.5
            }
          }
        }
        c_opt <- append(c_opt,par)
      }

      #Hourly streamflow series generation
      result<-"Generation"
      k<-1; Q1<-constrained_polynomial(c_opt[1])#Rising limb
      k<-3; Q2<-constrained_polynomial(c_opt[2])#Falling limb

      str_HH_flood <- c(Q1,Q2[-1])#Limbs merging (First and last streamflow of the two vectors are identical)

      #Negative streamflows: linearization
      #Rising limb
      slope     <- (str_GG_l[i-1,2] - str_GG_l[i-2,2])/24
      intercept <- str_GG_l[i-1,2] - slope*1
      str_HH_lin<- slope*c(1:pos_max)+intercept
      if(length(which(str_HH_flood[1:pos_max]<str_HH_lin))>0) str_HH_flood[which(str_HH_flood[1:pos_max]<str_HH_lin)]<-str_HH_lin[which(str_HH_lin>str_HH_flood[1:pos_max])]
      #Falling limb
      slope     <- (str_GG_l[i+2,2]-str_GG_l[i+1,2])/24
      intercept <-  str_GG_l[i+1,2] - slope*49
      str_HH_lin<- slope*c(pos_max:49)+intercept
      if(length(which(str_HH_flood[pos_max:49]<str_HH_lin))>0) str_HH_flood[which(str_HH_flood[pos_max:49]<str_HH_lin)+pos_max-1]<-str_HH_lin[which(str_HH_lin>str_HH_flood[pos_max:49])]

      #Update the hourly streamflow series obtained by linear interpolation with the flood hydrograph
      str_HH_l[which(lubridate::as_date(str_HH_l[,1])%in%str_GG_l[i-1,1]&lubridate::hour(str_HH_l[,1])%in%12):
                 which(lubridate::as_date(str_HH_l[,1])%in%str_GG_l[i+1,1]&lubridate::hour(str_HH_l[,1])%in%12),2]<-str_HH_flood
    }

    if(l==1)#First group of consecutive years
    {
      reg_hourly_str<-str_HH_l
    }else #Following groups
    {
      reg_hourly_str<-rbind.data.frame(reg_hourly_str,str_HH_l)
    }
  }

  #regional hourly streamflows dataframe
  return(reg_hourly_str)
}



