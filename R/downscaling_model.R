#'  Hourly synthetic streamflow series for the target site
#'
#' This function takes .
#'
#' @param target_section a character with the cross-section name where to estimate streamflows
#' @return  

#Hourly hydrograph in case of flood events: daily streamflows with duration < d_threshold and daily streamflow value higher than previous and next day
#1) Synthetic Design Hydrograph (Maione et al, Int. J. River Basin Management, 2003)
SDH_Maione<-function(k)
{
  if(k==1) 
  {
    t<--rd*durations
    Q_hh<-(drd_dd*durations*Qd+rd*Qd+rd*durations*dQd_dd)/(drd_dd*durations+rd)
  }  
  
  if(k==2 || k==3) 
  {
    t<-durations*(1-rd)
    Q_hh<-(-drd_dd*durations*Qd+(1-rd)*Qd+(1-rd)*durations*dQd_dd)/(-drd_dd*durations+1-rd)
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
  
  if (k==1 || k==2)#rising limb for "ordinary" and "particular" cases (k=1, see later) o falling limb for "particular" cases (k=2, see later) 
  {
    t_hh <- seq(from=1,to=pos_max,by=1)#time vector of the event
    Q_24h <- P[1]*t_hh^2 + P[2]*t_hh + P[3]#hourly streamflows of the flood event
  }
  
  if (k==3)#falling limb for "ordinary" and "particular" cases
  {
    t_hh <- seq(from=pos_max,to=49,by=1)#time vector of the event
    Q_24h <- P[1]*t_hh^2 + P[2]*t_hh + P[3]#hourly streamflows of the flood event
  }
  
  if(result=="Calibration")#parameters calibration: returns a weighted average between differences with daily flood volumes and SDH
  {
    Q_sdh<-SDH_Maione(k)
    Volume_hh <- trapz(t_hh,Q_24h)
    
    Delta_Q<-sum(abs((P[1]*Q_sdh$x^2 + P[2]*Q_sdh$x + P[3])-Q_sdh$y)/Q_sdh$y)/length(Q_sdh$y)
    Delta_V <- abs(Volume_hh/Volume_gg-1)
    
    Delta<-Delta_Q*λ_Q+Delta_V*λ_V
    return(Delta)
  }
  
  if(result=="Prediction")#generation of hourly streamflow series of the event with calib parameters
  {
    return(Q_24h)  
  }
}  

#' @export
Df_target<-30
downscaling_model<-function(reg_daily_str,target_section,abc_Tc_pars,ToC,tau_Tc_pars,lambda_Tc_pars,Df_target)
{
  ToC_target<-as.numeric(ToC[names(ToC)%in%target_section])
  durations<-c(.001,1:31)
}


