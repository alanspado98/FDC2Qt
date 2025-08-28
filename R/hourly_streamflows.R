#' Calculate hourly streamflows of the target station
#'
#' This function takes .
#'
#' @param target_section a character with the cross-section name where to estimate streamflows
#' @return  
#' @export

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

#Hourly hydrograph in case duration < d_threshold and daily streamflow above previous and next day
# Two solution:
# - donator hourly streamflow ->  
#y = a*x^3 + b*x^2 + c*x + d

polynomial_1 <- function(c)
{
  #Matrice dei coefficienti (imponendo i punti 1),2) e 3))
  # A <- matrix(c(x[k]^3-3*x[k]*x[which.max(y)]^2, x[k]^2-2*x[k]*x[which.max(y)],
  #               -2*x[which.max(y)]^3,-x[which.max(y)]^2),
  #             2, 2, byrow=TRUE)
  A <- matrix(c(x[k]^2,x[k],
                x[which.max(y)]^2,x[which.max(y)]),
              2, 2, byrow=TRUE)
  
  colnames(A) <- c("a","b")#Parametri matrice coefficienti (a,b)
  
  #Vettore dei termini noti
  b <- c(y[k]-c,y[which.max(y)]-c)#Parametro vettore termini noti (d) 
  
  #showEqn(A, b) #per mostrare l'equazione
  #Matrice dei coefficienti e dei termini noti devono avere stesso rango
  #c(R(A),R(cbind(A,b))) #rango
  #all.equal( R(A), R(cbind(A,b))) TRUE -> E' consistente
  
  Par <- solve(A, b) #Restituisce a, b
  
  #c<-as.numeric(-3*Par[1]*x[which.max(y)]^2-2*Par[2]*x[which.max(y)]) #c dall'equazione che impone max nell'ora del colmo
  
  P <- c(Par,c) #a,b,c,d
  
  if (k==1 || k==2)#Ramo ascendente per CASI ORDINARI e PARTICOLARI (k=1) o discendente per CASI PARTICOLARI (k=2) 
  {
    t_hh <- seq(from=1,to=pos_max,by=1)#Ore dell'evento di piena
    Q_24h <- P[1]*t_hh^2 + P[2]*t_hh + P[3]#P[1]*t_hh^3 + P[2]*t_hh^2 + P[3]*t_hh + P[4]#Portate orarie dell'evento di piena
  }
  
  if (k==3)#Ramo discendente per CASI ORDINARI e PARTICOLARI
  {
    t_hh <- seq(from=pos_max,to=49,by=1)#Ore dell'evento di piena
    Q_24h <- P[1]*t_hh^2 + P[2]*t_hh + P[3]#P[1]*t_hh^3 + P[2]*t_hh^2 + P[3]*t_hh + P[4]#Portate orarie dell'evento di piena
  }
  
  if(result=="Calibration")#Restituizione differenza volumetrica
  {
    Q_sdh<-SDH_Maione(k)
    Volume_hh <- trapz(t_hh,Q_24h)#Integrazione dell'idrogramma orario
    
    Delta_Q<-sum(abs((P[1]*Q_sdh$x^2 + P[2]*Q_sdh$x + P[3])-Q_sdh$y)/Q_sdh$y)/length(Q_sdh$y)
    Delta_V <- abs(Volume_hh/Volume_gg-1)#Differenza tra i volumi dell'idrogramma giornaliero e di quello orario
    
    Delta<-Delta_Q*λ_Q+Delta_V*λ_V
    return(Delta)
  }
  
  if(result=="Series")#Restituizione della serie oraria d'evento
  {
    return(Q_24h)  
  }
}  

hourly_streamflows<-function(target_section)
{
  #Observed hourly and daily streamflow series
  hourly_str<-use_internal_data(list_name="streamflows",dataset_name ="HH")
  daily_str<-use_internal_data(list_name="streamflows",dataset_name ="GG")
  
  #Regional daily streamflow series
  reg_daily_str<-use_internal_data(list_name="reg_streamflows",dataset_name ="GG")
}


