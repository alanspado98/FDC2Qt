#Hourly synthetic streamflow series
load("data/basin_descriptors.rda")
#Time of concentration, Tc (Giandotti, 1934)
ToC_Giandotti<-function(descr)
{
  tc<-round((4*sqrt(descr[,2])+1.5*descr[,3])/(0.8*sqrt(descr[,4]-descr[,5])),digits=2)
  names(tc)<-descr[,1]
  return(tc)
}

#' @export
pars_downscaling<-function(basin_descriptors,target_section)
{
  #Observed hourly and daily streamflow series
  streamflows_HH<-use_internal_data(list_name="streamflows",dataset_name ="HH")
  streamflows_HH<-streamflows_HH[-which(is.na(streamflows_HH[,3])),]
  streamflows_GG<-use_internal_data(list_name="streamflows",dataset_name ="GG")
  colnames(streamflows_HH)[1:2]<-c("Name","Datetime")
  
  #Regional daily streamflow series
  reg_streamflows_GG<-use_internal_data(list_name="reg_streamflows",dataset_name ="GG")

  basin_descriptors<-basin_descriptors[basin_descriptors[,1]%in%unique(streamflows_HH[,1]),]
  
  #Time of concentration: basin descriptors for hourly streamflow series
  Index_ToC <- which(sub("^[^_]*_","",colnames(basin_descriptors))%in%"h")
  if(length(Index_ToC)==0)
  {
    descriptors_tc<-basin_descriptors[,1:5]
    ToC<-ToC_Giandotti(descriptors_tc)
  }else 
  {
    ToC<-basin_descriptors[,Index_ToC]
    names(ToC)<-basin_descriptors[,1]
  }  
    
  #Initialization of the exponential model parameter dataframe for fitting Qhh/Qd=a+b*exp(-c*duration)
  qq_d_pars<- data.frame(matrix(NA,nrow=length(basin_descriptors[,1]),ncol=4,dimnames=list(x=NULL,y=c("Name","a","b","c")))) 
  qq_d_pars[,1]<-descriptors_tc[,1]
  
  #Exponential model function: y=a+b*exp(-c*x) with one free parameter (first and last point must be equivalent with smooth_median ones)
  exponential_model_1 <- function(c)
  {
    b <- (y[1]-y[length(x)])/(exp(-c*x[1])-exp(-c*x[length(x)]))#a,b model parameters
    a <- y[1] - b*exp(-c*x[1])
    est <- a + b*exp(-c*x)#Model estimates 
    sumsq<-sum((y - est)^2) #Sum of squared residuals
    return(sumsq)
  }

  #Estimation of exponential model parameters 
  for (r in unique(streamflows_HH[,1]))
  {
    #Hourly and daily streamflows for r-th site
    streamflows_HH_r <- streamflows_HH[streamflows_HH[,1]%in%r,]
    
    `%>%` <- magrittr::`%>%`
    streamflows_GG_r <- streamflows_HH_r %>% 
      dplyr::group_by(Date=substring(Datetime,1,10)) %>% 
      dplyr::summarize(Value = mean(Value, na.rm=TRUE)) %>%
      as.data.frame()
    
    #r-th site PORFDC 
    Q<-unique(sort(streamflows_GG_r[,2],decreasing = TRUE))
    porfdc_r <- cbind.data.frame(p=c(1:length(Q)/(length(Q)+1)),Q)
    
    #Initialization of the dataframe used by the fitting model:
    # x: daily streamflow durations
    # y: maximum hourly streamflow/daily streamflow
    data <- data.frame(matrix(NA,nrow=length(streamflows_GG_r[,1]),ncol=2,dimnames = list(x=NULL,y=c("x","y"))),check.names = FALSE)
    
    for(Q in unique(streamflows_GG_r[,2])) data[streamflows_GG_r[,2]%in%Q,1] <- porfdc_r[porfdc_r[,2]%in%Q,1] #daily streamflow durations
    
    max_HH_daily<-streamflows_HH_r %>% #maximum hourly streamflow
      dplyr::group_by(Name = Name, Date = substring(Datetime,1,10)) %>%
      dplyr::summarise(Value = max(Value,nr.rm=T),.groups="drop") %>%
      as.data.frame()

    data[,2]<-max_HH_daily[,3]/streamflows_GG_r[,2] #maximum hourly streamflow/daily streamflow

    #Only durations below the threshold, corresponding to ratios between daily and long term streamflows above 1 
    data <- data[data<=d_threshold,]
    
    data<- data[order(data[,1],na.last=NA),] #Dataframe sorting according to durations
    
    #Moving Median: 30-day step
    smooth_median <- caTools::runquantile(data[,2],30,probs = 0.5,type=7,endrule="NA")#We could try a spinline?
    smooth_median <- data.frame(x=data[,1],y=smooth_median)
    smooth_median <- smooth_median[-which(is.na(smooth_median[,2])),]
    
    x<-smooth_median[,1]
    y<-smooth_median[,2]
    #Exponential model fitting: y = a+b*exp(-c*x)
    exp_model <- optim(1,exponential_model_1,method="Brent",upper=2000,lower=-100)
    
    #Salvataggio dei parametri del sito i-esimo
    qq_d_pars$c[qq_d_pars$Name%in%r] <- exp_model$par
    qq_d_pars$b[qq_d_pars$Name%in%r] <- (smooth_median$y[1] - smooth_median$y[length(smooth_median[,1])])/(exp(-exp_model$par*smooth_median$x[1])-exp(-exp_model$par*smooth_median$x[length(smooth_median$x)])) 
    qq_d_pars$a[qq_d_pars$Name%in%r] <- smooth_median$y[1] - qq_d_pars$b[qq_d_pars$Name%in%r]*exp(-smooth_median$x[1]*exp_model$par)
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
  
  qq_d_pars<-qq_d_pars[-which(qq_d_pars$Name%in%"Reno@Pracchia"),]
  qq_d_pars<-qq_d_pars[-which(qq_d_pars$Name%in%"Reno@CasalecchioChiusa"),]
  Pippo<-ToC[-which(names(ToC)%in%c("Reno@Pracchia","Reno@CasalecchioChiusa"))]
  plot(Pippo,qq_d_pars$b,pch=16)
  
  
  #Initialization of the exponential model parameters dataframe for fitting a=f(Tc), b=f(Tc), c=f(Tc) relationships
  abc_Tc_pars <- data.frame(matrix(NA,nrow=3,ncol=5,dimnames = list(x=NULL,y=c("Parameter","Formula","P1","P2","P3"))))
  abc_Tc_pars$Parameter <- c("a","b","c")
  
  #a=f(Tc)
  y <- qq_d_pars$a 
  x <- as.numeric(tc)
  exp_model<-optim(c(1,1),exponential_model_2o3)#2 parameters, sufficient to characterize the trend
  if(length(exp_model$par)==2) 
  {
    abc_Tc_pars$Formula[1]<-"y=P1*exp(-P2*x)"
    abc_Tc_pars[1,3:4]<-exp_model$par
    
  }else
  {
    abc_Tc_pars$Formula[1]<-"y=P1+P2*exp(-P3*x)"
    abc_Tc_pars[1,3:5]<-exp_model$par
  }  
  
  #b=f(Tc)  
  y <- qq_d_pars$b  
  x <- as.numeric(tc)
  exp_model<-optim(c(1,1,1),exponential_model_2o3)#3 parameters
  if(length(exp_model$par)==2) 
  {
    abc_Tc_pars$Formula[2]<-"y=P1*exp(-P2*x)"
    abc_Tc_pars[2,3:4]<-exp_model$par
    
  }else
  {
    abc_Tc_pars$Formula[2]<-"y=P1+P2*exp(-P3*x)"
    abc_Tc_pars[2,3:5]<-exp_model$par
  } 
  
  #c=f(Tc)  
  y <- log(qq_d_pars$c)#logarithmic transformation (more variance of c for the gauged section)
  x <- as.numeric(tc)
  exp_model<-optim(c(1,1),exponential_model_2o3)#2 parameters, sufficient to characterize the trend
  if(length(exp_model$par)==2) 
  {
    abc_Tc_pars$Formula[3]<-"y=P1*exp(-P2*x)"
    abc_Tc_pars[3,3:4]<-exp_model$par
    
  }else
  {
    abc_Tc_pars$Formula[3]<-"y=P1+P2*exp(-P3*x)"
    abc_Tc_pars[3,3:5]<-exp_model$par
  }
  
  # Store regional daily streamflows inside a list
  .internal_env$reg_streamflows[["abc_Tc"]] <- abc_Tc_pars
}
