#' @export
#2 reviews based on:
# - Climate elasticity of streamflows
#Visto il clima caldo temperato dell'area di indagine, valori tipici di elasticita' sono positivi e non
# superiori a 5. Si rimuovono quindi gli anni delle stazioni che non rispettano questa regola.
# - PORFDC trends

#' Title
#'
#' @param mean_annual_precipitation 
#' @param LowerElasticity 
#' @param UpperElasticity 
#'
#' @return
#' @export
#'
#' @examples

check_streamflow_quality<-function(streamflows_GG,streamflows_HH,annual_precipitation,LowerElasticity,UpperElasticity) 
{
  ##check 1: precipitation elsticity of streamflows
  `%>%` <- magrittr::`%>%`
  
  #MAF evaluation
  sites_to_remove<-c()#removal of sites with less then 5 years
  for(r in unique(streamflows_GG[,1]))
  {
    if(length(unique(lubridate::year(streamflows_GG$Datetime[streamflows_GG$Name%in%r])))<5)
      sites_to_remove<-append(sites_to_remove,r)
  }
  if(length(sites_to_remove)>0) streamflows_GG<-streamflows_GG[-which(streamflows_GG$Name%in%sites_to_remove),]  
    
  mean_annual_streamflows<-streamflows_GG %>%
    dplyr::group_by(Name =Name, Year = lubridate::year(Datetime)) %>%
    dplyr::summarise(Value = mean(Value, na.rm = TRUE),.groups="drop") %>%
    as.data.frame()
  
  MAP<-annual_precipitations
  MAS<-mean_annual_streamflows
  
  #Remove NA values
  if(length(which(is.na(MAS[,3])))>0) MAS<-MAS[-which(is.na(MAS[,3])),]
  if(length(which(is.na(MAP[,3])))>0) MAP<-MAP[-which(is.na(MAP[,3])),]
  
  #Common sections
  MAS<-MAS[MAS[,1]%in%unique(MAP[,1]),]
  MAP<-MAP[MAP[,1]%in%unique(MAS[,1]),]

  #MAP sorting with the respect to MAS 
  MAP<-MAP[order(factor(MAP[,1],levels=unique(MAS[,1]))),]
  
  #Common years
  index_MAS<-NULL;index_MAP<-NULL
  for(r in unique(MAS[,1]))
  {
    index_MAS<-append(index_MAS,which(!MAS[MAS[,1]%in%r,2]%in%MAP[MAP[,1]%in%r,2])+which(MAS[,1]%in%r)[1]-1)
    index_MAP<-append(index_MAP,which(!MAP[MAP[,1]%in%r,2]%in%MAS[MAS[,1]%in%r,2])+which(MAP[,1]%in%r)[1]-1)
  }
  if(length(index_MAS)>0) MAS<-MAS[-index_MAS,]
  if(length(index_MAP)>0) MAP<-MAP[-index_MAP,]
  
  #Unique dataset with MAP, MAS and elasticity
  #
  loop=T
  attempt_1<-list(MAS=MAS,MAP=MAP)
  attempt_2<-list(MAS=MAS,MAP=MAP)
  while(loop)
  {
    attempt<-prec_elasticity(attempt_2$MAS,attempt_2$MAP,LowerElasticity,UpperElasticity)
    attempt_1<-attempt_2
    attempt_2<-attempt
    if(dim(attempt_1$MAS)[1]==dim(attempt_2$MAS)[1]) loop=F
  }
  LT_MAS<-attempt$LT_MAS
  Elasticity<-attempt$Elasticity
  
  index_to_keep<-c()
  for(r in unique(Elasticity[,1]))
  {
    index_to_keep<-append(index_to_keep,
                          which(lubridate::year(streamflows_GG[streamflows_GG[,1]%in%r,2])%in%
                                  Elasticity[Elasticity[,1]%in%r,2])+which(streamflows_GG[,1]%in%r)[1]-1)
  }
  
  streamflows_GG    <- streamflows_GG[index_to_keep,]
  streamflows_HH    <- streamflows_HH[which(paste(streamflows_HH[,1],substring(streamflows_HH[,2],1,10),sep=" ")
                                  %in%paste(streamflows_GG[,1],streamflows_GG[,2],sep=" ")),]
  
  
  
  ##check 2: number of inflection points in the logarithmic period of record

  log_streamflows_GG<-streamflows_GG
  log_streamflows_GG[,3]<-log(streamflows_GG[,3])
  if(length(which(log_streamflows_GG[,3]==-Inf))>0) log_streamflows_GG<-log_streamflows_GG[-which(log_streamflows_GG[,3]==-Inf),]
  
  #Definition of same number of points of each PORFDC
  p_inf<--Inf
  p_sup<-Inf
  N<--Inf
  for(r in unique(log_streamflows_GG[,1]))
  {
    str_r<-log_streamflows_GG[log_streamflows_GG[,1]%in%r,]
    p=c(1,length(str_r[,1]))/(length(str_r[,1])+1)
    if(p[1]>p_inf) p_inf<-p[1] 
    if(p[2]<p_sup) p_sup<-p[2] 
    if(length(str_r[,1])>N) N<-length(str_r[,1])
  }
  p_0<-1:N/(N+1)
  
  p_0<-p_0[which(p_0>p_inf & p_0<p_sup)]
  
  #evaluation of inflection points (change of sign in the second derivative)
  change_sign<-data.frame(Name=unique(log_streamflows_GG[,1]),Value=rep(NA,length(unique(log_streamflows_GG[,1]))))
  for(r in unique(log_streamflows_GG[,1]))
  {
    str_r<-log_streamflows_GG[log_streamflows_GG[,1]%in%r,]
    
    p=1:length(str_r[,1])/(length(str_r[,1])+1)
    Q=sort(log_streamflows_GG[log_streamflows_GG[,1]%in%r,3],decreasing = T)
    Q<-approx(x=p,y=Q,xout=p_0,method="linear")$y
    
    # Compute first derivative (slope)
    dy_dx <- diff(Q) / diff(p_0)
    
    #smooth_dy_dx<-caTools::runquantile(dy_dx,150,probs = 0.5,type=7,endrule="NA")
    # Compute second derivative (rate of slope change)
    d2y_d2x <- diff(dy_dx) / diff(p_0[-1]) 
  }
  index_cs<-which(change_sign[,2]>N/2)
  if(length(index_cs)>0) 
  {
    cat("Error:\nPORFDCs with too many inflection points\n",paste(change_sign[index_cs,1],collapse="\n"))
    ask_yes_no <- function(prompt = "Do you want to delete these sections?[Y/N]:") 
    {
      repeat {
        response <- toupper(trimws(readline(prompt)))
        
        if (response %in% c("Y", "N")) {
          return(response)
        } else {
          cat("Invalid input. Please enter 'Y' or 'N'.\n")
        }
      }
    }
    user_input<-ask_yes_no()
    if(user_input=="Y") 
    {  
      streamflows_GG<- streamflows_GG[-which(streamflows_GG[,1]%in%change_sign[index_cs,1]),]
      streamflows_HH<- streamflows_HH[-which(streamflows_HH[,1]%in%change_sign[index_cs,1]),]
      LT_MAS     <- LT_MAS[-index_cs,]
    }
  }
  
  # Store all dataframes inside a list
  streamflows_check <- list(HH=streamflows_HH,
                            GG=streamflows_GG,
                            LT=rbind.data.frame(LT_MAS,streamflows[["LT"]]))
  
  return(streamflows_check)
}