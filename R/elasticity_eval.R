#' Title
#'
#' @param MAS 
#' @param MAP 
#' @param LowerElasticity 
#' @param UpperElasticity 
#'
#' @return
#' @export
#'
#' @examples
elasticity_eval<-function(MAS,MAP,LowerElasticity,UpperElasticity)#input: mean annual precipitation, mean annual streamflows
{
  #Long term streamflows
  LT_MAS<-MAS %>%
    dplyr::group_by(Name = Name) %>%
    dplyr::summarise(Value = mean(Value, na.rm = TRUE),.groups="drop") %>%
    as.data.frame()
  
  #Long term precipitations
  LT_MAP<-MAP %>%
    dplyr::group_by(Name = Name) %>%
    dplyr::summarise(Value = mean(Value, na.rm = TRUE),.groups="drop") %>%
    as.data.frame()
  
  #Deviations evaluation
  dev_MAS<-NULL;dev_MAP<-NULL
  for(r in unique(MAS[,1]))
  {
    dev_MAS<-append(dev_MAS,MAS[MAS[,1]%in%r,3]/LT_MAS[LT_MAS[,1]%in%r,2]-1)
    dev_MAP<-append(dev_MAP,MAP[MAP[,1]%in%r,3]/LT_MAP[LT_MAP[,1]%in%r,2]-1)
  }
  
  #Elasticity evaluation
  index_to_keep<-NULL
  for(i in 1:length(dev_MAS))
  {
    if(abs(dev_MAP[i])<.1 && abs(dev_MAS[i])<.2)
      index_to_keep<-append(index_to_keep,i)
  }
  
  Elasticity<-MAS
  Elasticity[,3]<-dev_MAS/dev_MAP
  index_to_remove<-which(Elasticity[,3]<LowerElasticity | Elasticity[,3]>UpperElasticity)
  index_to_remove<-index_to_remove[!index_to_remove%in%index_to_keep]
  
  if(length(index_to_remove)>0)
  {
    MAS<-MAS[-index_to_remove,]
    MAP<-MAP[-index_to_remove,]
    Elasticity<-Elasticity[-index_to_remove,]
  }  
  
  output<-list(MAS=MAS,MAP=MAP,Elasticity=Elasticity,LT_MAS=LT_MAS)
  
  return(output)
}
