#' @export
#2 reviews based on:
# - Climate elasticity of streamflows
#Visto il clima caldo temperato dell'area di indagine, valori tipici di elasticita' sono positivi e non
# superiori a 5. Si rimuovono quindi gli anni delle stazioni che non rispettano questa regola.
# - PORFDC trends

load("data/annual_precipitations.rda")
LowerElasticity<-0
UpperElasticity<-5

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
streamflow_reviews<-function(mean_annual_precipitation,LowerElasticity,UpperElasticity) 
{
  # - Climate elasticity of streamflows
  streamflows_GG<-use_internal_data(list_name="streamflows",dataset_name="GG")
  streamflows_HH<-use_internal_data(list_name="streamflows",dataset_name="HH")

  `%>%` <- magrittr::`%>%`
  mean_annual_streamflows<-streamflows_GG %>%
    dplyr::group_by(Name =Name, Year = lubridate::year(Date)) %>%
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
  
  #Unico dataset con MAP, MAS e elasticity
  #
  loop=T
  attempt_1<-list(MAS=MAS,MAP=MAP)
  attempt_2<-list(MAS=MAS,MAP=MAP)
  while(loop)
  {
    attempt<-elasticity_eval(attempt_2$MAS,attempt_2$MAP,LowerElasticity,UpperElasticity)
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
  streamflows_HH <-streamflows_HH[which(paste(streamflows_HH[,1],substring(streamflows_HH[,2],1,10),sep=" ")
                                  %in%paste(streamflows_GG[,1],streamflows_GG[,2],sep=" ")),]
  par(mfrow = c(1, 2),mar=c(4.2,4.2,3,1),xpd=F)  # 1 row, 2 columns

  log_streamflows_GG<-streamflows_GG
  log_streamflows_GG[,3]<-log(streamflows_GG[,3])
  if(length(which(log_streamflows_GG[,3]==-Inf))>0) log_streamflows_GG<-log_streamflows_GG[-which(log_streamflows_GG[,3]==-Inf),]
  
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
    
  
    plot(p_0,Q,ylab="log(Q)",type="l",xlab="p",col="black",lwd=2,cex.axis=1.3,cex.lab=1.5)
    grid()
    lines(p_0,Q,col="black",lwd=2)
    
    #plot(p_0[-1],dy_dx,type="l",xlab="p",ylab="Q'",lwd=2,cex.axis=1.5,cex.lab=1.5,main=r,cex.main=2)
    plot(p_0[-1:-2],d2y_d2x,type="l",xlab="p",ylab="log(Q)''",lwd=2,cex.axis=1.3,cex.lab=1.5)
    grid()
    lines(p_0[-1:-2],d2y_d2x,col="black",lwd=2)
    title(r,cex=3)
    change_sign[change_sign[,1]%in%r,2]<-sum(diff(sign(d2y_d2x)) != 0, na.rm = TRUE)
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

  streamflows_GG<-streamflows_GG[-which(streamflows_GG$Name%in%"Rubicone@Savignano"),]
  .internal_env$streamflows[["GG"]] <- streamflows_GG
  .internal_env$streamflows[["HH"]] <- streamflows_HH
  .internal_env$streamflows[["LT"]] <- rbind.data.frame(LT_MAS,.internal_env$streamflows[["LT"]])
}

library(ggplot2) #Per realizzare il plot
library(pals) #Per utilizzare la palette di colori kelly() 
install.packages("ggExtra")
library(ggExtra)
library(dplyr)
#Plot delle Elasticita', a valle della Revisione per Elasticita'
colore <- kelly() #stessa palette di colori impiegata al PUNTO 1
colore <- colore[c(3:7,2,8:10,12,16,19,20,11,21)]
colore <- c("black",colore[c(-1,-2,-13,-14,-15,-17,-18,-22)])#Si rimuovono colori che possono confondersi con quelli scelti, fino ad averne 15

unique_chars<-unique(Elasticity$Name)
codes<-c(10,28,6,8,25,15,24,12,2,32,30,3,7,9,5)
num_vec<-codes
char_to_num <- setNames(num_vec[1:length(unique_chars)], unique_chars)
result <- paste0(char_to_num[Elasticity$Name],". ",Elasticity$Name)
Elasticity$Name<-result
Elasticity<-Elasticity %>% arrange(as.numeric(sub("\\..*","",Elasticity$Name)))

A <- ggplot(Elasticity,aes(x = factor(Year), y = Value, fill=Name))+#Si sfrutta il dataframe creato: l'Anno come ascissa, l'Elasticita' come ordinata, il Colore specifico per ciascuna Sezione
  labs(y = expression(paste(epsilon," (-)")), x=" ")+ #Titoli degli assi 
  scale_y_continuous(breaks = seq(from=-5,to=10,by=1))+
  theme(axis.title.x = element_text(size=15, vjust = -6),
        axis.title.y = element_text(size=18,vjust = 8),
        axis.text.x = element_text(size=15,angle = 90, vjust = .5),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.margin = margin(1,1,1,1.1, "cm"),
        legend.position = "top")+
  guides(fill = guide_legend(direction = "horizontal", nrow = 4))+
  coord_cartesian(ylim=c(-5,10))+
  removeGrid(x=TRUE,y=TRUE)+
  geom_bar(stat = "identity", position = "dodge",na.rm=TRUE)+ 
  geom_hline(yintercept=0,linetype="solid",col="grey")+
  geom_vline(xintercept = seq(from=1.5,to=20.5,by=1),linetype="solid",col="grey")+
  scale_fill_manual(NULL,values=colore,limits=c(unique(Elasticity$Name)))#Si utilizzano i colori scelti per realizzare il riempimento degli istogrammi relativi a ciascuna delle 15 sezioni
A









