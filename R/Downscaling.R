#################################################################################
#                       DOWNSCALING DELLE PORTATE GIORNALIERE                   #
#                                                                               #
# A partire dall'idrogramma sintetico giornaliero sulla sezione target, si      #
# convertono in portate orarie secondo due modalita':                           #
# - Interpolazione lineare per d>d_asterisco                                    #                                                                  
# - Funzione polinomiale per d<=d_asterisco (piene), con colmo stimato tramite  #
# relazione esponenziale negativa Qcolmo/Qgg=f(d) i cui parametri funzione di Tc#     
# d_asterisco = minima durata tale che Q/muQ>1 per i siti strumentati           #  
#################################################################################

dev.off() #Pulizia Plots
rm(list=ls()) #Pulizia Environment
cat("\014") #Pulizia Console

install.packages("DEoptim")
library(tools)#funzione file_ext() per estrarre uno o piu' file con la stessa estensione
library(lubridate) #per lavorare con date-ore
library(stringr) #per inserimento delle ore 
library(pracma) #pacchetto per il calcolo integrale (funzione "Idrogramma_hh.r")
library(matlib) #per risolvere sistemi di equazioni lineari (funzione "Idrogramma_hh.r")
library(ggplot2)#per plot
library(scales) #per assi log nel plot
library(DEoptim)

#IMPOSTARE "Algoritmo di generazione" COME CARTELLA DI LAVORO
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..");setwd("..")

#Funzione per eliminare ridonandanze di portata nella curva di durata
source("04_Metodo_Deflusso_Indice/04_2_FDC_dimensionale/Remove_redundancy_FDC.R")
#Funzione per il calcolo della funzione polinomiale delle portate orarie sintetiche (d<=d_asterisco)
source("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/Idrogramma_hh.R")
'%!in%'=Negate('%in%')#Negazione "%in%"

#Ai fini della Convenzione, si e' scelto:
#- Savio a San Carlo e Marecchia a Rimini SS16 (2009-2019) donatori per serie CALIBRAZIONE del Modello idro-morfodinamico
#- Reno a Casalecchio Chiusa (1921-2021), donatore per costruzione serie SCENARI IDROCLIMATICI (TIPICO e GRAVOSO, vedi "Scelta_scenari.r")

#Comandi per importare le serie di portata giornaliera:
Calibrazione   <- F #Calibrazione (donatori Savio a San Carlo e Marecchia a Rimini SS16)
#Se non e' TRUE: donatore Reno a Casalecchio Chiusa

tau_Tc_pars    <- read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/Parametri_tauTc.csv",sep=";",dec=",")
lambdas_Tc_pars <- read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/Parametri_lambdaTc.csv",sep=";",dec=",")

#Caricamento delle relazioni dei parametri, della funzione esponenziale negativa Qcolmo/Qgg=f(Durata), in funzione del tempo di corrivazione
Par_Tc_FUN <- read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/Parametri_abcTc.csv",sep=";",dec=",")
#Caricamento dei parametri della relazione esponenziale negativa e dei tempi di corrivazione per le sezioni impiegate nella stima dei parametri
Par_Tc <- read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/abc_sezioni_chiusura.csv",sep=";",dec=",")
#Caricamento dei tempi di corrivazione per ciascuna sezione delle sezioni non strumentate
Tc_ungauged<-read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/Tc_INMarecchia.csv",sep=";",dec=",")
#Unione dei tempi di corrivazione
Tc<-rbind.data.frame(Tc_ungauged,rep(NA,3));Tc$Code[9]<-8;Tc$Name[9]<-"SS16";Tc$Value[9]<-Par_Tc$Tc[which(Par_Tc$Code%in%8)]
remove(Tc_ungauged)

#Nomi file CSV delle serie sintetiche di portata giornaliera, ordinati per aree crescenti (i nomi sono gli stessi, qualunque sia il sito donatore)
Names <- dir("04_Metodo_Deflusso_Indice/04_3_Idrogramma_giornaliero/Qest_GG/SavioSCarlo_donatore")#Si sceglie "/RenoChiusa_donatore"
Names <- Names[which(file_ext(Names)%in%"csv")]#Si considerano solo i CSV
Names<-Names[order(match(sub("\\..*","",sub("^[^_]*_","",sub("^[^_]*_", "",Names))),Tc$Name))]#Ordinamento, dalla sezione piu' a monte

#Caricamento delle modifiche da apportare ai parametri dell'esponenziale negativo, in base al confronto tra massimi annuali di portata oraria sintetica
#e la curva di frequenza (WP3+Topkriging) di ciascun sito (v. ".../04_4_Idrogramma_orario/Revisione_parametri_WP3")
if(Calibrazione) Par_modifiche<-read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Revisione_Parametri_WP3/Modifiche_Parametri_Calibrazione_cfrWP3.csv",sep=";",dec=",")#Calibrazione  
if(!Calibrazione) Par_modifiche<-read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Revisione_Parametri_WP3/Modifiche_Parametri_RenoChiusa_cfrWP3.csv",sep=";",dec=",")#Scenari

#Caricamento d_asterisco
d_asterisco <- as.numeric(read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/Durata_soglia_Downscaling.csv",sep=";",dec=",",check.names = F)) 
λ_V <- .5
λ_Q <- .5
V_threshold<-.25
max_iter <- 20

#Posizione del picco: idrogramma orario Savio a S.Carlo
Qobs_hh_don<-read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/CSV_HH/2002-2021_SavioSCarlo.csv",sep=";",dec=",")

#Calcolo del rapporto Qcolmo/Qgg per la durata d_asterisco sulla sezione target "SS16"
# Index_target <- which(Par_Tc$Name=="MarecchiaSS16")
# Par_Tc$a[Index_target]+Par_Tc$b[Index_target]*exp(-Par_Tc$c[Index_target]*d_asterisco)# = 1.08, quindi per durate >= d_asterisco si puo' tranquillamente confondere
# #la portata giornaliera con la massima portata oraria in quel giorno
durations<-c(.001,1:31)

for(isite in 1:length(Names))#Ciclo su ciasun sito non strumentato di interesse
{
  isite<-9
  #Caricamento della/e serie sintetica/che di portata media giornaliera stimata/e per il sito isite-esimo
  Qest_gg_all<-c()

  # if(Calibrazione)
  # {
  #   #Sito donatore Savio a S.Carlo
  #   Qest_gg_all<-append(Qest_gg_all,list(read.csv(paste0("04_Metodo_Deflusso_Indice/04_3_Idrogramma_giornaliero/Qest_GG/SavioSCarlo_donatore/",Names[isite]),sep=";",dec=",")))
  #   names(Qest_gg_all)<-"Savio"
  # 
  #   #Sito donatore Marecchia a Rimini SS16: per la sezione "SS16" si prendono le osservazioni disponibili
  #   if(isite!=9)
  #   {
  #     Qest_gg_all <- append(Qest_gg_all,list(read.csv(paste0("04_Metodo_Deflusso_Indice/04_3_Idrogramma_giornaliero/Qest_GG/MarecchiaSS16_donatore/",Names[isite]),sep=";",dec=",")))
  #     names(Qest_gg_all)[2]<-"Marecchia"
  #   }else
  #   {
  #     QMG <- read.csv("02_Riempimento&Revisione_QMG/02_2_Revisione_QMG/Dataset_QMG_definitivo.csv",sep=";",dec=",",check.names = F)
  #     QMG<-QMG[which(QMG$Name%in%"MarecchiaSS16"),]
  #     Qest_gg_all<-append(Qest_gg_all,list(data.frame(DATA=paste0(sort(rep(QMG$Year,365)),"-",substring(colnames(QMG)[-1:-3],4,5),"-",substring(colnames(QMG)[-1:-3],1,2)),VALORE=c(t(QMG[,-1:-3])))))
  #     names(Qest_gg_all)[2]<-"Marecchia"
  #   }
  # }else Qest_gg_all<-read.csv(paste0("04_Metodo_Deflusso_Indice/04_3_Idrogramma_giornaliero/Qest_GG/RenoChiusa_donatore/",Names[isite]),sep=";",dec=",")#Sito donatore Reno a Casalecchio Chiusa
  
  Qest_gg_all <- read.csv(paste0("04_Metodo_Deflusso_Indice/04_3_Idrogramma_giornaliero/Qest_GG/SavioSCarlo_donatore/",Names[isite]),sep=";",dec=",")

  #Curva di durata delle portate e raggruppamento anni consecutivi (infatti, per il downscaling e' necessario il dato di portata del giorno antecedente e successivo)
  #Per Calibrazione (anni 2009-2019), 2 siti donatori -> 2 curve di durata; sapendo a priori che il periodo temporale e' interamente coperto dal dato di portata, non e' necessario un raggruppamento ma conoscere gli anni in cui abbiamo dato dell'uno o dell'altro sito donatore
  PORFDC_est<-c()
  if(Calibrazione)#Calibrazione: siti donatori Savio a S. Carlo e Marecchia a SS16   
  {
    Anni<-c()
    for(i in 1:2)
    {
      PORFDC_est_i <- data.frame(p=c(1:length(Qest_gg_all[[i]]$VALORE)/(length(Qest_gg_all[[i]]$VALORE)+1)),Q=sort(Qest_gg_all[[i]]$VALORE,decreasing = TRUE))
      if(length(unique(PORFDC_est_i[,2]))!=length(PORFDC_est_i[,2])) PORFDC_est_i <- Remove_redundancy_FDC(PORFDC_est_i)
      PORFDC_est <- append(PORFDC_est,list(PORFDC_est_i))
      Anni<-append(Anni,list(unique(year(Qest_gg_all[[i]]$DATA))))
    }
    names(PORFDC_est)<-c("Savio","Marecchia")
    names(Anni)<-c("Savio","Marecchia")
    remove(PORFDC_est_i)
    Group_Anni<-0
  }else#Scenari: Sito donatore Reno a Chiusa di Casalecchio
  {
    Group_Anni<-c()
    PORFDC_est <- data.frame(p=c(1:length(Qest_gg_all$VALORE)/(length(Qest_gg_all$VALORE)+1)),Q=sort(Qest_gg_all$VALORE,decreasing = TRUE))
    if(length(unique(PORFDC_est[,2]))!=length(PORFDC_est[,2])) PORFDC_est <- Remove_redundancy_FDC(PORFDC_est)
    Group_Anni<-split(as.numeric(unique(substring(Qest_gg_all$DATA,first=1,last=4))),cumsum(c(1, diff(as.numeric(unique(substring(Qest_gg_all$DATA,first=1,last=4))))!= 1)))
  }
  
  #Calcolo dei parametri dell'esponenziale negativo a,b e c per la stima di Qcolmo, noto il tempo di corrivazione (Tc)
  x<-Tc$Value[which(Tc$Name==sub("^[^_]*_[^_]*_","",sub("\\..*", "",Names))[isite])]#Tc
  Par_target <- c()#Inizializz. vettore dei parametri
  for(i in 1:length(Par_Tc_FUN[,1]))#Calcolo dei parametri a, b e c, in base alle relazioni costruite per gli stessi 
  {
    if(length(which(is.na(Par_Tc_FUN[i,])==T))==0) {P1<-Par_Tc_FUN$P1[i];P2<-Par_Tc_FUN$P2[i];P3<-Par_Tc_FUN$P3[i]} else {P1<-Par_Tc_FUN$P1[i];P2<-Par_Tc_FUN$P2[i]}
    Par_target <- append(Par_target,eval(parse(text=substring(Par_Tc_FUN$Formula[i],first=3,last=nchar(Par_Tc_FUN$Formula[i])))))
  }
  
  #Modifiche a posteriori dei tre parametri dell'esponenziale negativo in base al confronto dei massimi annuali sintetici con curve di frequenza WP3
  #Par_target[1]<-Par_target[1]*Par_modifiche$a[which(Par_modifiche$Denominazione%in%sub("^[^_]*_[^_]*_","",sub("\\..*", "",Names))[isite])]
  #Par_target[2]<-Par_target[2]*Par_modifiche$b[which(Par_modifiche$Denominazione%in%sub("^[^_]*_[^_]*_","",sub("\\..*", "",Names))[isite])]
  #Par_target[3]<-Par_target[3]*Par_modifiche$c[which(Par_modifiche$Denominazione%in%sub("^[^_]*_[^_]*_","",sub("\\..*", "",Names))[isite])]
  
  lambda1<-lambdas_Tc_pars$q[1]+lambdas_Tc_pars$m[1]*x
  lambda2<-lambdas_Tc_pars$q[2]+lambdas_Tc_pars$m[2]*x
  
  tau    <-tau_Tc_pars$q+tau_Tc_pars$m*x
  
  norm_durations<-durations/x
  rd<-.5*exp(-lambda1*norm_durations/(1+lambda2*norm_durations))
  epsd<-tau/norm_durations*(1-exp(-norm_durations/tau))  
  drd_dd<--((lambda1/x)/(1+lambda2*norm_durations)^2)*rd
  depsd_dd<-1/norm_durations*(-tau/(x*norm_durations)*(1-exp(-norm_durations/tau))+(1/x*exp(-norm_durations/tau)))
  
  ##Disaggregazione temporale delle serie sintetiche di portata giornaliera - DOWNSCALING
  if(isite==9) Data_piene<-c()#Inizializzazione vettore date degli eventi di piena (d<=d_asterisco e portata giornaliera maggiore di quella dei giorni contigui); servira' per rappresentare i diagrammi di dispersione sugli eventi di piena per il sito target SS16
  
  lower <- -5000
  upper <- 5000
  c_vals <- seq(lower, upper, by = 10)
  for(l in 1:length(Group_Anni))#Ciclo su ciasun raggruppamento di anni
  {
    #Serie di portata giornaliera sintetica per il raggruppamento di anni l-esimo, cui si associa la durata (in base alla curva di durata delle portate corrispondente) 
    if(Calibrazione)#Calibrazione: siti donatori Savio a S. Carlo e Marecchia a SS16
    {
      #Periodo di calibrazione: 2009-2019 (buffer di un anno, per prendere anche il dato giornaliero degli anni antecedente e successivo)
      DATA=seq.Date(from=as.Date("2008-01-01","%Y-%m-%d"),to=as.Date("2020-12-31","%Y-%m-%d"),by="1 day")
      DATA=DATA[-which(substring(DATA,6,10)%in%"02-29")]
      Qest_gg <- data.frame(DATA,VALORE=rep(NA,length(DATA)))
      
      Anni[["Marecchia"]]<-Anni[["Marecchia"]][which(Anni[["Marecchia"]]>=min(year(Qest_gg$DATA))&Anni[["Marecchia"]]<=max(year(Qest_gg$DATA)))]
      Anni[["Savio"]]<-Anni[["Savio"]][which(Anni[["Savio"]]>=min(year(Qest_gg$DATA))&Anni[["Savio"]]<=max(year(Qest_gg$DATA)) & Anni[["Savio"]]%!in%Anni[["Marecchia"]])]
      
      for(i in 1:length(Anni[["Marecchia"]])) Qest_gg$VALORE[which(year(Qest_gg$DATA)%in%Anni[["Marecchia"]][i])]<-Qest_gg_all[["Marecchia"]]$VALORE[which(year(Qest_gg_all[["Marecchia"]]$DATA)%in%Anni[["Marecchia"]][i])]
      for(i in 1:length(Anni[["Savio"]])) Qest_gg$VALORE[which(year(Qest_gg$DATA)%in%Anni[["Savio"]][i])]<-Qest_gg_all[["Savio"]]$VALORE[which(year(Qest_gg_all[["Savio"]]$DATA)%in%Anni[["Savio"]][i])]
      
      #Buffer piu' fine sulla serie giornaliera: interesse sul 2008 per il 31/12 e sul 2020 per il 01/01
      Qest_gg<-Qest_gg[which(Qest_gg$DATA>=as.Date("2008-12-31","%Y-%m-%d") & Qest_gg$DATA<=as.Date("2020-01-01","%Y-%m-%d")),]
      
      Qest_gg <- cbind.data.frame(Qest_gg,p=rep(NA,length(Qest_gg[,1])))
      for(i in 1:length(Qest_gg[,1])) 
      {
        if(length(which(year(Qest_gg$DATA[i])%in%Anni[["Marecchia"]]))==1) Qest_gg$p[i] <- PORFDC_est[["Marecchia"]]$p[which(PORFDC_est[["Marecchia"]]$Q==Qest_gg$VALORE[i])]#durata  
        if(length(which(year(Qest_gg$DATA[i])%in%Anni[["Savio"]]))==1) Qest_gg$p[i] <- PORFDC_est[["Savio"]]$p[which(PORFDC_est[["Savio"]]$Q==Qest_gg$VALORE[i])]#durata
      }
    }else#Scenari: sito donatore Reno a Chiusa di Casalecchio
    {
      Qest_gg<-Qest_gg_all[which(as.numeric(substring(Qest_gg_all$DATA,first=1,last=4))%in%Group_Anni[[l]]),]#Serie l-esima
      rownames(Qest_gg)<-NULL
      Qest_gg <- cbind.data.frame(Qest_gg,p=rep(NA,length(Qest_gg[,1])))
      for (i in 1:length(Qest_gg[,1])) Qest_gg$p[i] <- PORFDC_est$p[which(PORFDC_est$Q==Qest_gg$VALORE[i])] #durata
    }
    
    Qest_gg$DATA <- as.character(Qest_gg$DATA)#Conversione a carattere della data (altrimenti convertita a numero nel dataframe "Qest_hh")
    Qest_hh<-rep(NA,length(Qest_gg$DATA)*24-12-11)#Serie oraria
    
    #Interpolazione lineare della serie giornaliera (intero idrogramma)
    for(i in 1:(length(Qest_gg$DATA)-1)) Qest_hh[((i-1)*24+1):(i*24+1)] <- approx(x=i:(i+1),y=Qest_gg$VALORE[c(i,i+1)],xout=seq(i,i+1,1/24),method="linear")$y
    
    #Si associa l'ascissa temporale all'idrogramma orario
    Qest_hh <- data.frame(DATA=rep(NA,length(Qest_hh)),ORA=rep(NA,length(Qest_hh)),VALORE=Qest_hh)
    
    #Inserimento dei giorni
    Qest_hh$DATA[1:12]<-rep(Qest_gg$DATA[1],12);Qest_hh$DATA[(length(Qest_hh[,1])-12):length(Qest_hh[,1])]<-rep(Qest_gg$DATA[length(Qest_gg$DATA)],13)#Primo giorno (inizia alle 12) e ultimo giorno (finisce alle 12) 
    for (i in 1:(length(Qest_gg[,1])-2)) Qest_hh$DATA[((i-1)*24+13):((i-1)*24+36)] <- rep(Qest_gg$DATA[i+1],24)#Restanti giorni 
    
    #Inserimento delle ore
    Qest_hh$ORA[1:12]<-paste0(12:23,":00:00");Qest_hh$ORA[(length(Qest_hh[,1])-12):length(Qest_hh[,1])]<-paste0(str_pad(0:12,2, pad = "0"),":00:00")#Ore del primo e ultimo giorno
    Qest_hh$ORA[13:(length(Qest_hh[,1])-13)] <- rep(paste0(str_pad(0:23,2, pad = "0"),":00:00"),length(Qest_gg[,1])-2)#Restanti ore
    
    Qest_gg$DATA <- as.Date(Qest_gg$DATA,"%Y-%m-%d")#Riconversione a data
    Qest_hh$DATA <- as.Date(Qest_hh$DATA,"%Y-%m-%d")#Conversione a data
    
    #Eventi di piena (d<=d_astersico e portata giornaliera maggiore di quella dei giorni contigui): si applica l'esponenziale negativo per ricavare Qcolmo del giorno dell'evento, fissata alle 12, 
    #e si applicano due polinomi del tipo: y=a*t^3+b*t^2+c*t+ (dalle 12 del giorno precedente alle 12 del giorno dell'evento e dalle 12 del giorno dell'evento alle 12 del giorno successivo); i parametri dei due polinomi
    #sono ricavati imponendo il passaggio per la Qcolmo e per le portate giornaliere dei giorni antecedente e successivo (fissate sempre alle 12), il massimo in corrispondenza di Qcolmo e il rispetto dei volumi d'acqua tra idrogramma orario e giornaliero
    #691-693,698

    min_DeltaV<-NULL
    c_50<-NULL
    c_min_DeltaV<-NULL

    
    time<-system.time(for(i in 1:length(Qest_gg$DATA))#Ciclo su ciasun dato della serie sintetica di portata giornaliera
    {
      Qmax<-(Par_target[1]+Par_target[2]*exp(-Par_target[3]*Qest_gg$p[i]))*Qest_gg$VALORE[i]
      Qd<-epsd*Qmax
      dQd_dd<-depsd_dd*Qmax
      
      #CASI PARTICOLARI: eventi di prima per giorni estremi della serie
      #1) evento di piena il primo giorno della serie (25 ore simulate)
      # if (i==1)#Prima giorno
      # {
      #   if(Qest_gg$VALORE[i]>Qest_gg$VALORE[i+1] && Qest_gg$p[i]<=d_asterisco)
      #   {
      #     if(isite==9) Data_piene<-append(Data_piene,Qest_gg$DATA[i])
      #     
      #     pos_max<-which.max(Qobs_hh_don$VALORE[Qobs_hh_don$DATA%in%Qest_gg$DATA[i]])
      #     
      #     Volume_gg<-12*Qest_gg$VALORE[i]+12*Qest_gg$VALORE[i+1]#Volumi idrogramma giornaliero
      #     
      #     x<-c(1,25)#Ore per cui e' nota la portata oraria
      #     y<-c(Qmax,Qest_gg$VALORE[i+1])#Portate orarie note (Qcolmo tramite esponenziale e portata giornaliera)       
      #     
      #     #Calcolo dei parametri del polinomio di terzo grado, rappresentativo del ramo discendente dell'idrogramma di piena (1a-25a ora)  
      #     result<-"Calibration"#Affinche' "Idrogramma_hh" restituisca la differenza volumetrica tra idrogramma orario e giornaliero  
      #     k<-2#Indice di x e y relativo alla portata oraria nota oltre al colmo
      #     Delta_Vol <- optim(1,polynomial_1,method="Brent",lower=-10,upper=10)#Minima differenza volumetrica tra gli idrogrammi per d nel dominio [-10;10] (d parametro del polinomio y = a*x^3+b*x^2+c*x+d) 
      #     
      #     if(round(Delta_Vol$par)==10)#Raggiunto il limite superiore di d -> si modifica aggiungendo un'unita' fintanto che non si raggiunge una differenza volumetrica inferiore o pari al m3/s 
      #     {                           
      #       j=1
      #       while(Delta_Vol$value[1]>V_threshold)#Soglia: 10% su volumi e portate
      #       {
      #         Delta_Vol <- optim(1,polynomial_1,method="Brent",lower=-10,upper=10+j)
      #         j=j+1
      #       } 
      #     }
      #     
      #     if(round(Delta_Vol$par)==-10)#Raggiunto il limite inferiore di d -> si modifica togliendo un'unita' fintanto che non si raggiunge una differenza volumetrica inferiore o pari al m3/s 
      #     {
      #       j=1
      #       while(Delta_Vol$value[1]>V_threshold)#Soglia: 10% su volumi e portate
      #       {
      #         Delta_Vol <- optim(1,polynomial_1,method="Brent",lower=-10-j,upper=10)
      #         j=j+1
      #       } 
      #     }
      #     
      #     result<-"Series"#Affinche' "Idrogramma_hh" restituisca la serie di portata oraria dell'evento, con i parametri del polinomio calibrati
      #     Q_hh<-polynomial_1(Delta_Vol$par)
      #     
      #     Slope     <- (Qest_gg$VALORE[i+2]-Qest_gg$VALORE[i+1])/24#pendenza
      #     Intercept <-  Qest_gg$VALORE[i+1] - Slope*25#intercetta
      #     Q_hh_lin<- Slope*c(1:25)+Intercept
      #     if(length(which(Q_hh<Q_hh_lin))>0) Q_hh[which(Q_hh<Q_hh_lin)]<-Q_hh_lin[which(Q_hh_lin>Q_hh)]
      #     
      #     #Inserimento del nuovo idrogramma d'evento nella serie oraria ottenuta per interpolazione lineare
      #     Qest_hh$VALORE[which(Qest_hh$DATA%in%Qest_gg$DATA[i]&Qest_hh$ORA%in%"12:00:00"):
      #                      which(Qest_hh$DATA%in%Qest_gg$DATA[i+1]&Qest_hh$ORA%in%"12:00:00")]<-Q_hh        
      #   }
      # }
      # 
      # #2) evento di piena l'ultimo giorno della serie (25 ore simulate)
      # if (i==length(Qest_gg$DATA))#Ultimo giorno
      # {
      #   if(Qest_gg$VALORE[i]>Qest_gg$VALORE[i-1] && Qest_gg$p[i]<=d_asterisco)
      #   {
      #     if(isite==9) Data_piene<-append(Data_piene,Qest_gg$DATA[i])
      #     Volume_gg<-12*Qest_gg$VALORE[i-1]+12*Qest_gg$VALORE[i]#Volumi idrogramma giornaliero (24 ore per l'evento di piena che diventano 12 in "Idrogramma_HH", stessi comandi di CASI ORDINARI)
      #     x<-c(1,25)#Ore per cui e' nota la portata oraria
      #     y<-c(Qest_gg$VALORE[i-1],Qmax)#Portate orarie note (portata giornaliera e Qcolmo tramite esponenziale)            
      #     
      #     #Calcolo dei parametri del polinomio di terzo grado, rappresentativo del ramo ascendente dell'idrogramma di piena (1a-25a ora) 
      #     result<-"Calibration"#Affinche' "Idrogramma_hh" restituisca la differenza volumetrica tra idrogramma orario e giornaliero  
      #     k<-1#Indice di x e y relativo alla portata oraria nota oltre al colmo
      #     Delta_Vol <- optim(1,polynomial_1,method="Brent",lower=-10,upper=10)#Minima differenza volumetrica tra gli idrogrammi per d nel dominio [-10;10] (d parametro del polinomio y = a*x^3+b*x^2+c*x+d) 
      #     
      #     if(round(Delta_Vol$par)==10)#Raggiunto il limite superiore di d -> si modifica aggiungendo un'unita' fintanto che non si raggiunge una differenza volumetrica inferiore o pari al m3/s 
      #     {
      #       j=1
      #       while(Delta_Vol$value>V_threshold)#Soglia: 10% su volumi e portate
      #       {
      #         Delta_Vol <- optim(1,polynomial_1,method="Brent",lower=-10,upper=10+j)
      #         j=j+1
      #       } 
      #     }
      #     
      #     if(round(Delta_Vol$par)==-10)#Raggiunto il limite inferiore di d -> si modifica togliendo un'unita' fintanto che non si raggiunge una differenza volumetrica inferiore o pari al m3/s 
      #     {
      #       j=1
      #       while(Delta_Vol$value>V_threshold)#Soglia: 10% su volumi e portate
      #       {
      #         Delta_Vol <- optim(1,polynomial_1,method="Brent",lower=-10-j,upper=10)
      #         j=j+1
      #       } 
      #     }
      #     
      #     result<-"Series"#Affinche' "Idrogramma_hh" restituisca la serie di portata oraria dell'evento, con i parametri del polinomio calibrati 
      #     Q_hh<-polynomial_1(Delta_Vol$par)
      #     
      #     #Linearizzazione delle portate
      #     Slope     <- (Qest_gg$VALORE[i-1] - Qest_gg$VALORE[i-2])/24#pendenza
      #     Intercept <- Qest_gg$VALORE[i-1] - Slope*1#intercetta
      #     Q_hh_lin<- Slope*c(1:25)+Intercept
      #     if(length(which(Q_hh<Q_hh_lin))>0) Q_hh[which(Q_hh<Q_hh_lin)]<-Q_hh_lin[which(Q_hh_lin>Q_hh)]
      #     
      #     #Inserimento del nuovo idrogramma d'evento nella serie oraria ottenuta per interpolazione lineare
      #     Qest_hh$VALORE[which(Qest_hh$DATA%in%Qest_gg$DATA[i-1]&Qest_hh$ORA%in%"12:00:00"):
      #                      which(Qest_hh$DATA%in%Qest_gg$DATA[i]&Qest_hh$ORA%in%"12:00:00")]<-Q_hh 
      #   }
      # }
      
      
      #CASI ORDINARI: evento di piena per giorno "interno" alla serie (49 ore simulate)
      if(i!=1 && i!=length(Qest_gg$DATA))
      {
        if(Qest_gg$VALORE[i]>Qest_gg$VALORE[i-1] && Qest_gg$VALORE[i]>Qest_gg$VALORE[i+1] && Qest_gg$p[i]<=d_asterisco)
        {
          
          #if(isite==9) Data_piene<-append(Data_piene,Qest_gg$DATA[i])
          pos_max<-25#+which.max(Qobs_hh_don$VALORE[Qobs_hh_don$DATA%in%Qest_gg$DATA[i]])
          x<-c(1,pos_max,49)#Ore per cui e' nota la portata oraria
          y<-c(Qest_gg$VALORE[i-1],Qmax,Qest_gg$VALORE[i+1])#Portate orarie note (portate giornaliere e Qcolmo tramite esponenziale)

          #Calcolo dei parametri dei polinomi di terzo grado, rappresentativi del ramo ascendente (1a-25a ora) e del ramo discendente (25a-49a ora) dell'idrogramma di piena  
          result<-"Calibration"#Affinche' "Idrogramma_hh" restituisca la differenza volumetrica tra idrogramma orario e giornaliero (valutata separatemente sui due rami)
          d_opt<-c() #Variabile per salvataggio dei d calibrati (d parametri dei polinomi y = a*x^3+b*x^2+c*x+d)  
          for(k in c(1,3))#Indice di x e y relativo alle portate orarie note oltre al colmo
          {
            if(k==1) Volume_gg<-12*Qest_gg$VALORE[i-1]+(pos_max-13)*Qest_gg$VALORE[i]#Volumi idrogramma giornaliero
            if(k==3) Volume_gg<-12*Qest_gg$VALORE[i+1]+(37-pos_max)*Qest_gg$VALORE[i]

            #Delta_Vol <- optim(1,polynomial_1,method="Brent",lower=-70,upper=4050)#Minima differenza volumetrica tra gli idrogrammi per d nel dominio [-10;10] (d parametro del polinomio y = a*x^3+b*x^2+c*x+d) 
            #Delta_Vol$value
            
            dV_above_threshold <- TRUE
            lower_bound <- -10
            upper_bound <- 10
            v0 <- 1
            cont <- 0
            
            while (dV_above_threshold && cont < max_iter)
            {
              cont<-cont+1
              Delta_Vol <- optimize(polynomial_1, lower=lower_bound, upper=upper_bound)
              
              if (Delta_Vol$objective < V_threshold) {
                dV_above_threshold <- FALSE
              } else 
              {
                par <- Delta_Vol$minimum
                # Expand bounds adaptively
                if (abs(par - upper_bound) <= abs(par - lower_bound)) {
                  lower_bound <- upper_bound
                  upper_bound <- upper_bound + abs(upper_bound) * 0.5
                } else {
                  upper_bound <- lower_bound
                  lower_bound <- lower_bound - abs(lower_bound) * 0.5
                }
                # New guess in the middle
                v0 <- (lower_bound + upper_bound) / 2
              }
            }
            
            # DEoptim(fn=polynomial_1,lower = -10000,upper = 1000000000)
            min_DeltaV<-append(min_DeltaV,Delta_Vol$objective)#append(min_DeltaV,Delta_Vol$optim$bestval)
            # errors <- sapply(c_vals,polynomial_1)
            
            # c_min_DeltaV<-append(c_min_DeltaV,c_vals[which.min(errors)])
            # c_50<-append(c_50,polynomial_1(50))
            # # plot(c_vals,errors,lwd=3,type = "l", main = "Delta vs c",
            # #      xlab = "c", ylab = "Delta V")
            # 
            # err_perc<-(c_50-min_DeltaV)/min_DeltaV*100

          
            # wrapper_function <- function(theta) {
            #   err <- compute_error(theta)
            #   if (err > threshold) return(Inf)
            #   return(err)
            # }
            
            
            
            
            # min(c_min_DeltaV)
            # max(c_min_DeltaV)
          #   if(round(Delta_Vol$par)==10)#Raggiunto il limite superiore di d -> si modifica aggiungendo un'unita' fintanto che non si raggiunge una differenza volumetrica inferiore o pari al m3/s 
          #   {
          #     j=1
          #     while(Delta_Vol$value>V_threshold)#Soglia: 1 m3; e' una soglia "severa" per compensare la linearizzione dei giorni precedente e successivo l'evento di piena (v. sotto)
          #     {
          #       Delta_Vol <- optim(1,polynomial_1,method="Brent",lower=-10,upper=10+j)
          #       j=j+1
          #     } 
          #   }
          #   
          #   if(round(Delta_Vol$par)==-10)#Raggiunto il limite inferiore di d -> si modifica togliendo un'unita' fintanto che non si raggiunge una differenza volumetrica inferiore o pari al m3/s
          #   {
          #     j=1
          #     while(Delta_Vol$value>V_threshold)#Soglia: 1 m3
          #     {
          #       Delta_Vol <- optim(1,polynomial_1,method="Brent",lower=-10-j,upper=10)
          #       j=j+1
          #     } 
          #   }
          #   d_opt <- append(d_opt,Delta_Vol$par)#Salvataggio d calibrati per i due rami dell'idrogramma di piena 
          }  
          # 
          # result<-"Series"#Affinche' "Idrogramma_hh" restituisca la serie di portata oraria dell'evento, con i parametri dei polinomi calibrati
          # k<-1; Q1<-polynomial_1(d_opt[1])#Ramo ascendente
          # k<-3; Q2<-polynomial_1(d_opt[2])#Ramo discendente
          # 
          # Q_hh <- c(Q1,Q2[-1])#Unione dei rami (ultima e prima portata oraria delle due serie coincidono)
          # 
          # #Linearizzazione delle portate
          # Slope     <- (Qest_gg$VALORE[i-1] - Qest_gg$VALORE[i-2])/24#pendenza
          # Intercept <- Qest_gg$VALORE[i-1] - Slope*1#intercetta
          # Q_hh_lin<- Slope*c(1:pos_max)+Intercept
          # if(length(which(Q_hh[1:pos_max]<Q_hh_lin))>0) Q_hh[which(Q_hh[1:pos_max]<Q_hh_lin)]<-Q_hh_lin[which(Q_hh_lin>Q_hh[1:pos_max])]
          # 
          # Slope     <- (Qest_gg$VALORE[i+2]-Qest_gg$VALORE[i+1])/24#pendenza
          # Intercept <-  Qest_gg$VALORE[i+1] - Slope*49#intercetta
          # Q_hh_lin<- Slope*c(pos_max:49)+Intercept
          # if(length(which(Q_hh[pos_max:49]<Q_hh_lin))>0) Q_hh[which(Q_hh[pos_max:49]<Q_hh_lin)+pos_max-1]<-Q_hh_lin[which(Q_hh_lin>Q_hh[pos_max:49])]
          # #plot(Q_hh,type="l")
          # #Inserimento del nuovo idrogramma d'evento nella serie oraria ottenuta per interpolazione lineare
          # Qest_hh$VALORE[which(Qest_hh$DATA%in%Qest_gg$DATA[i-1]&Qest_hh$ORA%in%"12:00:00"):
          #                  which(Qest_hh$DATA%in%Qest_gg$DATA[i+1]&Qest_hh$ORA%in%"12:00:00")]<-Q_hh
        }
      }
    })
    plot(c_50,type="l",lwd=2,col="red")
    lines(min_DeltaV,lwd=2,col="black")
    p_weibull<-1:length(min_DeltaV)/(length(min_DeltaV)+1)
    length(which(min_DeltaV<0.5))
    plot(sort(min_DeltaV,decreasing = T),p_weibull,col="black",ylab="p (-)",xlab="min_DeltaV (m3)",
         main=paste0("# iter:",max_iter,", V_threshold:",V_threshold,", time:",time[3]))
    max(min_DeltaV)
    median(min_DeltaV)
    if(l==1)#Primo Lancio
    { 
      Qest_hh_all<-Qest_hh #Inizializzione dataframe della serie oraria completa 
    }else #Lanci successivi al primo 
    {
      Qest_hh_all<-rbind.data.frame(Qest_hh_all,Qest_hh)#Si accodano le serie orarie al dataframe precedente inizializzato
    }
  }
  save(Data_piene,file="Data_piene.rda")
  plot(Q_hh,type="l")
  Q_maione<-SDH_Maione(3)
  lines(Q_maione$x+12,Q_maione$y,col="red")

  #Si verifica qualitativamente la corrispondenza fra idrogrammi orario e giornaliero
  if(Calibrazione)
  {
    t<-as.Date(Qest_hh_all$DATA,"%Y-%m-%d")+as.numeric(substring(Qest_hh_all$ORA,first=1,last=2))/24#Vettore tempi
    plot(t,Qest_hh_all$VALORE,type="l",col="red",lwd=2,ylim=c(min(Qest_hh_all$VALORE),max(Qest_hh_all$VALORE)),xlab="Anni",ylab=expression(paste("Portata ",m^3,"/s")),main=paste0(sub("\\..*","",sub("^[^_]*_[^_]*_","",Names[isite]))," - downscaling piene"))
    lines(Qest_gg[,c(1,2)],col="black",lwd=2,lty="dotted")
    legend("topleft",legend=c("Serie giornaliera","Serie oraria"),lwd=c(2,2),lty=c("dotted","solid"),col=c("black","red"))   
  }else
  {
    Qest_gg_all$DATA <- as.Date(Qest_gg_all$DATA,"%Y-%m-%d")
    t<-as.Date(Qest_hh_all$DATA,"%Y-%m-%d")+as.numeric(substring(Qest_hh_all$ORA,first=1,last=2))/24#Vettore tempi
    plot(t,Qest_hh_all$VALORE,type="l",col="red",lwd=2,ylim=c(min(Qest_hh_all$VALORE),max(Qest_hh_all$VALORE)),xlab="Anni",ylab=expression(paste("Portata ",m^3,"/s")),main=paste0(sub("\\..*","",sub("^[^_]*_[^_]*_","",Names[isite]))," - downscaling piene"))
    lines(Qest_gg_all[,c(1,2)],col="black",lwd=2,lty="dotted")
    legend("topleft",legend=c("Serie giornaliera","Serie oraria"),lwd=c(2,2),lty=c("dotted","solid"),col=c("black","red"))
  }
  
  #Salvataggio delle serie di portate oraria sintetica
  if(Calibrazione) write.table(Qest_hh_all,file=paste0("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/Qest_HH/Calibrazione/conMaione/",sub("gg","hh",Names[isite])),quote = FALSE,row.names = FALSE,col.names = TRUE,sep=";",dec=",")
  if(!Calibrazione) write.table(Qest_hh_all,file=paste0("04_Metodo_Deflusso_Indice/04_4_Idrogramma_oraria/Downscaling/Qest_HH/RenoChiusa_donatore/conMaione/",sub("gg","hh",Names[isite])),quote = FALSE,row.names = FALSE,col.names = TRUE,sep=";",dec=",")
}
plot(Qest_hh$VALORE[which(Qest_hh$DATA%in%c("2019-05-13","2019-05-12","2019-05-14"))],type="l")
min(min_Delta_V)
length(min_Delta_V)
##Idrogrammi orari (completo e di alcuni eventi di piena) e diagramma di dispersione (eventi di piena) per Rimini SS16:
# - serie sintetica (donatori Savio a S.Carlo e Marecchia a Rimini SS16, poiche' Reno a Casalecchio Chiusa ha poca sincronia con i deflussi di Rimini SS16)
# - serie osservata

Calibrazione<-NA#NA: per utilizzare la sola serie oraria con solo Savio a S.Carlo donatore 
if(Calibrazione) Qest_hh<-read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/Qest_HH/Calibrazione/conMaione/Qest_hh_SS16.csv",sep=";",dec=",")#Portate orarie sintetiche (Savio e Marecchia)
if(is.na(Calibrazione)) Qest_hh<-read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/Qest_HH/Calibrazione/Qest_hh_SS16_SavioSCarlo.csv",sep=";",dec=",")#Portate orarie sintetiche (Savio)
Qest_hh<-Qest_hh_all
Qobs_hh<-read.csv("04_Metodo_Deflusso_Indice/04_4_Idrogramma_orario/Downscaling/CSV_HH/2005-2020_MarecchiaSS16.csv",sep=";",dec=",")#Portate orarie osservate 
Qobs_hh<-Qobs_hh[-which(substring(Qobs_hh$DATA,first=6,last=10)%in%"02-29"),]#Rimozione di 29/02 (assenti, per costruzione, nella serie sintetica)
load("FDC2Qt_datasets/Modello AD/Streamflows/Q_hh_20_reg_SSQ_1813.rda")
Qsim_hh<-Qh_sim

#Accorpamento DATA+ORA in unico vettore (TEMPO)
Qest_hh$DATA<- parse_date_time(paste(Qest_hh$DATA,Qest_hh$ORA,sep=" "),orders="%Y-%m-%d %H:%M:%S",quiet=FALSE)
colnames(Qest_hh)[1]<-"TEMPO"
Qest_hh<-Qest_hh[,-2]
Qobs_hh$DATA<- parse_date_time(paste(Qobs_hh$DATA,Qobs_hh$ORA,sep=" "),orders="%Y-%m-%d %H:%M:%S",quiet=FALSE)
colnames(Qobs_hh)[1]<-"TEMPO"
Qobs_hh<-Qobs_hh[,-2]
colnames(Qsim_hh)<-c("TEMPO","VALORE")

t<-seq(min(Qest_hh$TEMPO,Qobs_hh$TEMPO,Qsim_hh$TEMPO),max(Qest_hh$TEMPO,Qobs_hh$TEMPO,Qsim_hh$TEMPO),by="hour")
#Vettore dei tempi 
if(Calibrazione) t<-seq(min(Qest_hh$TEMPO),max(Qest_hh$TEMPO),by="hour")
if(is.na(Calibrazione)) t<-seq(min(Qest_hh$TEMPO,Qobs_hh$TEMPO),max(Qest_hh$TEMPO,Qobs_hh$TEMPO),by="hour")

t<-t[-which(substring(t,first=6,last=10)%in%"02-29")]

#Inizializzazione dataframe per il plot degli idrogrammi osservati e simulati
df_plot<-data.frame(matrix(NA,length(t),4,dimnames=list(x=NULL,y=c("TEMPO","OBS","EST","SIM"))))

df_plot$TEMPO<-t#tempi
df_plot$OBS[which(df_plot$TEMPO%in%Qobs_hh$TEMPO)]<-Qobs_hh$VALORE[which(Qobs_hh$TEMPO%in%df_plot$TEMPO)]#portate osservate
df_plot$EST[which(df_plot$TEMPO%in%Qest_hh$TEMPO)]<-Qest_hh$VALORE[which(Qest_hh$TEMPO%in%df_plot$TEMPO)]#portate sintetiche
df_plot$SIM[which(df_plot$TEMPO%in%Qsim_hh$TEMPO)]<-Qsim_hh$VALORE[which(Qsim_hh$TEMPO%in%df_plot$TEMPO)]#portate simulate

Hydro_all<-F#Idrogramma completo
if(!Hydro_all)#Idrogramma di sngoli eventi di piena 
{
  #df_plot_05<-df_plot[which(substring(df_plot$TEMPO,1,19)%in%"2005-11-24"):which(substring(df_plot$TEMPO,1,19)%in%"2005-11-29 23:00:00"),]#26/11/2005
  #df_plot<-df_plot[which(substring(df_plot$TEMPO,1,19)%in%"2019-05-05"):which(substring(df_plot$TEMPO,1,19)%in%"2019-05-07 23:00:00"),]#06/05/2019
  df_plot_05<-df_plot[which(year(df_plot$TEMPO)%in%2005 & month(df_plot$TEMPO)%in%11),]#Maggio 2019 
}
#df_plot_05[,which(df_plot_05$OBS%in%max(df_plot_05$OBS))]
t<-df_plot_05$TEMPO
#Plot idrogrammi orari simulato e osservato
par(mar=c(4,4,2,1))
par(mgp=c(2.5,1,0))
plot(df_plot_05[,c(1,2)],col="forestgreen",lwd=3.5,type="l",
     xlab="",
     ylab=expression(paste("streamflow (",m^3,"/s)",sep="")),
     #xaxt="n",
     yaxt="n",
     main=if(Hydro_all|length(unique(day(df_plot_05$TEMPO)))>27) "Marecchia@Rimini SS16, 11/2005" else paste0("Marecchia@RiminiSS16, ",format(strptime(unique(substring(df_plot_05$TEMPO,6,10))[2],"%Y-%m-%d"),"%d/%m/%Y")),
     xlim=c(min(df_plot_05$TEMPO),max(df_plot_05$TEMPO)),
     ylim=c(0,max(c(df_plot_05$OBS,df_plot_05$EST),na.rm=T)),cex.lab=1.2,cex.axis=1.2)
#axis(1,at=t[c(1,which(month(t)==1 & day(t)==1 & hour(t)==0))],labels=t[],cex.axis=1.2)
axis(2,at=seq(0,max(c(Qest_hh$VALORE,Qobs_hh$VALORE),na.rm=T),100),labels=seq(0,max(c(Qest_hh$VALORE,Qobs_hh$VALORE),na.rm=T),100),cex.axis=1.1)
abline(h=seq(0,max(c(Qest_hh$VALORE,Qobs_hh$VALORE),na.rm=T),100),col="lightgrey",lwd=.001,lty="solid")
lines(df_plot[,c(1,3)],col="black",lwd=4,lty=3)
lines(df_plot[,c(1,4)],col="red",lwd=4,lty=if(Hydro_all) "dotted" else 3)
if(!Hydro_all) points(df_plot[,c(1,3)],pch=16,col="black")
legend("topleft",legend=c("observed","GR5H","regional"),lty=c(if(Hydro_all) "dotted" else "solid","dotted","dotted"),lwd=c(3.5,4,4),col=c("forestgreen","red","black"),pch=c(if(Hydro_all) NA else NA),cex=1.2,box.col="transparent",bg="transparent")

#Diagrammi di dispersione eventi di piena
df_plot<-df_plot[-which(is.na(df_plot$OBS)==T|is.na(df_plot$EST)==T),]#Taglio sui soli istanti comuni alle due serie
df_plot<-df_plot[which(strptime(substring(df_plot$TEMPO,1,10),"%Y-%m-%d")%in%Data_piene),]#Taglio sui soli eventi di piena
R2<-cor(df_plot$OBS,df_plot$EST,use="complete.obs")^2#0.579

dataframe<-data.frame(matrix(NA,nrow=length(df_plot$TEMPO),ncol=2))
A<-ggplot(dataframe,aes(x=df_plot$OBS,y=df_plot$EST))+
  geom_bin2d(bins = 50)+
  scale_fill_continuous(type="viridis",trans="reverse")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels=function(x) sprintf("%.1f", x),limits=c(0.05,max(c(df_plot$OBS,df_plot$EST),na.rm=T)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels=function(x) sprintf("%.1f", x),limits=c(0.05,max(c(df_plot$OBS,df_plot$EST),na.rm=T)))+
  geom_abline(intercept = 0, slope = 1,colour="red",lwd=1.4)+
   ggtitle("")+
  theme(axis.text=element_text(size=17),
        axis.title=element_text(size=17))+
  theme(plot.title = element_text(hjust = 0.5,color = "black", size =30, face = "bold"),legend.text = element_text(size=17),legend.title = element_text(size=17))+
  labs(x=expression(paste("observed (",m^3,"/s)",sep="")),y=expression(paste("regional (",m^3,"/s)",sep="")))+
  annotate("text",x=0.5,y=max(c(df_plot$OBS,df_plot$EST),na.rm=T)*0.95,label= bquote(R^2*"="*~ .(round(c(R2),2))*" (0.52)"),size=8)
A

save(Data_piene,file="Data_piene.rda")
