#Multiregression model
target_section<-"Marecchia@RiminiSS16"
load("data/basin_descriptors.rda")

reg_LT_streamflow<-function(basin_descriptors,target_section)
{
  LT_streamflows<-use_internal_data(list_name="streamflows",dataset_name="LT")#Long term streamflows
  #basin_descriptors<-use_internal_data(list_name="input_data",dataset_name="basin_descriptors")#Basin descriptors
  
  #Common and sorted sections
  LT_streamflows<-LT_streamflows[LT_streamflows[,1]%in%unique(basin_descriptors[,1]),]
  basin_descriptors<-basin_descriptors[basin_descriptors[,1]%in%unique(LT_streamflows[,1]),]
  basin_descriptors<-basin_descriptors[order(factor(basin_descriptors[,1],levels=unique(LT_streamflows[,1]))),]
  
  #Power relationship between long term streamflows and basin descriptors: logarithmic transformation in order to apply linear regression
  #Removing descriptors with at least one negative or NA value
  Index_c<-which(sapply(basin_descriptors,function(column) any(column<=0)))#Negative values
  Index_c<-append(Index_c,which(sapply(basin_descriptors,function(column) any(is.na(column)))))#NA values
  if(length(Index_c)>0) basin_descriptors <- basin_descriptors[,-Index_c] 

  best_descriptors<-function(streamflows,descriptors,N)
  {
    #Log transformation
    y=log(streamflows[,-1])
    x=log(descriptors[,-1])
    
    #R2 adjusted metrics temporary variable 
    R2_temp_1  <- 0  #Highest R2 adjusted of the k-th while cycle
    R2_temp_2  <- -1 #Highest R2 adjusted of the (k-1)-th while cycle
    #lm summary temporaty variables
    Model_best_sum_1 <- list() #Summary of the highest R2 adjusted of the k-th while cycle
    Model_best_sum_2 <- list() #Summaries of the highest R2 adjusted of the (k-1)-th while cycle
    #Model descriptors
    Best_descriptors <- matrix(,nrow=length(x[,1]),ncol=0,dimnames=list(x=descriptors[,1],y=NULL))
    
    while(R2_temp_1 > R2_temp_2 && length(Best_descriptors[1,])<N)#repeat while adding another descriptor increase R2 adjusted metrics 
    {
      R2_temp_2    <- R2_temp_1 #Highest adjusted R2 (best model) of the previous cycle
      R2_temp_1    <- 0 #Temporary variable update
      
      #Searching for a new descriptor between the residual ones
      for (i in 1:length(x[1,]))
      {
        if(R2_temp_2 == 0)#Enter in the for cycle for the first time
        {
          Model <- lm(y ~ x[,i])#Linear model 
        }else if(R2_temp_2 != 0)#Enter in the for cycle from the second time
        {
          Model <- lm(y ~ Best_descriptors + x[,i])#Multiple linear regression
        }
        Model_sum <- summary(Model) #Model summary
        
        if(is.na(Model_sum$adj.r.squared)) break #Esc while cycle in case R2 adjusted is NA 
        if(Model_sum$adj.r.squared > R2_temp_1)
        {
          R2_temp_1 <- Model_sum$adj.r.squared #Update R2 adjusted
          Model_best_sum_1 <- Model_sum #Update summary 
          Pos_best_descriptor <- i #index i-th descriptor (in order to save it in "Best_descriptors" and remove it from "x")
        }
      }
      
      #Adding the best descriptor found in the previous iteration 
      Best_descriptors <- cbind(Best_descriptors,x[,Pos_best_descriptor])
      colnames(Best_descriptors)[which(is.null(colnames(Best_descriptors))==TRUE)]<- colnames(x)[Pos_best_descriptor]#First column 
      colnames(Best_descriptors)[which(colnames(Best_descriptors)=="")] <- colnames(x)[Pos_best_descriptor]#Column consecutive to the first one: "" 
      x<-x[,-Pos_best_descriptor]#Removal of the best descriptors from the dataset
      Model_best_sum_2 <- append(Model_best_sum_2,Model_best_sum_1)#Appending best model summary of the previous while cycle
    }
    
    #Output parameters
    output<-Model
    output$descriptors<-Best_descriptors
    
    return(output)
  }
  
  ask_which_descriptor <- function(prompt = "Which descriptor/s?:")
  {
    repeat
    {
      response <- trimws(readline(prompt))
      response <- strsplit(response,",")[[1]]
      if(length(which(response %in% colnames(output$descriptors)))>0) return(response) 
      else cat("Invalid input. Please enter the proper descriptor/s name/s.\n")
    }
  }
  
  ask_yes_no <- function(prompt = "Do you want to delete one of these descriptors?[Y/N]:") 
  {
    repeat 
    {
      response <- toupper(trimws(readline(prompt)))
      if(response %in% c("Y", "N")) return(response)
      else cat("Invalid input. Please enter 'Y' or 'N'.\n")
    }
  }
  
  N<-as.numeric(readline("Set number of parameters: "))
  user_input_1<-"Y"
  while(user_input_1=="Y")
  {
    output<-best_descriptors(LT_streamflows,basin_descriptors,N)
    cat("Best descriptors:\n",paste(colnames(output$descriptors),collapse="\n"))
    user_input_1<-ask_yes_no()
    if(user_input_1=="Y") 
    {
      user_input_2<-ask_which_descriptor()
      basin_descriptors<-basin_descriptors[,-which(colnames(basin_descriptors)%in%user_input_2)]
    }  
  }
  
  R2<-round(summary(output)$adj.r.squared,4)
  dataframe<-cbind.data.frame(Code=as.character(1:length(LT_streamflows[,1])),observations=LT_streamflows[,2],estimates=exp(output$fitted.values))
  ggplot2::ggplot(dataframe, #dataframe vuoto
                  ggplot2::aes(x=observations,y=estimates,color=Code,label=Code)) +#valori: asse x -> osservazioni; asse y --> valori stimati dal modello 
    ggplot2::geom_point(size=4)+ #si rappresentano punti
    #ggplot2::geom_point(aes(x=,y=),colour="red",size=4)+ #si mette in evidenza la sezione di Marecchia SS16    
    #ggplot2::geom_point(aes(x=10.4,y=exp(Modello_4a$fitted.values)[6]),colour="orange",size=4)+ #si aggiunge per Marecchia SS16 il confronto con l'osservato da Annali Idrologici
    ggplot2::scale_color_manual(name=NULL,values=rep("black",length(LT_streamflows[,1])),labels=paste(dataframe$Code,LT_streamflows[,1],sep=" "))+
    ggplot2::guides(color=ggplot2::guide_legend(override.aes = list(shape=NA,color=NA)))+
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels=function(x) sprintf("%.0f", x),limits=c(min(dataframe[,-1]),max(dataframe[,-1])))+
    ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels =function(x) sprintf("%.0f", x),limits=c(min(dataframe[,-1]),max(dataframe[,-1])))+
    ggplot2::geom_abline(intercept = 0,slope = 1,colour="blue",lwd=1.2)+ #bisettrice primo/terzo quadrante
    ggplot2::geom_text(size=5,hjust=-0.15, vjust=-0.5)+ #etichette sezioni
    ggplot2::ggtitle(paste(N,"-parameter regression model",sep=""),subtitle = bquote("Adjusted R"^2*"="~.(R2)))+#titolo
    ggplot2::theme(axis.text=ggplot2::element_text(size=15),#dimesioni titolo e etichette assi
                  axis.title=ggplot2::element_text(size=15),
                  plot.title = ggplot2::element_text(hjust = 0.5,color = "black", size =20, face = "bold"),
                  plot.subtitle = ggplot2::element_text(hjust = 0.5,color = "black", size =15),
                  legend.key = ggplot2::element_blank(),
                  legend.background = ggplot2::element_rect(fill = "white", color = "black"))+
    #ggplot2::annotate("text", x=10.4+0.9, y=exp(Modello_4a$fitted.values)[6]+1.1, label= "8",size=5)+ #etichetta per confronto con Marecchia SS16 con Annali Idrologici
    ggplot2::labs(x=expression(paste("observations (",m^3,"/s)",sep="")),y=expression(paste("estimates (",m^3,"/s)",sep="")))#titoli assi

  # Store long term streamflow inside a list
  .internal_env$regPORFDC[["LT"]] <- exp(as.numeric(output$fitted.values[rownames(output$descriptors)%in%target_section]))
}
