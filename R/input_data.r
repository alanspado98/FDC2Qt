#' Evaluate and Visualize Hydrometric Data Availability
#'
#' Computes stream stages, rating curves and streamflows availability matrices
#' and plot them as heatmaps
#'
#' @param stream_stages a data frame with three columns: station name
#'   (character), time istant of observation (character), stream stage
#'   observation (numeric)
#'
#' @param rating_curves a data frame with six columns: station name
#'   (character), start date of applicability (character), end date of
#'   applicability (character), minimum applicable stage (numeric), maximum
#'   applicable stage (numeric), rating curve equation (character)
#'
#' @param obs_streamflows a data frame with three columns: station name
#'   (stations in common with stream_stages must match), time istant of
#'   observation (character), streamflow observation (numeric)
#'
#' @details Parse a date-time or date string in common international format,
#' match rating curve to stream stages by extracting `location` from
#' `station@location`, calculates annual data availability for stream stages and
#' streamflows as proportions of expected values based on timestamp (5 min, 30
#' min, 1 hour), constructs three availability matrices (stream stages, rating
#' curves, streamflows) by stations and years and plots stream stages and daily
#' streamflow availability in grayscale color ("white": no data, "black":
#' complete year)
#'
#' @return store output in the internal environment `.internal_env`:
#' \itemize{
#'   \item \code{.internal_env$input_data}: List of data frames: `stream_stages`,`rating_curves`, `obs_streamflows`
#'   \item \code{.internal_env$DA_matrices}: List of matrices: `stream_stages` and `obs_streamflow` availability matrix,
#'   `rating_curves` binary matrix
#' }
#'
#' @note This function uses `lubridate` and `plot.matrix` assumes a global
#'   environment variable `.internal_env` exists
#'
#' @importFrom lubridate year
#' @importFrom graphics plot par axis title points legend
#'
#' @examples
#' # Assuming you have stream_stages, rating_curves and obs_streamflows dataframes loaded
#' input_data(stream_stages, rating_curves, obs_streamflows)
#'
#' # Access availability matrices
#' .internal_env$DA_matrices$stream_stages
#' .internal_env$DA_matrices$rating_curves
#' .internal_env$DA_matrices$obs_streamflows
#' @export
load("data/obs_streamflows.rda")
load("data/stream_stages.rda")
load("data/rating_curves.rda")

input_data<-function(stream_stages,rating_curves,obs_streamflows)
{
  #obs_streamflows<-use_internal_data("streamflows","GG")
  stream_stages[,2] <- std_datetime(stream_stages[,2])
  
  rating_curves[,2]<-std_datetime(rating_curves[,2])
  rating_curves[,3]<-std_datetime(rating_curves[,3])
  
  obs_streamflows[,2]<-std_datetime(obs_streamflows[,2])
  
  #stream stages availability matrix and binary matrix
  years<-min(lubridate::year(stream_stages[,2]),na.rm=T):max(lubridate::year(stream_stages[,2]),na.rm=T)
  stations<-sort(unique(stream_stages[,1]))
  ss_DA_matrix<-matrix(0,length(stations),length(years),dimnames=list(x=stations,y=years))
  rc_DA_matrix<-ss_DA_matrix
  
  for (r in stations)
  {
    r_as_rc<-gsub("([a-z])([A-Z])","\\1 \\2",strsplit(r,"@")[[1]][2])
    Index_r_rc<-grep(r_as_rc,rating_curves[,1],ignore.case = T)
    if(length(Index_r_rc)>0)
    {   
      first_year_rc<-lubridate::year(rating_curves[Index_r_rc,2])
      last_year_rc<-lubridate::year(rating_curves[Index_r_rc,2])+round(as.numeric(rating_curves[Index_r_rc,3]-rating_curves[Index_r_rc,2])/365)-1
      years_rc<-unique(unlist(Map(seq,first_year_rc,last_year_rc)))
      rc_DA_matrix[rownames(rc_DA_matrix)%in%r,colnames(rc_DA_matrix)%in%years_rc]<-1
    }
  
    ss_r<-stream_stages[stream_stages[,1]%in%r,]
    
    for(c in unique(lubridate::year(ss_r[,2])))
    {
      ss_r_c<-ss_r[lubridate::year(ss_r[,2])%in%c,]
      ss_ts<-diff.POSIXt(ss_r_c[,2])
      n_ss_ts<-c(length(which(ss_ts==1)),length(which(ss_ts==30)),length(which(ss_ts==5)))
      n_ss_ts[n_ss_ts!=0]<-n_ss_ts[n_ss_ts!=0]+1
      
      if(c%%4==0) days_per_year=366 else days_per_year=365#Check for leap years
      N=(n_ss_ts/length(ss_r_c[,1]))%*%c(days_per_year*24,days_per_year*24*2,days_per_year*24*12)
      
      ss_DA_matrix[rownames(ss_DA_matrix)%in%r,colnames(ss_DA_matrix)%in%c]<-length(ss_r_c[which(!is.na(ss_r_c[,3])),3])/N
    }
  }
  
  library(plot.matrix)
  #greyscale plot (black=100%,white=0%)
  par(mar=c(4,14,6,8),mgp = c(2,.7,0),xpd=F)
  plot(ss_DA_matrix,
       symbols(1,1,squares=5),
       col=paste("gray",100:0,sep=""),
       axis.col = NULL,
       cex.axis=.9,
       axis.row = list(side=2,las=1),
       key=NULL,
       main="",
       xlab="",ylab="")
  axis(side=1,las=2,at=1:length(ss_DA_matrix[1,]),labels=colnames(ss_DA_matrix),cex.axis=1)
  axis(side=3,las=2,at=1:length(ss_DA_matrix[1,]),labels=colnames(ss_DA_matrix),cex.axis=1)

  title("Stream stage and rating curve",line=1,cex=2)
  points(col(rc_DA_matrix)[rc_DA_matrix==1],nrow(rc_DA_matrix)-row(rc_DA_matrix)[rc_DA_matrix==1]+1,pch=21,col="black",bg="yellow",cex=1)  # Adjust size/color if needed
  par(xpd=T)
  legend(x=36,y=25,title="Stream stage\navailability (-)",legend=c("100%","75%","50%","25%","0%"),fill=c("grey0","grey25","grey50","grey75","grey100"),ncol=1,cex=1)
  legend(x=36,y=12,legend="Rating curve\n availability",pch=21,col="black",pt.bg="yellow",cex=1)
  
  #streamflows avaliability matrix
  years<-min(lubridate::year(obs_streamflows[,2]),na.rm=T):max(lubridate::year(obs_streamflows[,2]),na.rm=T)
  stations<-sort(unique(obs_streamflows[,1][!is.na(obs_streamflows[,2])]))
  str_DA_matrix<-matrix(0,length(stations),length(years),dimnames=list(x=stations,y=years))

  for (r in stations)
  {
    str_r<-obs_streamflows[obs_streamflows[,1]%in%r,]
    for(c in unique(lubridate::year(str_r[,2])))
    {
      if(c%%4==0) N=366 else N=365#Check for leap years
      str_r_c<-str_r[lubridate::year(str_r[,2])%in%c,]
      str_DA_matrix[rownames(str_DA_matrix)%in%r,colnames(str_DA_matrix)%in%c]<-length(which(!is.na(str_r_c[,3])))/N
    }
  }
library(plot.matrix)
  str_DA_matrix<-str_DA_matrix[-c(4,12),]
  #greyscale plot (black=100%,white=0%)
  par(mar=c(4,14,3,8),mgp = c(2,.7,0),xpd=F)
  plot(str_DA_matrix,
       col=paste("gray",100:0,sep=""),
       axis.col = list(side=1,las=2),cex.axis=1,
       axis.row = list(side=2,las=1),
       key=NULL,
       main="",
       xlab="",ylab="")
  axis(side=3,las=2,at=1:length(str_DA_matrix[1,]),labels=colnames(str_DA_matrix),cex.axis=1)
  title("Daily streamflows",line=1,cex=2)
  
  par(xpd=T)
  legend(x=105,y=7,title="Streamflow\navailability (-)",legend=c("100%","75%","50%","25%","0%"),fill=c("grey0","grey25","grey50","grey75","grey100"),ncol=1,cex=1)
  
  # Store all dataframes inside a list
  .internal_env$input_data  <- list(stream_stages=stream_stages,rating_curves=rating_curves,obs_streamflows=obs_streamflows)
  .internal_env$DA_matrices <- list(stream_stages=ss_DA_matrix,rating_curves=rc_DA_matrix,obs_streamflows=str_DA_matrix)
}
