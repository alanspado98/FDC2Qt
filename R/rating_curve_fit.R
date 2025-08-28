#' Rating Curves Fit for Unrated Gauging Stations
#'
#' Fit a power-law rating curves to stations that have contemporaneity of stream stages and streamflow data
#' but lack rating curves
#'
#' @param streamflows a `data.frame` a data frame with three columns: station name
#'   (character), time istant of observation (POSIXct), streamflow observation (numeric)
#'
#' @param stream_stages a data frame with three columns: station name
#'   (character), time istant of observation (POSIXct), stream stage
#'   observation (numeric)
#'
#' @return store output in the internal environment `.internal_env`:
#'  \item \code{.internal_env$input_data}: List of data frames: `rating_curves` by running 
#'
#' @details 
#' Perform the calibration using a smoothed spline and piecewise optimization
#' to identify breakpoints and fit segmented rating curve equations. The resulting rating curves are appended
#' to the `rating_curves` dataset and passed to the `input_data()` function for updating internal data.
#' The function:
#' \itemize{
#'   \item Identifies stations that have both streamflows and hydrometric levels but no associated rating curves.
#'   \item Aggregates hourly hydrometric levels to daily averages.
#'   \item Matches dates with both hydrometric level and streamflow observations.
#'   \item Fits a smoothed spline to the matched dataset, and detects breakpoints using a median filter on the spline’s slope.
#'   \item Fits a power-law model \eqn{Q = a + b(H + \delta)^c} to each segment using non-linear optimization (`optim`).
#'   \item Adds newly calibrated rating curves to the existing `rating_curves` dataset with fields: 
#'         \code{Name}, \code{StartDate}, \code{EndDate}, \code{LowerLevel}, \code{UpperLevel}, and \code{UserDefinedEquation}.
#' }
#'
#' The fitted equations are written in the format:  
#' \deqn{Q = a + b(H + \delta)^c}
#' where \eqn{\delta} is a shift applied if the minimum hydrometric level is negative.
#'
#' @note 
#' - Internal data is accessed using `use_internal_data(list_name = "DA_matrices", ...)`.
#' - The function assumes existence of `input_data()` to register updated rating curves.
#' - Breakpoints are identified from changes in the slope of a smoothed version of the rating curve.
#' - Only stations with valid data and no existing rating curve are processed.
#'
#' @importFrom lubridate year
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' @importFrom stats smooth.spline optim
#' @importFrom caTools runquantile
#'
#' @examples
#' \dontrun{
#' fit_rating_curve(streamflows, stream_stages)
#' }
#'
#' @export


fit_rating_curve<-function(streamflows,stream_stages)
{
  hl_DA_matrix<-use_internal_data(list_name = "DA_matrices",dataset_name = "stream_stages")
  rc_DA_matrix<-use_internal_data(list_name = "DA_matrices",dataset_name = "rating_curves")
  str_DA_matrix<-use_internal_data(list_name = "DA_matrices",dataset_name = "obs_streamflows")
  
  #Section with no rating curve but with hydrometric levels and streamflows
  stations<-rownames(str_DA_matrix)[rownames(str_DA_matrix)%in%rownames(hl_DA_matrix)][rownames(str_DA_matrix)[rownames(str_DA_matrix)%in%rownames(hl_DA_matrix)]%in%rownames(rc_DA_matrix)[which(apply(rc_DA_matrix,1,sum)==0)]]
  `%>%` <- magrittr::`%>%`
  for(r in stations)
  {
    #Same station 
    hl_r <-stream_stages[stream_stages[,1]%in%r,]
    hl_r<- hl_r %>%#Daily hydrometric levels
      dplyr::group_by(Date = substring(Datetime,1,10)) %>%
      dplyr::summarise(Value = mean(Value, na.rm = TRUE)) %>%
      as.data.frame() 
    str_r<-obs_streamflows[obs_streamflows[,1]%in%r,-1]
    
    #Remove NA values
    hl_r<-hl_r[!is.na(hl_r[,2]),]
    str_r<-str_r[!is.na(str_r[,2]),]
      
    #years with hydrometric levels but no streamflows
    hl_index<-which(lubridate::year(hl_r[,1])%in%names(which(str_DA_matrix[rownames(str_DA_matrix)%in%r,]==0)))#years with no streamflows
    
    if(length(hl_index)>0)
    {
      #Same days
      hl_r_c<-hl_r[which(hl_r[,1]%in%str_r[,1]),]
      str_r_c<-str_r[str_r[,1]%in%hl_r_c[,1],]
      
      rc_data<-data.frame(hl=hl_r_c[,2],str=str_r_c[,2])
      rc_data<-rc_data[order(rc_data$hl),]
      if(min(rc_data$hl)<0) hl_shift<-abs(min(rc_data$hl)) else hl_shift<-0
      rc_data$hl<-rc_data$hl+hl_shift
        
      #Find breakpoints
      #Apply smoothing (adjust spar for smoothing level)
      smooth_fit <- smooth.spline(rc_data,spar=.8)#3000  # Adjust spar for more/less smoothnes
      index_min<-max(which(smooth_fit$x%in%min(smooth_fit$x[smooth_fit$x>0]))[1],which(smooth_fit$y%in%min(smooth_fit$y[smooth_fit$y>0]))[1])
      index_max<-min(which(smooth_fit$x%in%max(smooth_fit$x[smooth_fit$x>0]))[1],which(smooth_fit$y%in%max(smooth_fit$y[smooth_fit$y>0]))[1])

      breakpoints_x<-c(smooth_fit$x[index_min],smooth_fit$x[index_max])
      breakpoints_y<-c(smooth_fit$y[index_min],smooth_fit$y[index_max])
      
      #Median filter to uniform the smoothing
      #Major breakpoints 
      breakpoint_indexes<-c()
      smooth_median<-caTools::runquantile(smooth_fit$y,length(rc_data[,1])/3,probs = 0.5,type=7,endrule="NA")
      dy_dx<-diff(smooth_median)/smooth_fit$x[-length(smooth_fit$x)]
      relative_error <- abs(diff(dy_dx)/dy_dx[-length(dy_dx)])
      breakpoint_indexes<-append(breakpoint_indexes,which(relative_error>50&relative_error<Inf)+2)
  
      options(warn = -1)  # Disable warnings
      breakpoints_x<-c(breakpoints_x,smooth_fit$x[breakpoint_indexes])
      breakpoints_y<-c(breakpoints_y,smooth_fit$y[breakpoint_indexes])
      breakpoints_x<-sort(breakpoints_x)
      breakpoints_y<-sort(breakpoints_y)
      
      power_model <- function(c) 
      {
        b <- (breakpoints_y[i_bp+1] - breakpoints_y[i_bp]) / 
          (breakpoints_x[i_bp+1]^c - breakpoints_x[i_bp]^c)
        a <- breakpoints_y[i_bp+1] - b * breakpoints_x[i_bp+1]^c
        Est <- a + b * x^c
        sumsq <- sum((y - Est)^2, na.rm = TRUE)
        return(sumsq)  
      }
      
      for(i_bp in 1:(length(breakpoints_x)-1))
      {
        lowerlevel<-breakpoints_x[i_bp]
        upperlevel<-breakpoints_x[i_bp+1]
        
        if(i_bp==1) lowerlevel<-rc_data[1,1]
        if(i_bp==length(breakpoints_x)-1) upperlevel<-rc_data[length(rc_data[,1]),1]
        
        x<-rc_data$hl[rc_data$hl>=lowerlevel & rc_data$hl<=upperlevel]
        y<-rc_data$str[rc_data$hl>=lowerlevel & rc_data$hl<=upperlevel]
        
        pow_model<-optim(par=1,fn=power_model,method="Brent",lower=-10,upper=10)
        
        c <- pow_model$par
        b <- (breakpoints_y[i_bp+1]-breakpoints_y[i_bp])/(breakpoints_x[i_bp+1]^c-breakpoints_x[i_bp]^c)
        a <- breakpoints_y[i_bp]-b*breakpoints_x[i_bp]^c
        
        #Add the new rating curve
        Name         <- paste(gsub("([a-z])([A-Z])","\\1 \\2",strsplit(r,"@")[[1]]),collapse = " a ")
        StartDate    <- paste0(hl_r[hl_index,1][1], " 00:00:00")#First day of the year in which there are hydrometric levels but neither streamflows nor rating curves
        EndDate      <- paste0(hl_r[hl_index,1][length(hl_index)], " 23:00:00")#Last day of the year in which there are hydrometric levels but neither streamflows nor rating curves
        LowerLevel <- lowerlevel-hl_shift
        UpperLevel <- upperlevel-hl_shift
        UserDefinedEquation <- paste0(round(b,3),"*(X+",hl_shift,")^",round(c,3),if(a>=0) "+" else "",round(a,3))#Controllo sul segno
        rating_curves<-rbind.data.frame(rating_curves,
                                        cbind.data.frame(Name,StartDate,EndDate,LowerLevel,UpperLevel,UserDefinedEquation))
      }
      
    }
    options(warn = 0)   # Re-enable warnings
  }
  input_data(rating_curves)#It should be revaluated only for rating curves
}