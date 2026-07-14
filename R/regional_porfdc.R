#' Regional period-of-record flow duration curve computation
#'
#' Estimates dimensionless period-of-record flow duration curve for a target site
#'
#' @param daily_streamflows A named list of zoo time series. Each element must be
#'   a zoo object whose time index is of class `Date` or `POSIXct`. The list names
#'   represent site names
#' @param basin_descriptors A data frame whose rows represents sites, columns
#'   multiple columns: station name in the first (`character`) and morpho‑climatic
#'   descriptors in the remaining columns (`numeric`).
#' @param target site name of the unagauged site (`character`)
#' @param psi exponent applied for at-site weights computation: the higher the
#'   exponent, more severe is the hydrological similiarity measure (`numeric`)
#'
#' @return A list whose first element is a data.frame of the computed flow
#'    duration curve and the second a numeric representing duration threshold
#'    used in downscaling
#'
#' @examples
#' # Assuming you have daily_streamflows, basin_descriptors, target, psi
#' regional_porfdc(daily_streamflows,basin_descriptors,target,psi)
#'
#' # This evaluates the period-of-record flow duration curve in `target`

#' @export
regional_porfdc<-function(daily_streamflows,basin_descriptors,target,psi)
{
  basin_descriptors<-basin_descriptors[basin_descriptors[,1]%in%names(daily_streamflows),]
  rownames(basin_descriptors)<-basin_descriptors[,1]
  basin_descriptors<-basin_descriptors[,-1]

  #Principal Component Analysis
  #Removing descriptors with at least one NA value
  index_c<-which(sapply(basin_descriptors,function(column) any(is.na(column))))
  if(length(index_c)>0) basin_descriptors <- basin_descriptors[,-index_c]

  pca <- stats::prcomp(basin_descriptors,scale = TRUE)#scale=TRUE: normalization (variance=1)
  summary_pca<-summary(pca)
  #Number of principal components; they have to explain at least 75% of the
  #cumulative proportion of variance
  n_pc<-min(which(summary_pca$importance[3,]>.75))
  principal_components<-pca$x[,1:n_pc]

  #principal components of the target section
  pc_target <- principal_components[rownames(principal_components)%in%target,]
  #all the other principal components
  pc_all   <- principal_components[-which(rownames(principal_components)%in%target),]

  #hydrological distance: Euclidean distance of principal components in n_pc-dimensional space
  hydro_distance <- 0
  for (i in 1:n_pc) hydro_distance <- hydro_distance + (pc_target[i] - pc_all[,i])^2
  hydro_distance<- sqrt(hydro_distance)

  #weights to be assigned for each catchment
  weights <- 1/hydro_distance^psi/sum(1/hydro_distance^psi)

  #Regional porfdcs: weighting avarage of dimensionless flow duration curves of
  # each gauged section
  #Weights must be applied to PORFDC with the same length: linear interpolation to the minimum length
  min_length<-min(lengths(daily_streamflows))
  porfdcs <- matrix(NA,dim(basin_descriptors)[1],min_length,
                    dimnames=list(x<-rownames(basin_descriptors),y=NULL))

  #Weibull exceedance probability (final p-vector in the interpolation)
  p_0<-1:(min_length)/(min_length+1)

  #Threshold duration: duration of the lowest streamflow value | Q/MAF>=1
  d_threshold<-1

  for(s in rownames(basin_descriptors))
  {
    porfdc <- porfdc_eval(daily_streamflows[[s]])
    dimless_porfdc<-porfdc/mean(porfdc)

    if(d_threshold>time(dimless_porfdc)[which.min(dimless_porfdc>=1)])
      d_threshold <- time(dimless_porfdc)[which.min(dimless_porfdc>=1)]

    #Resampling of the dimensionless PORFDC, so that all series have the same length
    porfdcs[rownames(porfdcs)%in%s,]<-approx(x=time(dimless_porfdc),y=dimless_porfdc,
                                             xout=p_0,method="linear")$y
  }

  #RoI approach for weighting the dimensionless flow duration curves in order to
  #derive the dimensionless curve of the target catchment
  porfdcs <- porfdcs[-which(rownames(porfdcs)==target),]
  regional_porfdc <- as.numeric(t(porfdcs)%*%weights)

  regional_porfdc<-unique(regional_porfdc) #remove streamflows redundancies
  regional_porfdc<-data.frame(p=1:length(regional_porfdc)/(length(regional_porfdc)+1),
             Q=regional_porfdc,row.names = NULL)

  #Extension of regional porfdc in order to be applicable for every durations of the
  #donor site
  regional_porfdc<- rbind(c(0,max(porfdcs[,1])),#maximum regional streamflow
                                                #associated to duration=0
                          regional_porfdc,
                          c(1,min(regional_porfdc$Q)))#minimum streamflow of the target
                                                      #associated to duration=1

  # Store dimensionless period-of-record flow duration curve for the target site
  #and duration
  output<-list(regional_porfdc=regional_porfdc,d_threshold=d_threshold)
  return(output)
}
