#' Time of concentration calculator
#'
#' Evaluates time of concentration for catchments, based on morphological and
#' climatic descriptors
#'
#' @param basin_descriptors A data frame with one row per catchment and
#'   multiple columns: station name in the first (`character`) and morpho‑climatic
#'   descriptors in the remaining columns (`numeric`).
#' @param method A character string indicating which time of concentration is to be
#'   computed. It has to match with one of the following: `giandotti`, `kirpich` or
#'   `pezzoli`
#'
#' @return A named numeric vector of the computed times of concentration

#' @export
time_of_concentration<-function(basin_descriptors,method)
{
  #Giandotti formulation
  if(method=="giandotti")
    tc<-(4*sqrt(basin_descriptors[,"area"])+1.5*basin_descriptors[,"LMC"])/
        (0.8*sqrt(basin_descriptors[,"elev_mean"]-basin_descriptors[,"elev_min"]))

  #Kirpich formulation
  if(method=="kirpich")
  tc<-.0195*(basin_descriptors[,"LLDP"]*1000)^.77/
      (basin_descriptors[,"LLDP_slope"]/100)^.385*1/60

  #Pezzoli formulation
  if(method=="pezzoli")
    tc<-.055*basin_descriptors[,"LLDP"]/sqrt(basin_descriptors[,"LLDP_slope"]/100)

  #aggiungere descrittori mancanti rispetto a FOCA
  tc<-round(tc,digits=2)
  names(tc)<-basin_descriptors[,1]

  return(tc)
}

