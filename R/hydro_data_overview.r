#' Evaluate and Visualize Hydrometric Data Availability
#'
#' Computes stream stages, rating curves and streamflows availability matrices,
#' plotting them as heatmaps
#'
#' @param stream_stages A named list of zoo time series (one per site). Each
#'    element must be a zoo object whose time index is of class `Date` or
#'    `POSIXct`. The list names represent site names
#'
#' @param rating_curves A six column-data frame: station name (`character`),
#'   start date of applicability and end date of applicability (`Date`or`POSIXct`),
#'   minimum and maximum applicable stage (`numeric`), rating curve equation
#'  (`character`)
#'
#' @param streamflows A named list of zoo time series (one per site). Each
#'    element must be a zoo object whose time index is of class `Date` or
#'    `POSIXct`. The list names represent site names
#'
#' @param coordinates A dataframe contaning the outlet coordinates in wgs84.
#'
#' @param boundaries A feature of the catchment boundaries (one per site).
#'
#' @details Calculates annual data availability for stream stages and streamflows
#'   as proportions of expected values based on timestamp (30 minutes, 1 hour,
#'   1 day), constructs three availability matrices (stream stages, rating curves,
#'   streamflows) by stations and years and plots the study area along with
#'   stream stages and streamflow availability in grayscale color ("white": no
#'   data, "black":complete year)
#'
#' @return A List of matrices: `stream_stages` and `obs_streamflow` availability matrix,
#'   `rating_curves` binary matrix
#'
#' @note This function uses `stringr`, `zoo`, `lubridate`, `sf`,  `ggplot2` and  `ggrepel`
#'
#' @importFrom lubridate year
#' @importFrom graphics plot par axis title points legend
#'
#' @examples
#' # Assuming you have stream_stages, rating_curves, streamflows, outlet coordinates dataframes
#' and catchment boundaries loaded
#' availability_matrices<-hydro_data_overview(stream_stages,rating_curves,streamflows,
#'                                            mean_annual_flows,coordinates,boundaries)
#'
#' # Access availability matrices
#' availability_matrices$stream_stages
#' availability_matrices$rating_curves
#' availability_matrices$streamflows

get_years <- function(x)
{
  yrs <- lubridate::year(time(x))
  if (length(yrs) == 0) return(NA) else return(yrs)
}

coverage_eval<- function(x)
{
  #Most common time-steps: day, hour, 1/2 hour
  n_mins<-length(which(format(time(x),"%M")%in%c("00","30")&format(time(x),"%S")=="00"))
  n_hours<-length(which(format(time(x),"%M")=="00"&format(time(x),"%S")=="00"))
  n_days<-length(unique(as.Date(time(x))))
  if(y%%4==0) days_per_year=366 else days_per_year=365#leap or no leap year
  coverage=c(n_mins,n_hours,n_days)*c(1/c(days_per_year*48,days_per_year*24,days_per_year))
  names(coverage)<-c("half-hour","hour","day")
  return(coverage)
}

#' @export
hydro_data_overview<-function(stream_stages,rating_curves,streamflows,mean_annual_flows,coordinates,boundaries)
{
  sites_plot<-sort(unique(c(names(stream_stages),names(streamflows),mean_annual_flows[,1])))
  names(sites_plot)<-stringr::str_pad(as.character(1:length(sites_plot)),2,pad="0")

  #Stream stages availability matrix and rating curves binary matrix
  # min and max year of each series
  first_year <- min(unlist(sapply(stream_stages,get_years)),na.rm=T)
  last_year  <- max(unlist(sapply(stream_stages,get_years)),na.rm=T)
  years<-first_year:last_year
  stations<-names(stream_stages)
  stages_avail_matrix<-matrix(0,length(stations),length(years),dimnames=list(x=stations,y=years))
  rating_avail_matrix<-stages_avail_matrix

  for (s in stations)
  {
    index_s_rating<-grep(s,rating_curves[,1],ignore.case = T)
    if(length(index_s_rating)>0)
    {
      first_years<-lubridate::year(rating_curves[index_s_rating,2])
      last_years<-lubridate::year(rating_curves[index_s_rating,2])+
        round(lubridate::year(rating_curves[index_s_rating,3])+
             -lubridate::year(rating_curves[index_s_rating,2]))-1
      years_rating<-unique(unlist(Map(seq,first_years,last_years)))
      rating_avail_matrix[rownames(rating_avail_matrix)%in%s,
                      colnames(rating_avail_matrix)%in%years_rating]<-1
    }

    stages_s <- stream_stages[[s]]
    for(y in unique(lubridate::year(stages_s)))
    {
      stages_s_y <- stages_s[lubridate::year(time(stages_s))%in%y]
      stages_avail_matrix[rownames(stages_avail_matrix)%in%s,colnames(stages_avail_matrix)%in%y]<-coverage_eval(stages_s_y)["hour"]
    }
  }

  stages_df<-data.frame(name=rep(names(sites_plot)[sites_plot%in%rownames(stages_avail_matrix)],
                                 times=dim(stages_avail_matrix)[2]),
                        year=rep(colnames(stages_avail_matrix),
                                 each=dim(stages_avail_matrix)[1]),
                        value=c(stages_avail_matrix))

  rating_df<-data.frame(name=rep(names(sites_plot)[sites_plot%in%rownames(rating_avail_matrix)],
                                           times=dim(rating_avail_matrix)[2]),
                                  year=rep(colnames(rating_avail_matrix),
                                           each=dim(rating_avail_matrix)[1]),
                                  value=c(rating_avail_matrix))

  rating_df<-rating_df[rating_df[,3]==1,]
  rating_df[,3]<-as.character(rating_df[,3])

  #Greyscale plot (black=100%,white=0%)
  stages_plot<-ggplot2::ggplot(stages_df,ggplot2::aes(year,name))+
    ggplot2::geom_tile(ggplot2::aes(fill=value),
                       color = "black", linewidth = 0.1)+
    ggplot2::geom_point(data=rating_df,ggplot2::aes(x=year,y=name,colour=value),
                        shape=21,fill="yellow",size=3.5,stroke=.7) +
    ggplot2::scale_fill_gradient(name="stream stages",
                                 low = "white", high = "black",
                                 labels = scales::label_percent())+
    ggplot2::scale_color_manual(name = "rating curves",values = "black",labels="") +
    ggplot2::ggtitle("Availability matrices")+
    ggplot2::labs(x = NULL, y = NULL, fill = "Value") +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = .5),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  #Streamflows avaliability matrix
  first_year <- min(unlist(sapply(streamflows,get_years)),na.rm=T)
  last_year <- max(unlist(sapply(streamflows,get_years)),na.rm=T)
  years<-first_year:last_year
  stations<-names(streamflows)
  flows_avail_matrix<-matrix(0,length(stations),length(years),dimnames=list(x=stations,y=years))

  for (s in stations)
  {
    flows_s<-streamflows[[s]]
    for(y in unique(lubridate::year(flows_s)))
    {
      flows_s_y <- flows_s[lubridate::year(time(flows_s))%in%y]
      flows_avail_matrix[rownames(flows_avail_matrix)%in%s,colnames(flows_avail_matrix)%in%y]<-coverage_eval(flows_s_y)["day"]
    }
  }

  flows_df<-data.frame(name=rep(names(sites_plot)[sites_plot%in%rownames(flows_avail_matrix)],
                                times=dim(flows_avail_matrix)[2]),
                        year=rep(colnames(flows_avail_matrix),each=dim(flows_avail_matrix)[1]),
                        value=c(flows_avail_matrix))

  #Greyscale plot (black=100%,white=0%)
  flows_plot<-ggplot2::ggplot(flows_df,ggplot2::aes(year,name)) +
    ggplot2::geom_tile(ggplot2::aes(fill=value),
                       color = "black", linewidth = 0.1)+
    ggplot2::scale_x_discrete(labels = as.character(seq(1921,2022,2)),breaks=seq(1921,2022,2))+
    ggplot2::scale_fill_gradient(low = "white", high = "black",
                                 labels = scales::label_percent()) +
    ggplot2::labs(x = NULL, y = NULL, fill = "streamflows") +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = .5),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  #Map
  map <- rworldmap::getMap(resolution="low")
  pts_wgs84 <- sf::st_as_sf(coordinates, coords = colnames(coordinates), crs = 4326)# Convert to sf points with correct CRS
  basins_wgs84 <- sf::st_transform(boundaries, 4326)
  map <- sf::st_as_sf(map)

  map_plot<-ggplot2::ggplot()+
    ggplot2::geom_sf(data=map,fill = "grey92", color = "grey40", linewidth = 0.2)+
    ggplot2::geom_sf(data=basins_wgs84, color = "black", fill= "grey50",linewidth = .4)+
    ggplot2::geom_sf(data=pts_wgs84,color="black",size=1)+
    ggplot2::coord_sf(xlim = range(sf::st_coordinates(basins_wgs84)[,1])+c(-.5,.5),
                      ylim = range(sf::st_coordinates(basins_wgs84)[,2])+c(-.5,.5),
                      expand=FALSE)+
    ggrepel::geom_label_repel(ggplot2::aes(x=sf::st_coordinates(pts_wgs84)[,1],
                                              y=sf::st_coordinates(pts_wgs84)[,2],
                              label=names(sites_plot)[sites_plot%in%rownames(coordinates)]),
                              fill = "white",
                              segment.size=.4,
                              size = 3,
                              max.overlaps = length(sites_plot))+
  ggplot2::ggtitle("Gauged sites")+
    ggspatial::annotation_scale(
      location   = "br",
      width_hint = 0.4,        # fraction of panel width
      pad_x      = ggplot2::unit(.1, "cm"),
      pad_y      = ggplot2::unit(.1, "cm"),
      text_cex = 1.2,
      bar_cols = c("grey20", "white")
    ) +
    ggspatial::annotation_north_arrow(
      location   = "tr",
      which_north = "grid",
      pad_x      = ggplot2::unit(.1, "cm"),
      pad_y      = ggplot2::unit(.1, "cm"),
      style      = ggspatial::north_arrow_fancy_orienteering(text_size = 10),
      height = ggplot2::unit(1.7,"cm"),
      width= ggplot2::unit(1.7,"cm")
    )+
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = .5),
                   panel.grid = ggplot2::element_line(color = "grey80", linewidth = 0.3),
                   panel.background=ggplot2::element_rect(fill="white"))+
    ggplot2::labs(x = "Longitude", y = "Latitude")

  n<-length(sites_plot)
  nhalf<-ceiling(n/2)
  legend_plot<- ggplot2::ggplot()+
    ggplot2::ggtitle("Names")+
    ggplot2::annotate("text",
                      x = 0, y = 0.5,
                      label = paste(paste(names(sites_plot)[1:nhalf],sites_plot[1:nhalf]," "), collapse = "\n"),
                      size = 4.4,
                      hjust=0,
                      lineheight = 0.95)+
    ggplot2::annotate("text",
                      x = 0.5, y = 0.5,
                      label = paste(paste(names(sites_plot)[(nhalf+1):n],sites_plot[(nhalf+1):n]," "), collapse = "\n"),
                      size = 4.4,
                      hjust=0,
                      lineheight = 0.95)+
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)+
    ggplot2::theme_void()+
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = .5,size=18))

  a1 <- patchwork::area(t = 1, l = 1, b = 2, r = 2)  # p1 full height, left column
  a2 <- patchwork::area(t = 1, l = 3, b = 2, r = 6)  # p1 full height, left column
  a3 <- patchwork::area(t = 3, l = 1, b = 4, r = 2)  # p1 full height, left column
  a4 <- patchwork::area(t = 3, l = 3, b = 4, r = 6)  # p1 full height, left column

  #four panel plot: map, legend and availability matrices
  plot<-map_plot + legend_plot + stages_plot + flows_plot + patchwork::plot_layout(design=c(a1,a3,a2,a4))
  x11(height=5, width=8.4, pointsize=10)
  plot

  #Store all dataframes inside a list
  output <- list(stream_stages=stages_avail_matrix,
                 rating_curves=rating_avail_matrix,
                 streamflows=flows_avail_matrix)

  return(output)
}
