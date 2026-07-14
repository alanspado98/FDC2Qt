#' Dataset: North-Eastern Italy stream stages
#'
#' Hourly and semi-hourly stream stages of some gauged cross-sections in
#' Emilia-Romagna region from 1986 to 2021; values are expressed in meters
#' above a zero level (specific to each cross-section)
#'
#' @format A named list of 28 time series, where the names correspond to cross-
#' section names. Each element is a \code{zoo} object with:
#' \describe{
#'   \item{index}{Date and time (POSIXct)}
#'   \item{value}{Stream stages (numeric)}
#' }
#' @source https://simc.arpae.it/dext3r/
"stream_stages"

#' Dataset: North-Eastern Italy rating-curves
#'
#' Stream stage to streamflow equations of gauged cross-sections in Emilia-
#' Romagna region from 1986 to 2021; hydrometric levels are expressed in meters,
#' streamflows in cubic meters per second
#'
#' @format A data frame with 2,420 rows and 6 columns:
#' \describe{
#'   \item{Name}{Cross-section name (character)}
#'   \item{StartDate}{Starting validity period (POSIXct)}
#'   \item{EndDate}{Ending validity period (POSIXct)}
#'   \item{LowerLevel}{Lower validity level (numeric)}
#'   \item{UpperLevel}{Upper validity level (numeric)}
#'   \item{UserDefinedEquation}{Equation definition (character)}
#' }
#' @source https://www.arpae.it/it/temi-ambientali/meteo/report-meteo/annali-idrologici
"rating_curves"

#' Dataset: North-Eastern Italy streamflow data
#'
#' Daily streamflow of some Emilia-Romagna region gauged cross-sections
#' from 1921 to 2022; values are expressed in cubic meters per second
#'
#' @format A named list of 8 time series, where the names correspond to cross-
#' section names. Each element is a \code{zoo} object with:
#' \describe{
#'   \item{index}{Date and time (POSIXct)}
#'   \item{value}{Streamflows (numeric)}
#' }
#' @source https://simc.arpae.it/dext3r/
"streamflows"

#' Dataset: North-Eastern Italy mean annual precipitation data
#'
#' Mean annual precipitation of some Emilia-Romagna region gauged basins
#' from 2000 to 2021; values are expressed in millimetrs
#'
#' @format A named list of 15 time series, where the names correspond to cross-
#' section names. Each element is a \code{zoo} object with:
#' \describe{
#'   \item{index}{Year (numeric)}
#'   \item{value}{Total annual precipitation (numeric)}
#'#' }
#' @source https://simc.arpae.it/dext3r/
"annual_precipitation"

#' @format
#' A data frame with 11 rows and 2 columns:
#' \describe{
#'   \item{Name}{Basin cross-section name (character)}
#'   \item{Value}{Long-term average daily streamflow (numeric)}
#'#' }
#' @source https://www.arpae.it/it/temi-ambientali/meteo/report-meteo/annali-idrologici
"mean_annual_flows"

#' Dataset: North-Eastern Italy Basin Descriptors
#'
#' Morphological and climatic features of some Emilia-Romagna region gauged
#' and ungauged basins
#'
#' @format A data frame with multiple rows and 68 variables:
#' \describe{
#'   \item{name}{Corss-section name (character,"River@@Location" format)}
#'   \item{area}{Basin area in square kilometers (numeric)}
#'   \item{LMC}{Length of the main channel in kilometers (numeric)}
#'   \item{elev_mean}{Mean elevation of the basin in meters above sea level (numeric)}
#'   \item{elev_closure}{elevattion of the closure section in meters above sea level (numeric)}
#'   \item{x_g}{Easting coordinate of the basin centroid in meters (numeric)}
#'   \item{y_g}{Northing coordinate of the basin centroid in meters (numeric)}
#'   \item{LLDP}{Length of the longest drainage path in kilometers (numeric)}
#'   \item{drain_dens}{Drainage density in inverse kilometers (numeric)}
#'   \item{elev_max}{Maximum elevation in meters above sea level (numeric)}
#'   \item{elev_min}{Minimum elevation in meters above sea level (numeric)}
#'   \item{elev_2.5}{Elevation at 2.5% of the hypsographic curve (numeric)}
#'   \item{elev_5}{Elevation at 5% of the hypsographic curve (numeric)}
#'   \item{elev_10}{Elevation at 10% of the hypsographic curve (numeric)}
#'   \item{elev_25}{Elevation at 25% of the hypsographic curve (numeric)}
#'   \item{elev_50}{Elevation at 50% of the hypsographic curve (numeric)}
#'   \item{elev_75}{Elevation at 75% of the hypsographic curve (numeric)}
#'   \item{elev_90}{Elevation at 90% of the hypsographic curve (numeric)}
#'   \item{elev_95}{Elevation at 95% of the hypsographic curve (numeric)}
#'   \item{elev_97.5}{Elevation at 97.5% of the hypsographic curve (numeric)}
#'   \item{LLDP_slope}{Mean slope along the longest drainage path in percentage (numeric)}
#'   \item{dir_length}{Vector length from centroid to outlet in kilometers (numeric)}
#'   \item{orient}{Basin orientation in degrees (numeric)}
#'   \item{slope1}{Mean slope of type 1 in percentage (numeric)}
#'   \item{slope2}{Mean slope of type 2 in percentage (numeric)}
#'   \item{aspect}{Mean aspect in degrees (numeric)}
#'   \item{elong_r}{Elongation Ratio (numeric)}
#'   \item{shape_f}{Shape factor (numeric)}
#'   \item{width_mean}{Mean width function value (numeric)}
#'   \item{width_var}{Variance of the width function (numeric)}
#'   \item{width_skw{Skewness of the width function (numeric)}
#'   \item{width_kur}{Kurtosis of the width function (numeric)}
#'   \item{MHL}{Mean hillslope length in kilometers (numeric)}
#'   \item{topo_d}{Topological diameter in kilometers (numeric)}
#'   \item{HS_num_1}{Number of first-order streams (integer)}
#'   \item{HS_num_2}{Number of second-order streams (integer)}
#'   \item{HS_num_3}{Number of third-order streams (integer)}
#'   \item{R_b}{Bifurcation ratio of stream orders (numeric)}
#'   \item{R_l}{Stream length ratio (numeric)}
#'   \item{R_a}{Drainage area ratio (numeric)}
#'   \item{R_s}{Stream slope ratio (numeric)}
#'   \item{width_5}{Width function at 5% in meters (numeric)}
#'   \item{width_10}{Width function at 10% in meters (numeric)}
#'   \item{width_15}{Width function at 15% in meters (numeric)}
#'   \item{width_30}{Width function at 30% in meters (numeric)}
#'   \item{width_40}{Width function at 40% in meters (numeric)}
#'   \item{width_50}{Width function at 50% in meters (numeric)}
#'   \item{width_60}{Width function at 60% in meters (numeric)}
#'   \item{width70}{Width function at 70% in meters (numeric)}
#'   \item{width_85}{Width function at 85% in meters (numeric)}
#'   \item{width_95}{Width function at 95% in meters (numeric)}
#'   \item{HS_length_1}{Mean length of first-order streams in kilometers (numeric)}
#'   \item{HS_length_2}{Mean length of second-order streams in kilometers (numeric)}
#'   \item{HS_length_3}{Mean length of third-order streams in kilometers (numeric)}
#'   \item{HS_area_1}{Drainage area of first-order streams in square kilometers (numeric)}
#'   \item{HS_area_2}{Drainage area of second-order streams in square kilometers (numeric)}
#'   \item{HS_area_3}{Drainage area of third-order streams in square kilometers (numeric)}
#'   \item{HS_slope_1}{Slope of first-order streams in percentage (numeric)}
#'   \item{HS_slope_2}{Slope of second-order streams in percentage (numeric)}
#'   \item{HS_slope_3}{Slope of third-order streams in percentage (numeric)}
#'   \item{TSL}{Total stream length in kilometers (numeric)}
#'   \item{TI_mean}{Mean of the topographic index (numeric)}
#'   \item{TI_sd}{Standard deviation of the topographic index (numeric)}
#'   \item{MAP}{Mean annual precipitation over the last 20 years in millimeters (numeric)}
#'   \item{MAP_sd}{Standard deviation of annual precipitation over the last 20 years in millimeters (numeric)}
#'   \item{MAPE}{Mean annual potential evapotranspiration over the last 20 years in millimeters (numeric)}
#'   \item{MAPE_sd}{Standard deviation of annual potential evapotranspiration over the last 20 years in millimeters (numeric)}
#'   \item{MAEE}{Mean annual actual evapotranspiration over the last 20 years in millimeters (numeric)}
#'   \item{MAEE_sd}{Standard deviation of annual actual evapotranspiration over the last 20 years in millimeters (numeric)}
#' }
#' @source #Morphological descriptors
#'         https://www.academia.edu/5241352/Atlante_bacini_piemontesi_LR
#'         https://sdi.eea.europa.eu/catalogue/copernicus/api/records/66fa7dca-8772-4a5d-9d56-2caba4ecd36a
#'         #Climatic descriptors
#'         #https://www.isprambiente.gov.it/pre_meteo/idro/BIGBANG_ISPRA.html
"basin_descriptors"
