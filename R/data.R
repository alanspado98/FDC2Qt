#' Dataset: North-Eastern Italy hydrometric-level data 
#'
#' Hourly and semi-hourly hydrometric levels of some Emilia-Romagna region
#' gauged cross-sections from 1986 to 2021; values are expressed in meters 
#' above a zero level specific to each cross-section 
#'
#' @format ## `hydrometric_levels`
#' A data frame with 11,608,248 rows and 3 columns:
#' \describe{
#'   \item{Name}{Cross-section name (character, "River@@Location" format)}
#'   \item{Datetime}{Measurement date and time (character, "%Y-%m-%d %H:%M" format)}
#'   \item{Value_m}{Hydrometric-level values (numeric)}
#'#' }
#' @source https://simc.arpae.it/dext3r/
"hydrometric_levels"

#' Dataset: North-Eastern Italy rating-curves
#'
#' Hydrometric-level to streamflow equations of Emilia-Romagna region
#' gauged cross-sections from 1986 to 2021; hydrometric levels are expressed 
#' in meters, streamflows in cubic meters per second 
#' 
#' @format ## `rating_curves`
#' A data frame with 2,420 rows and 6 columns:
#' \describe{
#'   \item{Name}{Cross-section name (character,"Location" format)}
#'   \item{StartDate}{Starting validity period (character, "%Y-%m-%d %H:%M:%S" format)}
#'   \item{EndDate}{Ending validity period (character, "%Y-%m-%d %H:%M:%S" format)}
#'   \item{LowerLevel_m}{Lower validity level (numeric)}
#'   \item{UpperLevel_m}{Upper validity level (numeric)}
#'   \item{UserDefinedEquation}{Equation definition (character)}
#' }
#' @source https://www.arpae.it/it/temi-ambientali/meteo/report-meteo/annali-idrologici
"rating_curves"

#' Dataset: North-Eastern Italy streamflow data
#' 
#' Daily streamflow of some Emilia-Romagna region gauged cross-sections
#' from 1921 to 2022; values are expressed in cubic
#' meters per second
#'
#' @format ## `streamflows`
#' A data frame with 107,321 rows and 3 columns:
#' \describe{
#'   \item{Name}{Cross-section name (character, "River@@Location" format)}
#'   \item{Datetime}{Measurement date and time (character, "%Y-%m-%d" format)}
#'   \item{Values_m3/s}{Streamflow values (numeric)}
#' }
#' @source https://simc.arpae.it/dext3r/
"streamflow"

#' Dataset: North-Eastern Italy mean annual precipitation data 
#'
#' Mean annual precipitation of some Emilia-Romagna region gauged basins
#' from 2000 to 2021; values are expressed in millimetrs
#'
#' @format ## `mean_annual_precipitation`
#' A data frame with 11,608,248 rows and 3 columns:
#' \describe{
#'   \item{Name}{Basin cross-section name (character, "River@@Location" format)}
#'   \item{2001}{Mean annual values of each basins in 2001 (numeric)}
#'   \item{2002}{Mean annual values of each basins in 2002 (numeric)}
#'   ...
#'#' }
#' @source https://simc.arpae.it/dext3r/
"hydrometric_levels"

#' Dataset: North-Eastern Italy Basin Descriptors
#'
#' Morphological and climatic features of some Emilia-Romagna region gauged
#' and ungauged basins
#'
#' @format A data frame with multiple rows and 68 variables:
#' \describe{
#'   \item{Name}{Corss-section name (character,"River@@Location" format)}
#'   \item{Area_km2}{Basin area in square kilometers (numeric)}
#'   \item{CentroidEst_m}{Easting coordinate of the basin centroid in meters (numeric)}
#'   \item{CentroidNorth_m}{Northing coordinate of the basin centroid in meters (numeric)}
#'   \item{MeanElevation_masl}{Mean elevation of the basin in meters above sea level (numeric)}
#'   \item{LongestDrainagePathLength_km}{Length of the longest drainage path in kilometers (numeric)}
#'   \item{DrainageDensity_1/km}{Drainage density in inverse kilometers (numeric)}
#'   \item{MaxElevation_masl}{Maximum elevation in meters above sea level (numeric)}
#'   \item{MinElevation_masl}{Minimum elevation in meters above sea level (numeric)}
#'   \item{HypsographicCurve2.5%_masl}{Elevation at 2.5% of the hypsographic curve (numeric)}
#'   \item{HypsographicCurve5%_masl}{Elevation at 5% of the hypsographic curve (numeric)}
#'   \item{HypsographicCurve10%_masl}{Elevation at 10% of the hypsographic curve (numeric)}
#'   \item{HypsographicCurve25%_masl}{Elevation at 25% of the hypsographic curve (numeric)}
#'   \item{HypsographicCurve50%_masl}{Elevation at 50% of the hypsographic curve (numeric)}
#'   \item{HypsographicCurve75%_masl}{Elevation at 75% of the hypsographic curve (numeric)}
#'   \item{HypsographicCurve90%_masl}{Elevation at 90% of the hypsographic curve (numeric)}
#'   \item{HypsographicCurve95%_masl}{Elevation at 95% of the hypsographic curve (numeric)}
#'   \item{HypsographicCurve97.5%_masl}{Elevation at 97.5% of the hypsographic curve (numeric)}
#'   \item{MainChannelLength_km}{Length of the main channel in kilometers (numeric)}
#'   \item{MeanSlopeLDP_%}{Mean slope along the longest drainage path in percentage (numeric)}
#'   \item{CentroidOutletVectorLength_km}{Vector length from centroid to outlet in kilometers (numeric)}
#'   \item{Orientation_degr}{Basin orientation in degrees (numeric)}
#'   \item{MeanSlope1_%}{Mean slope of type 1 in percentage (numeric)}
#'   \item{MeanSlope2_%}{Mean slope of type 2 in percentage (numeric)}
#'   \item{MeanAspect_degr}{Mean aspect in degrees (numeric)}
#'   \item{ElongationRatio}{Elongation Ratio (numeric)}
#'   \item{ShapeFactor}{Shape factor (numeric)}
#'   \item{MeanWidthFunction}{Mean width function value (numeric)}
#'   \item{VarWidthFunction}{Variance of the width function (numeric)}
#'   \item{SkewWidthFunction}{Skewness of the width function (numeric)}
#'   \item{KurtoWidthFunction}{Kurtosis of the width function (numeric)}
#'   \item{MeanHillslopeLength_km}{Mean hillslope length in kilometers (numeric)}
#'   \item{TopologicalDiameter_km}{Topological diameter in kilometers (numeric)}
#'   \item{1stOrderStreams}{Number of first-order streams (integer)}
#'   \item{2ndOrderStreams}{Number of second-order streams (integer)}
#'   \item{3rdOrderStreams}{Number of third-order streams (integer)}
#'   \item{BifurcationRatio}{Bifurcation ratio of stream orders (numeric)}
#'   \item{StreamLengthRatio}{Stream length ratio (numeric)}
#'   \item{DrainageAreaRatio}{Drainage area ratio (numeric)}
#'   \item{StreamSlopeRatio}{Stream slope ratio (numeric)}
#'   \item{WidthFunction5%_m}{Width function at 5% in meters (numeric)}
#'   \item{WidthFunction10%_m}{Width function at 10% in meters (numeric)}
#'   \item{WidthFunction15%_m}{Width function at 15% in meters (numeric)}
#'   \item{WidthFunction30%_m}{Width function at 30% in meters (numeric)}
#'   \item{WidthFunction40%_m}{Width function at 40% in meters (numeric)}
#'   \item{WidthFunction50%_m}{Width function at 50% in meters (numeric)}
#'   \item{WidthFunction60%_m}{Width function at 60% in meters (numeric)}
#'   \item{WidthFunction70%_m}{Width function at 70% in meters (numeric)}
#'   \item{WidthFunction85%_m}{Width function at 85% in meters (numeric)}
#'   \item{WidthFunction95%_m}{Width function at 95% in meters (numeric)}
#'   \item{1stOrderStreamLengths_km}{Mean length of first-order streams in kilometers (numeric)}
#'   \item{2ndOrderStreamLengths_km}{Mean length of second-order streams in kilometers (numeric)}
#'   \item{3rdOrderStreamLengths_km}{Mean length of third-order streams in kilometers (numeric)}
#'   \item{1stOrderDrainageAreas_km2}{Drainage area of first-order streams in square kilometers (numeric)}
#'   \item{2ndOrderDrainageAreas_km2}{Drainage area of second-order streams in square kilometers (numeric)}
#'   \item{3rdOrderDrainageAreas_km2}{Drainage area of third-order streams in square kilometers (numeric)}
#'   \item{1stOrderStreamSlopes_%}{Slope of first-order streams in percentage (numeric)}
#'   \item{2ndOrderStreamSlopes_%}{Slope of second-order streams in percentage (numeric)}
#'   \item{3rdOrderStreamSlopes_%}{Slope of third-order streams in percentage (numeric)}
#'   \item{TotalStreamLength_km}{Total stream length in kilometers (numeric)}
#'   \item{MeanTopographicIndex}{Mean of the topographic index (numeric)}
#'   \item{StDevTopographicIndex}{Standard deviation of the topographic index (numeric)}
#'   \item{MeanMAPLast20_mm}{Mean annual precipitation over the last 20 years in millimeters (numeric)}
#'   \item{StDevMAPLast20_mm}{Standard deviation of annual precipitation over the last 20 years in millimeters (numeric)}
#'   \item{MeanMAPELast20_mm}{Mean annual potential evapotranspiration over the last 20 years in millimeters (numeric)}
#'   \item{StDevMAPELast20_mm}{Standard deviation of annual potential evapotranspiration over the last 20 years in millimeters (numeric)}
#'   \item{MeanMAAELast20_mm}{Mean annual actual evapotranspiration over the last 20 years in millimeters (numeric)}
#'   \item{StDevMAEELast20_mm}{Standard deviation of annual actual evapotranspiration over the last 20 years in millimeters (numeric)}
#' }
#' @source #Morphological descriptors
#'         https://www.academia.edu/5241352/Atlante_bacini_piemontesi_LR
#'         https://sdi.eea.europa.eu/catalogue/copernicus/api/records/66fa7dca-8772-4a5d-9d56-2caba4ecd36a
#'         #Climatic descriptors
#'         #https://www.isprambiente.gov.it/pre_meteo/idro/BIGBANG_ISPRA.html
"basin_descriptors"
