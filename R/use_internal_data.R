#' Extract lists in the internal environment
#'
#' @param list_name character string of name list
#' @param dataset_name  character string of name data frame 
#'
#' @return return the requested dataframe
#' if no matching is found (list or dataframe), the function stops with an error
#' 
#' @examples
#' # Assuming you have stored stream_stages dataframe inside DA_matrices list
#' # Access availability matrices
#' use_internal_data(list_name = "DA_matrices",dataset_name = "stream_stages")
#' 
#' @export

use_internal_data <- function(list_name,dataset_name) 
{
  if (!exists(list_name, envir = .internal_env))
    stop(paste("List",list_name,"has not been created yet."))
  
  datasets <- .internal_env[[list_name]]  # Retrieve the list
  
  if (!dataset_name %in% names(datasets))
    stop(paste("Dataset",dataset_name, "has not been found in list",list_name))

  return(datasets[[dataset_name]]) 
}
