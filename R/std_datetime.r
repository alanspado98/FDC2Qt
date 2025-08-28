#' Standardize Date-Time Strings
#'
#' Attempts to parse a vector of date or date-time strings in lubridate formats and convert them to 
#' standardized POSIXct objects in GMT timezone
#'
#' @param datetime_strings a character vector of date or date-time strings in potentially mixed formats. 
#'   Supported formats include:
#'   - `"yyyy-mm-dd"` (ISO-style)
#'   - `"dd-mm-yyyy"` (European-style)
#'   - `"mm-dd-yyyy"` (US-style)
#'   - `"yyyy-mm-dd HH:MM"` or `"yyyy-mm-dd HH:MM:SS"` (ISO-style)
#'   - `"dd-mm-yyyy HH:MM"` or `"dd-mm-yyyy HH:MM:SS"` (European-style)
#'   - `"mm-dd-yyyy HH:MM"` or `"mm-dd-yyyy HH:MM:SS"` (US-style)
#'
#' @return a `POSIXct` vector in `"GMT"` timezone, with `NA` for unrecognized or malformed date strings
#'
#' @details 
#' Use a predefined list of parsers from the `lubridate` package; parsing is silent with respect to warnings and errors, and the first successful match is returned
#'
#' @importFrom lubridate ymd ymd_hm ymd_hms dmy dmy_hm dmy_hms mdy mdy_hm mdy_hms
#'
#' @examples
#' std_datetime(c("2024-05-01 12:30", "01-05-2024 13:45:00", "May 1, 2024"))
#' std_datetime(c("2024-05-01", "01/05/2024", "05/01/2024"))
#'
#' @export

std_datetime <- function(datetime_strings) 
{
  output <- rep(as.POSIXct(NA, tz = "GMT"), length(datetime_strings))
  
  formats <- list(lubridate::ymd_hm, lubridate::ymd_hms,
                  lubridate::dmy_hm, lubridate::dmy_hms,
                  lubridate::mdy_hm, lubridate::mdy_hms,
                  lubridate::ymd,lubridate::dmy, lubridate::mdy)
  
  for (parser in formats) 
  {
    parsed <- try(suppressWarnings(parser(datetime_strings)), silent = TRUE)
    if (!inherits(parsed, "try-error")) 
    {
      matched <- is.na(output) & !is.na(parsed)
      output[matched] <- as.POSIXct(parsed[matched], tz = "GMT")
    }
  }
  
  return(output)
}
