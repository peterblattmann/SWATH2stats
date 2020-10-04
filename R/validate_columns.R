#' Validate columns for a data.frame
#'
#' This function looks at the different possible column names given and chooses 
#' the one present in a data.frame. If none of the column names fit or if 
#' multiple names fit the function stops with an appropriate error  message.
#' The functions returns a list with the column names existing that can be used.
#' 
#' @name validate_columns
#' @param data data.frame to check for columns.
#' @param columns List of column names to be checked if they exist. 
#' @param verbose Logical if message should be printed. Default = FALSE
#' @return Returns list of columns that are present
#' @author Peter Blattmann
#' @examples{
#'  validate_columns(cars, list(Speed = c("speed")))
#'  
#'  # if out of two possible column one exists
#'  validate_columns(cars, list(Speed = c("speed", "velocity")))
#'  validate_columns(cars, list(Speed = c("speed", "velocity")), verbose = TRUE)
#'  
#'  }
#' @export
validate_columns <- function(data, columns, verbose = FALSE){
  columns_checked <- columns
  for(i in seq_len(length(columns))){
    if(sum(columns[[i]] %in% colnames(data) == 1)){
      if(length(columns[[i]]) > 1){
        columns_checked[[i]] <- columns[[i]][columns[[i]] %in% colnames(data)]
        if(verbose){
          message(paste0(columns_checked[[i]], " is used as ", names(columns[i]), " column."))
        }
      } else (
        columns_checked[i] <- columns[i])
    } else if(sum(columns[[i]]  %in% colnames(data) == 0)){
      if(is.null(names(columns[i]))){
        stop(paste0("Column ", i, " does not exist in data"))
      } else (stop(paste0(names(columns[i]), " column does not exist in data")))
    } else if(sum(columns[[i]]  %in% colnames(data) > 1)){
      if(is.null(names(columns[i]))){
        stop(paste0("Both names from column ", i, 
                    " exist in data. Please specify which to use"))
      } else (stop(paste0("Both names from the ", names(columns[i]), 
                          " column exist in data. Please specify which to use")))
    }
  }
  return(columns_checked)
}
