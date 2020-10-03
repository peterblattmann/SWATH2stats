#' Select alternate m_score column in JPP data and avert user
#' 
#' The output from JPP (Rosenberger, Bludau et al. 2017) has not anymore the 
#' m_score column but an ProteinName_m_score and transition_group_id_m_score. 
#' To make the users aware this function tests if the m_score column still 
#' exists and selects as default the transition_group_id_m_score column.
#' 
#' @param data Data table that is produced by the OpenSWATH/pyProphet workflow
#' @param mscore_col Defines the column from which to retrieve the m_score. 
#'    If you use JPP (Rosenberger, Bludau et al. 2017) this can be used to 
#'    select between Protein and transition_group m_score.
#' @return Returns the mscore_col that might have been changed to 
#'    transition_group_id_m_score and gives a message to the user.
#' @author Peter Blattmann
#' @examples{
#' data("OpenSWATH_data", package="SWATH2stats")
#' JPP_update(OpenSWATH_data, "m_score")
#' }
#' @export
JPP_update <- function(data, mscore_col) {
    if (!(mscore_col %in% colnames(data))) {
        if ("transition_group_id_m_score" %in% colnames(data)) {
            mscore_col <- "transition_group_id_m_score"
            message("The transition_group_id_m_score was used as mscore. 
                    If you want to use Protein q-value please change option 
                    in the function.")
        } else {
            stop("Column for m_score column doesn't exist in data.")
        }
    }
    return(mscore_col)
}
