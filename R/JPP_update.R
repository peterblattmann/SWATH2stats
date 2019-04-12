#' Select alternate m_score column in JPP data and avert user.
#'
#' Thr output from JPP (Rosenberge, Bludau et al. 2017) no longer has the
#' m_score column, but a ProteinName_m_score and
#' transition_group_id_m_score.  To make users aware of this, this function
#' tests if the m_score still exists and selects as default the
#' transition_group_id_m_score column.
#'
#' @param data Data table produced by the OpenSWATH/pyProphet workflow.
#' @param mscore.col The column from which to retrieve the m_score data.  If one uses JPP
#' (Rosenberger, Bludau et al. 2017) this can be used to select between Protein
#' and transition_group m_score.
#' @return The mscore column that may have been changed to
#'   transition_group_id_m_score.  Also this returns a message to the user.
#' @author Peter Blattmann
#' @examples
#'   data("OpenSWATH_data", package="SWATH2stats")
#'   JPP_update(OpenSWATH_data, "m_score")
#' @export
JPP_update <- function(data, mscore.col) {
  if(!(mscore.col %in% colnames(data))) {
    if("transition_group_id_m_score" %in% colnames(data)) {
      mscore.col <- "transition_group_id_m_score"
      message("The transition_group_id_m_score was used as mscore.
If you want to use Protein q-value please change option in the function.")
    } else {
      stop("Column for m_score column doesn't exist in data.")
    }
  }
  return(mscore.col)
}
