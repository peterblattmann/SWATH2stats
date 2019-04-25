JPP_update <- function(data, mscore.col) {
    if (!(mscore.col %in% colnames(data))) {
        if ("transition_group_id_m_score" %in% colnames(data)) {
            mscore.col <- "transition_group_id_m_score"
            message("The transition_group_id_m_score was used as mscore. If you want to use Protein q-value please change option in the function.")
        } else {
            stop("Column for m_score column doesn't exist in data.")
        }
    }
    return(mscore.col)
}
