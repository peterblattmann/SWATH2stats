## I am assuming but am not certain that when R byte-compiles this directory, it
## will read this first and create its environment with this material first,
## since I prefixed it with '01'.  I should probably test this and see if it is
## true.  In any event, maybe it does not matter as the important stuff is in a
## .onLoad?

## The following was taken from ggplot2's ggplot2.r
## I presume it is a blanket importer cue for roxygen2 to add
## import statements to the NAMESPACE file so that when ggplot2 is used
## it will ensure that these libraries are available.
## I checked the roxygen documentation and it appears that
## imports are saved as the exclusive set, as a result repeating these
## at each function declaration serves to make explicit what each function
## requires while not (I think) adding excessive cruft to the NAMESPACE

#' SWATH2stats: a suite of tools to filter, plot, and convert DIA data.
#'
#' @docType package
#' @name SWATH2stats
#' @importFrom ggplot2 aes aes_string ggplot
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis grid legend lines mtext par plot points
#' @importFrom stats aggregate cor density median na.omit sd
#' @importFrom utils head select.list write.csv write.table
NULL

#' data.table's funky column assignment operator
#'
#' Shamelessly scabbed from Hadley: https://github.com/sckott/analogsea/issues/32
#'
#' @name :=
#' @rdname column_assignment
#' @keywords internal
#' @export
#' @importFrom data.table :=
NULL

#' The following is used to set the ggplot2 default text size.
base_size <- 16
