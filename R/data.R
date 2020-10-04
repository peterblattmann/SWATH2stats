#' @name MSstats_data
#' @docType data
#' @title Testing dataset in MSstats format. 
#' @description A small table with the column names corresponding to the MSstats format. 
#'  This data is intended only to test functions. 
#' @author Peter Blattmann
NULL

#' @name OpenSWATH_data
#' @aliases data
#' @docType data
#' @title Testing dataset from OpenSWATH.
#' @description A small selection of the data obtained from the iPortal pipeline for an 
#'  SWATH/DIA experiment with perturbations relating to cholesterol regulation. 
#'  Protein and Peptides have been anonymized as the data is unpublished.
#'  The FDR version of the test data contains modified (lowered) decoy peak 
#'  group m_scores to simulate FDR behaviour of a large dataset.
#' @author Peter Blattmann
NULL

#' @name Spyogenes
#' @docType data
#' @title S.pyogenes example data.
#' 
#' @description A table containing SWATH-MS data from S.pyogenes.
#'  This table was generated from the original data deposited on PeptideAtlas 
#'  (PASS00289, file "rawOpenSwathResults_1pcnt_only.tsv") by selecting only 
#'  the column necessary for the SWATH2stats. 
#'
#' @references Rost, H. L., et al. (2014). OpenSWATH enables automated, 
#'   targeted analysis of data-independent acquisition MS data. 
#'   Nat Biotechnol 32(3): 219-223. 
NULL


#' @name Study_design
#' @docType data
#' @title A table containing the meta-data defining the study design of the OpenSWATH data.
#' 
#' @description This table contains a unique identifier in the column "Filename"
#'  corresponding to the filename in the SWATH data.
#'  In the column "Condition" the perturbation performed is described.
#'  In the column "BioReplicate" the biological replicate is indicated.
#'  In the column "Run" a unique identifier for each injection.
#'  Technical injections would have different Run numbers but the same 
#'  BioReplicate number.
#' 
#' @author Peter Blattmann
NULL
