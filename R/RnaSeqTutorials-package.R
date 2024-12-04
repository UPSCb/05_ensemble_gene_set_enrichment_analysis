# RnaSeqTutorial05 package
#
# * To find the imported packages, in the terminal
#
# ---
# cd inst
# grep "library(" */*/*.Rmd | sed -e 's:.*library::g' | tr -d '()' | sort | uniq
# ---
#
# * To build the DESCRIPTION Imports string
#
# ---
# library(here)
# pkgs <- c("dplyr","edgeR","EGSEA","EGSEAdata","here","learnr","limma","org.Hs.eg.db","readr","tibble","UpSetR")
# write(paste0("    ",pkgs," (>= ",unlist(installed.packages()[pkgs,"Version"],use.names=FALSE),"),"),
#       file="Imports.tmp")
# ---
#
#' @title RnaSeqTutorials package
#' @section Tutorials:
#' This is the fifth in a series of tutorials
#' \itemize{
#' \item\code{05_ensemble_gene_set_enrichment_analysis} a tutorial on ensemble gene set enrichment analysis
#' }
#'
#' @name RnaSeqTutorial05 package
#' @rdname RnaSeqTutorial05-package
#' @author Nicolas Delhomme [aut,cre]
#' @keywords package
#' @description A simple description of the RnaSeqTutorial05 package
#' @seealso The vignette
#' @examples
#' 	\dontrun{
#' 	learnr::run_tutorial("05_ensemble_gene_set_enrichment_analysis", package = "RnaSeqTutorial05")
#' 	}
#' @keywords internal
"_PACKAGE"
#'
NULL
