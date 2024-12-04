# RnaSeqTutorialXX package
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
# pkgs <- c("learnr","palmerpenguins","tidyverse")
# write(paste0("    ",pkgs," (>= ",unlist(installed.packages()[pkgs,"Version"],use.names=FALSE),"),"),
#       file="Imports.tmp")
# ---
#
#' @title RnaSeqTutorials package
#' @section Tutorials:
#' This is the XX in a series of tutorials
#' \itemize{
#' \item\code{XXX} XXX
#' }
#'
#' @name RnaSeqTutorialXX package
#' @rdname RnaSeqTutorialXX-package
#' @author Nicolas Delhomme [aut,cre]
#' @keywords package
#' @description A simple description of the RnaSeqTutorialXX package
#' @seealso The vignette
#' @examples
#' 	\dontrun{
#' 	learnr::run_tutorial("XXX", package = "RnaSeqTutorialXX")
#' 	}
#' @keywords internal
"_PACKAGE"
#'
NULL
