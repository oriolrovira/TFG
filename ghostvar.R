#' @title Understading complex predective models with Ghost Variables
#' @description
#' Procedure for assigning a relevance measure to each explanatory variable in a complex predictive model.
#' We have a training set to fit the model and a test set to check the out of sample performance.
#' The individual relevance of each variable is computed by comparing in the test set the predictions of the model that includes all the variables,
#' with those of another model in which the variable of interest is substituted by its ghost variable,
#' defined as the prediction of this variable by using the rest of explanatory variables.
#' Furthermore, we check the joint relevance of the variables by using the eigenvalues of a relevance matrix
#' that is the covariance matrix of the vectors of individual effects.
#' @name ghostvar
#' @author \emph{Authors:} Pedro Delicado \email{pedro.delicado@@upc.edu}
#' and Daniel Pe\~na \email{daniel.pena@@uc3m.es}
#' \emph{Contributors:} Oriol Rovira \email{oriol-rovira@@hotmail.com}
#' \emph{Maintainer:} Oriol Rovira \email{oriol-rovira@@hotmail.com}
#' @references Delicado, P., Pe\~na, D. (2019).
#' \emph{Understanding complex predictive models with Ghost Variables}
#' \url{"https://arxiv.org/a/0000-0003-3933-4852.html"}
#' @keywords package
#' @import mgcv
#' @import ggplot2
#' @import grid
#' @import maptools
NULL
#' @export