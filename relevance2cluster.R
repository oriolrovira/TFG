#' @title Relevance cluster
#' @name relevance2cluster
#' 
#' @description \code{\link{relevance2cluster}} make a plot by using distance matrix from similarity matrix \code{V} of ward method.
#' 
#' @return
#' \code{relevance2cluster} Hierarchical Clustering. A dissimilarity structure as produced by dist (ward.D2)
#' 
#' @param V relevance matrix.
#' @param method value that will be passed to the function \code{hclust}.
#' @param ... additional graphical parameters passed to the function plot. 
#'
#' @seealso \code{\link{relev.ghost.var}} and \code{\link{relev.rand.perm}}
#' @rdname relevance2cluster
#' @export 

relevance2cluster <- function(V,method="ward.D2",...){
  dV <- diag(V)
  p<- length(dV)
  one <- 1+numeric(p)
  # Distance matrix from similarity matrix V:
  W <- as.dist(sqrt(one %*% t(dV) + dV %*% t(one) -2*V))
  hcl.W <- hclust(W, method = method)
  plot(hcl.W,...)
}