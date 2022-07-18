#' Initialize Reference-Free Cell Mixture Projection
#'
#' Initializes the methylome matrix "Mu" for RefFreeCellMix
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values
#' @param K Number of cell types
#' @param Y.Distance Distance matrix (object of class "dist") to use for clustering.
#' @param Y.Cluster Hiearchical clustering object (from hclust function)
#' @param largeOK OK to calculate distance matrix for large number of subjects? (See details.)
#' @param dist.method Method for calculating distance matrix
#' @param ... Additional parameters for hclust function
#' @details Initializes the methylome matrix "Mu" for RefFreeCellMix by computing the mean methylation (from Y)
#'          over K clusters of Y, determined by the Y.Cluster object.  If Y.Cluster object does not exist, it will be 
#'          created from Y.Distance (using additional clustering parameters if supplied).  If Y.Distance does not exist,
#'          it will be created from t(Y).  As a protection against attempting to fit a very large distance matrix, the
#'          program will stop if the number of columns of Y is > 2500, unless largeOK is explicitly set to TRUE.
#' @return An m x K matrix of mean methylation values.
#' @seealso `RefFreeCellMix`
