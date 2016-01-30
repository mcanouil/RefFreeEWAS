\name{RefFreeCellMixInitialize}
\alias{RefFreeCellMixInitialize}
\title{Initialize Reference-Free Cell Mixture Projection}
\description{Initializes the methylome matrix "Mu" for RefFreeCellMix}
\usage{RefFreeCellMixInitialize(Y,K=2,Y.Distance=NULL, Y.Cluster=NULL, 
    largeOK=FALSE, dist.method = "euclidean", ...)}
\arguments{
  \item{Y}{Matrix (m CpGs x n Subjects) of DNA methylation beta values}
  \item{K}{Number of cell types}
  \item{Y.Distance}{Distance matrix (object of class "dist") to use for clustering.}
  \item{Y.Cluster}{Hiearchical clustering object (from hclust function)}
  \item{largeOK}{OK to calculate distance matrix for large number of subjects? (See details.)}
  \item{dist.method}{Method for calculating distance matrix}
  \item{\dots}{Additional parameters for hclust function}
}
\details{Initializes the methylome matrix "Mu" for RefFreeCellMix by computing the mean methylation (from Y)
over K clusters of Y, determined by the Y.Cluster object.  If Y.Cluster object does not exist, it will be 
created from Y.Distance (using additional clustering parameters if supplied).  If Y.Distance does not exist,
it will be created from t(Y).  As a protection against attempting to fit a very large distance matrix, the
program will stop if the number of columns of Y is > 2500, unless largeOK is explicitly set to TRUE.}
\value{
  An m x K matrix of mean methylation values.
}
\references{Houseman EA, Kile ML, et al., Reference-free deconvolution of DNA methylation data and mediation by cell composition effects (2016). http://biorxiv.org/content/early/2016/01/23/037671.
}
\author{E. Andres Houseman}
\seealso{\code{\link{RefFreeCellMix}}}

