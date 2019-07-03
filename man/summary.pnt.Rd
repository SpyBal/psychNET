\name{summary.pnt}
\alias{summary.pnt}
\title{
Summary method for S3 class "pnt" objects.
}
\description{
Generic function that is used in order to summarize information from "pnt" class objects.
}
\usage{
\method{summary}{pnt}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{:
An object obtained by the \code{psychNET} function.
}
\item{...}{:
Not used in this version of the package.
}
}
\details{
This a generic function that summarize the information of the networks obtained after a time series model has been fitted to data. Since the main function \code{psychNET} is a wrapper of several models the summary methods of each method can also be used (if available) by typing \code{summary(object$fit)} where \code{object} is an object obtained by  the \code{psychNET} function.

}
\value{
The function \code{summary.pnt} returns a list with the following components:

\item{contemporaneous}{: When estimates of the covariance matrix (contemporaneous network) are available then contemporaneous is a list with two components named local and global that describe the contemporaneous network locally and globally. For global description of the network, transitivity, reciprocity, distance, density, and diameter are returned in a one column matrix. At the local (nodes) level the function calculates node transitivity, the degree centrality, step-1 and -2 node expected influence, betweeness centrality and closeness centrality.}

\item{temporal}{A list with two components named local and global that describe the temporal network locally and globally. At the global graph level the same descriptives as in the contemporaneous network are returned for each lag of the VAR model. At the local level, node transitivity, in and out degree centrality, step -1 and -2 expected influence centralities, betweeness, out and in closeness centralities are returned to the user.}
%% 
}
\author{Spyros E. Balafas (author, creator), Sanne Booij, Marco A. Grzegorczyk, Hanneke Wardenaar-Wigman, Ernst C. Wit 

Maintainer: Spyros E. Balafas (\email{s.balafas@rug.nl})
}

\references{
Lutkepohl, H. (2006). \emph{New Introduction to Multiple Time Series Analysis.} Springer, New York.

Basu, S., Michailidis, G. (2015). \emph{Regularized estimation in sparse high-dimensional time series models.} Ann. Statist. 43, no. 4, 1535-1567.

Abegaz, F., Wit, E. (2013). \emph{Sparse time series chain graphical models for reconstructing genetic networks.} Biostatistics. 14, 3: 586-599.

Haslbeck, J., Waldorp, L. J. (2016). \emph{mgm: Structure Estimation for time-varying Mixed Graphical Models in high-dimensional Data.} 

Nicholson, W. B., Bien, J., Matteson, D. S. (2017). \emph{High Dimensional Forecasting via Interpretable Vector Autoregression.}.

Wilms, I., Basu, S., Bien, J., Matteson D. S. (2017). \emph{Sparse Identification and Estimation of High-Dimensional Vector AutoRegressive Moving Averages.} 

Epskamp, S., Waldorp, L. J., Mottus, R.,  Borsboom, D. (2016). \emph{The Gaussian Graphical Model in Cross-sectional and Time-series Data.} 

}


