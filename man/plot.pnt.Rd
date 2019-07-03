\name{plot.pnt}
\alias{plot.pnt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for S3 class "pnt" objects.
}
\description{
Generic that plots an object of class "pnt".
}
\usage{
\method{plot}{pnt}(x, type, person, community, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{:
An object resulted from the \code{psychNET} function.
}
\item{type}{:
String argument, which controls the type of plot that will be returned. The available options are: "temporal" (for temporal network), "contemporaneous" (for contemporaneous network), "between" (for between subjects network), and "both" (both temporal and contemporaneous networks)
}
\item{person}{:
This can be a single number or vector of numbers that denotes the person index for which plots will be returned. This argument is used only when the model fitted is a Multi-level VAR (i.e. the argument \code{model} in the \code{psychNET} function equals to "MLVAR".)
}
\item{community}{:
Logical argument. When \code{TRUE}, the function fits a spinglass community detection algorithm with negative weights where only the present edges are taken into account. !!WARNING!! this can make the \code{plot.pnt} function to be very slow.
}
\item{...}{:
Other arguments to be passed on to the \code{plot.igraph} function. Use this with care since some arguments in are already specified internally.
}
}
\details{
For  This function is used to visualize networks estimated via the \code{psychNET} function. See additional details of the function \code{psychNET}.

}
\value{
The value returned by the plot is a list where its elements are two lists named temporal and contemporaneous respectively. These lists contain objects of class \code{igraph} and \code{qgraph} that can be used by the user to create tailor made plots.
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
\author{Spyros E. Balafas (author, creator), Sanne Booij, Marco A. Grzegorczyk, Hanneke Wardenaar-Wigman, Ernst C. Wit 

Maintainer: Spyros E. Balafas (\email{s.balafas@rug.nl})
}





