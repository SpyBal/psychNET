\name{psychNET-package}
\alias{psychNET-package}
\docType{package}
\title{psychNET : Psychometric network modelling for multivariate time series data.
}
\description{
This package has been designed in order to provide various psychometric network modelling techniques for multivariate time series data from the behavioral domain, in a simple wrapper function. The aim is to estimate the temporal and contemporaneous interaction structure in a (possibly) sparse fashion and thereby being able to visualize these interactions by means of conditional independence graphs. Depending on the model, the interactions are estimated at the individual level -- intra-individual dynamics -- or population level -- inter-individual dynamics.
}
\details{
In the last decade, time series data became popular in the behavioral sciences. These data, allow us to study the dynamics of complex behavioral systems both at the individual and the population level. Vector autoregression (VAR) is the cornerstone in the statistical modelling of multivariate time series data and various VAR extensions are available that can handle cases where additional complexities are imposed in the analysis.

This package introduces the main \code{psychNET} function that is used to fit various dynamic models. Models that can be fitted to time-series data from one person using the \code{psychNET} function are:
\itemize{
\item{Vector Autoregressive (VAR) models}{: Traditional VAR model or VAR model with exogenous variables (VARX). See the help file from the \code{VAR()} function from the R package \pkg{vars} for more details.}
\item{Sparse Vector Autoregressive (SVAR) models}{: SVAR model estimated using penalized multivariate least squares (see \code{fitVAR} function in the R package \pkg{sparsevar}). SVARs estimated using penalized least squares with simultaneous lag estimation (SVARHL) and exogenous variables (SVARHLX) implemented in the functions \code{sparseVAR} and \code{sparseVARX} respectively from the R package \pkg{bigtime}. SVAR model with simultaneous covariance estimation using penalized likelihood approach known as graphical VAR (GVAR) implemented in the function \code{graphicalVAR} from the package \pkg{graphicalVAR}. A SVAR model for mixed type of time series with marginal distributions in the exponential family (SMVAR), implemented in  the function \code{mvar} from the R package \pkg{mgm}. }

\item{Sparse Vector Error Correction Model (SVECM)}{: SVECM estimated using penalized multivariate least squares that is implemented in the function \code{fitVECM} from the R package \pkg{sparsevar}.
}

\item{Sparse Vector Autoregressive Moving Average (SVARMA) model}{: SVARMA model where the parameters are estimated using penalized least squares with simultaneous lag estimation (SVARMAHL) implemented in the function \code{sparseVARMA} from the R package \pkg{bigtime}.
}

\item{Sparse Time-Varying Vector Autoregressive (TV-SVAR) model}{: TV-SVAR model where the parameters are allowed to vary smoothly with time. This model can handle time series with marginal distributions in the exponential family. The model is implemented in the function \code{tvmvar}, which is available with the R package \pkg{mgm}. (This model is NOT yet available in the \pkg{psychNET} package. To be expected with the next update.)
}

\item{Dynamic factor model (DFM)}{: DFM with three different types of estimates for the latent series. For details you can check the function \code{dfm} from the R package \pkg{dynfactoR} available at github.}
}

For network inference at the population level from nested time-series data, two models in the broad class of VAR models can be fitted via the \code{psychNET} function. These are:
\itemize{
\item{Multilevel Vector Autoregressive (MLVAR) model}{: This model estimates group as well as individual level effects and it is implemented in the R function \pkg{mlVAR} from the R package \pkg{mlVAR}.}

\item{Sparse Vector Autoregressive (SVAR) model}{: This model estimates sparse group level networks from nested time series with simultaneous covariance estimation. We call this model group graphical VAR (GGVAR) and it has been primarily implemented in the function \code{sparse.tscgm} available with R package \pkg{SparseTSCGM}.}

}
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

\examples{
# Load the R package psychNET
library(psychNET)
}


\keyword{package}
