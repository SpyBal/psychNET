\name{psychNET}
\alias{psychNET}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Psychometric networks estimated by multivariate time series methods.
}
\description{
Wrapper function of multivariate time-series models, which can be used to obtain symptom networks (i.e., networks of symptom-symptom interactions).
}
\usage{
psychNET(data, model, lag, criterion, nFact, penalty, lambda1, lambda2, 
covariates, impute, transform, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{: Usually a \code{"matrix"}, \code{"data.frame"}, or  \code{"longitudinal"} object of dimension \emph{(n*T) x p} where \emph{n} is the number of individuals. For time-series data from multiple persons, an additional numeric column explicitly named \code{"ID"} indicating from which person the measurements are coming from should be included in the data. Additionally, with experience sampling data that are nested within days, additional columns named \code{"DAY"} and \code{"BEEP"} can be included in \code{data} from which we calculate the consecutiveness of the measurements. Alternatively, a column named \code{"TIME"} denoting the measuremnt consecutiveness can be provided in \code{data}. See details. 

}
\item{model}{: This argument controls the  model to be fitted to the data. The available options are: \code{"VAR"}, \code{"SVAR"}, \code{"SMVAR"}, \code{"SVECM"}, \code{"GVAR"}, \code{"SVARHL"}, \code{"SVARHLX"}, \code{"SVARHLMA"}, \code{"DFM"}, \code{"MLVAR"}, \code{"GGVAR"}.
}
\item{lag}{: Specifies the lag order of the process. When \code{model} is \code{"DFM"} \code{lag}  corresponds to the order of the factor process. When \code{model} is \code{"SVARHL"}, \code{"SVARHLX"}, or \code{"SVARHLMA"} \code{lag} can be \code{NULL} and the model chooses the optimal lag via regularized hierarchical lag structures.
}
\item{criterion}{: The information criterion to be used in order to tune the penalty parameters of regularized (sparse) VAR models. Available options are: \code{"CV"}, \code{"AIC"}, \code{"BIC"}, \code{"EBIC"}, \code{"GIC"}, and \code{"MBIC"}. This argument depends on the \code{model} argument (see details).
}

\item{nFact}{: This argument is used only when \code{model="DFM"} and controls the number of the static factors in the dynamic factor model.
}
\item{penalty}{: This argument controls the type of regularization to be used and depends on the \code{model} argument (see details). Available options are \code{"LASSO"}, \code{"ENET"}, \code{"SCAD"}, \code{"MCP"}, and \code{"HLag"}.
}
\item{lambda1}{: A numeric vector with length greater than one that specifies the values of the penalty parameter for regularized VARs. Typically, this corresponds to the penalty parameter for the lagged effects. When \code{model} is \code{"GVAR"} or \code{"GGVAR"} then \code{lambda1} is used as the penalty value for the precision matrix.
}
\item{lambda2}{: A numeric vector with length greater than one. This argument is used only when \code{model} is \code{"GVAR"} or \code{"GGVAR"} and controls the value of the penalty parameter for the lagged effects. 
}
\item{covariates}{: This argument is used only when \code{model} is \code{"VAR"} or \code{"SVARHLX"} in order to specify covariates for VARX or SVARX models.
}
\item{impute}{: String that is used to specify the imputation method when missing values exist in the data or when missing values are inserted to transform non-equidistant measurements into equidistant. It can be a single string or a vector of strings with length equal to the number of symptoms in the data. In the latter case, each symptom can be imputed using an explicit method. Available options are: \code{"Kalman.arima"}, \code{"Kalman.struct"}, \code{"Interpol.linear"}, \code{"Interpol.spline"}, \code{"Interpol.stine"}, \code{"LOCF"}, \code{"NOCB"}, \code{"MA.simple"}, \code{"MA.linear"}, \code{"MA.exponent"}, \code{"mean"}, \code{"mode"}, \code{"median"}, \code{"random"}, \code{"Season.int.linear"}, \code{"Season.int.spline"}, \code{"Season.int.stine"}, \code{"Season.LOCF"}, \code{"Season.NOCB"},\code{"Season.mean"}, \code{"Season.median"}, \code{"Season.mode"}, \code{"Season.MA"}, \code{"Season.kalman"}, \code{"Season.random"}.
}

\item{transform}{: String that is used to specify the transormation function. It can be a single string or a vector of strings with length equal to the number of symptoms in the data. In the latter case, each symptom can be transformed using an explicit function. Available options are: \code{"log"}, \code{"log10"}, \code{"Copula_discr"}, \code{"Copula_skew"}, \code{"Zero.mean"}, \code{"Standardize"}, \code{"Power"}, \code{"Logit"}, \code{"Square.root"}, \code{"Power2"}, \code{"Power3"}, \code{"Cube.root"}.
}

\item{...}{Argument that depends on the value of the \code{model} argument. It is used to pass additional arguments to the model fitting function. 
}

}




\details{
\strong{The data}

The data to be used in the \code{psychNET} function must strictly contain symptom expression data where the measurements are ordered by time. The only additional variables that are allowed in the data except symptom expression are a variable explicitly coded as "ID" when multiple persons need to be analyzed together, a variable coded as "DAY" when multiple measurements per day are taken and a variable called "BEEP" to indicate the measurement intensity within a given day. The variables "DAY" and "BEEP" are used internally to construct a time variable that denotes the consecutiveness of the measurements. The time variable can also be provided in the data explicitly named as "TIME". The only models in the package that can handle non-equidistant observations naturally are the \code{"GVAR"}, \code{"SMVAR"}, and \code{"MLVAR"}. In any other model, non-equidistant measurements first are transformed to equidistant by inserting missing values, which are then imputed using the imputation method specified in the \code{impute} argument.

\strong{The models}

Here we provide some more details with respect to the models implemented in the package. 

\strong{\emph{VAR model} (VAR)}

The traditional VAR model is fitted with the \code{psychNET} function by setting \code{model="VAR"}. This model is originally implememnted in the function \code{VAR} that comes with the R package \pkg{vars}. In this implementation the unknown model parameters are estimated by ordinary least squares (OLS) per equation. The arguments \code{lag} and \code{covariates} of the \code{psychNET} function correspond to the parameter \code{p} and \code{exogen} respectively of the original \code{VAR} function. Other parameters that are used in the original function can be passed through the \code{...} structure (technically known as ellipsis). For details on arguments that can be passed to the three dots when fitting a VAR the user can type \code{?vars::VAR} into the console in order to access its help file.


\strong{\emph{Sparse VAR model} (SVAR)}

A sparse VAR model is fitted with the \code{psychNET} function by setting \code{model="SVAR"}. This model is implememnted in the function \code{fitVAR} from the R package \pkg{sparsevar} and the unknown model parameters are estimated via penalized multivariate least squares. Three different types of penalization are available that can be passed to the function through the \code{penalty} argument. If you set \code{penalty="ENET"} an elasticnet penalty is applied to the parameters with other option being \code{"SCAD"} for SCAD type of penalization and \code{"MCP"} for MC+ penalty. The only option available for tuning the penalty parameter is k-fold cross validation. The parameter \code{lag} of the \code{psychNET} function corresponds to the parameter \code{p} of the \code{fitVAR} function, the parameter \code{penalty} is used in the same way as in the \code{fitVAR} function while \code{criterion} substitutes the argument \code{method} of the \code{fitVAR} function. Additional arguments can be passed to the \code{fitVAR} function by using the \code{...} structure of the \code{psychNET} function. For details on the original function and additional parameters that can be used you can look at \code{?sparsevar::fitVAR}.


\strong{\emph{Sparse VAR model} with hierarchical lag structures and additional covariates (SVARHL and SVARHLX)}

A sparse VAR that offers the possibility of simultaneous lag estimation is fitted by setting \code{model="SVARHL"}. Additional covariates can be used by setting \code{model="SVARHLX"} and providing the covariates data to the argument \code{covariates} of the \code{psychNET} function. These two models are available with the R package \pkg{bigtime} in the functions \code{sparseVAR} and \code{sparseVARX} respectively. These models are estimated by penalized least squares optimized by a proximal gradient algorithm. These models offer two types of penalization and the tuning of the penalty parameters is done by k-fold cross validation. One option is to penalize the VAR coefficients by hierarchical group LASSO regularization and the other option is the standard LASSO also known as L1 regularization. For hierarchical group LASSO the user must set \code{penalty="HLag"} and for standard LASSO \code{penalty="LASSO"}. In the \code{psychNET} function, \code{penalty} argument corresponds to the \code{VARpen} and \code{VARXpen} arguments of the functions \code{sparseVAR} and \code{sparseVARX} respectively, while \code{lag} is associated with the parameter \code{p} of these two functions. Typically, the function internally calculates a grid for the regularization parameter corresponding to sparse penalty. The user can provide her/his own grid of regularization parameters with the argument \code{lambda1}. The latter acts like the arguments \code{VARlseq} and \code{VARXlPhiseq} of the original functions \code{sparseVAR} and \code{sparseVARX} respectively. Other arguments of these functions can be passed via the \code{...} structure of the \code{psychNET} function. If you want the function to estimate the lag of the process simultaneously you must set \code{lag=NULL}. For details on the original functions and additional parameters that can be used you can look at \code{?bigtime::sparseVAR} and \code{?bigtime::sparseVARX}.


\strong{\emph{Sparse VAR model} with simultaneous covariance estimation (GVAR)}

A sparse VAR with simultaneous covariance estimation (also known as graphical VAR) is fitted by setting \code{model="GVAR"}. This model is available with the R package \pkg{graphicalVAR} in the function \code{graphicalVAR}. The model parameters are estimated by otimizing a penalized maximum likelihood. The only penalization option when \code{model="GVAR"} is LASSO, which means that \code{penalty} argument must be equal to \code{"LASSO"}. The tuning of the penalty parameters is done by the extended Bayesian information criterion (EBIC) by setting \code{criterion="EBIC"}. Model selection via standard BIC is possible by using an additional argument \code{gamma=0} in the \code{psychNET} function. This model uses two distinct penalties, one for the autoregressive parameters and one for the covariance parameters. Typically, the function calculates its own grid of penalty parameters for the regularization, however, the user can provide her/his own grid by using the parameters \code{lambda1} (for covariance penalization) and \code{lambda2} (for autoregressive penalization) respectively in the function \code{psychNET}. Other arguments of the original \code{graphicalVAR} function can be passed via the \code{...} structure of the \code{psychNET} function. 

\strong{!!WARNING!!} Do not use the \code{...} to provide the following arguments of the \code{graphicalVAR} function: \code{vars}, \code{scale}, \code{idvar}, \code{beepvar}, \code{dayvar}, \code{centerWithin}, \code{deleteMissings} since they are already specified internally. If you want to use \code{idvar}, \code{beepvar}, \code{dayvar} arguments of the original function \code{graphicalVAR}, these variables must be provided as columns in the \code{data} object explicitly named as "ID", "BEEP", and "DAY" respectively. For \code{scale} and \code{centerWithin} arguments you can use the argument \code{transform} of the \code{psychNET} function. For details on the original function and additional parameters that can be used you can look at \code{?graphicalVAR::graphicalVAR}.

\strong{\emph{Sparse VAR model} for mixed type of time series (SMVAR)}

A sparse VAR with simultaneous covariance estimation (also known as graphical VAR) is fitted by setting \code{model="SMVAR"}. This model is available with the R package \pkg{mgm} in the function \code{mvar}. The model parameters are estimated by otimizing a penalized least squares function per equation. The only penalization option when \code{model="SMVAR"} is elasticnet, which means that \code{penalty} argument must be equal to \code{"ENET"}. The tuning of the penalty parameters is done by the extended Bayesian information criterion (EBIC) by setting \code{criterion="EBIC"} or k-fold cross validation by setting \code{criterion="CV"}. Typically, the function internally calculates a grid for the regularization parameter. The user can provide her/his own grid of regularization parameters with the argument \code{lambda1}. The latter acts analogous to the argument \code{lambdaSeq} of the original \code{mvar} function. Model selection via standard BIC is possible by using an additional argument \code{gamma=0} in the \code{psychNET} function. Other arguments of the original \code{mvar} function can be passed via the \code{...} structure of the \code{psychNET} function.

\strong{!!WARNING 1!!} When \code{model="SMVAR"} the \code{data} must be given in a "data.frame" format with the following properties:
\enumerate{
\item{Gaussian variables}{: must be of class "numeric"}
\item{Poisson variables}{: must be of class "integer"}
\item{Categorical variables}{: must be of class "factor"}
\item{Ordinal variables}{: must be of class "ordered factor" (SMVAR for ordinal type of variables to be expected in a following version of the package \pkg{psychNET}).}
}

\strong{!!WARNING 2!!} Do not use the \code{...} to provide the following arguments of the \code{mvar} function: \code{type}, \code{level}, \code{scale}, \code{beepvar}, \code{dayvar} since they are already specified internally. If you want to use \code{beepvar} and \code{dayvar} arguments of the original function \code{mvar}, these variables must be provided as columns in the \code{data} object explicitly named as "BEEP", and "DAY" respectively. For \code{scale} you can use the argument \code{transform} of the \code{psychNET} function. For details on the original function and additional parameters that can be used you can look at \code{?mgm::mvar}.

\strong{\emph{Sparse VEC model} (SVECM)}

A sparse VEC (vector error correction) model is fitted with the \code{psychNET} function by setting \code{model="SVECM"}. This model is implememnted in the function \code{fitVECM} from the R package \pkg{sparsevar} and the unknown model parameters are estimated via penalized multivariate least squares. Three different types of penalization are available that can be passed to the function through the \code{penalty} argument. If you set \code{penalty="ENET"} an elasticnet penalty is applied to the parameters with other option being \code{"SCAD"} for SCAD type of penalization and \code{"MCP"} for MC+ penalty. The only option available for tuning the penalty parameter is k-fold cross validation. The parameter \code{lag} of the \code{psychNET} function corresponds to the parameter \code{p} of the \code{fitVECM} function, the parameter \code{penalty} is used in the same way as in the \code{fitVECM} function while \code{criterion} substitutes the argument \code{method} of the \code{fitVAR} function. Additional arguments can be passed to the \code{fitVAR} function by using the \code{...} structure of the \code{psychNET} function. For details on the original function and additional parameters that can be used you can look at \code{?sparsevar::fitVECM}. The argument \code{logScale} of the orgiginal function is set to \code{FALSE}, but you can use the argument \code{transform} of the \code{psychNET} function instead.


\strong{\emph{Sparse VARMA model} with hierarchical lag structures (SVARHLMA)}

A sparse VARMA (vector autoregressive moving average) model that offers the possibility of simultaneous lag estimation is fitted by setting \code{model="SVARHLMA"}. This model is available with the R package \pkg{bigtime} in the function \code{sparseVARMA}. The model parameters are estimated by penalized least squares optimized by a proximal gradient algorithm. The model offers two types of penalization and the tuning of the penalty parameters is done by k-fold cross validation. One option is to penalize the VAR coefficients by hierarchical group LASSO regularization and the other option is the standard LASSO. For hierarchical group LASSO the user must set \code{penalty="HLag"} and for standard LASSO \code{penalty="LASSO"}. In the \code{psychNET} function, the \code{penalty} argument corresponds to the \code{VARpen} argument of the function \code{sparseVARMA}, while \code{lag} is associated with the parameter \code{p} of this function. Typically, the function internally calculates a grid for the regularization parameter corresponding to sparse penalty. The user can provide her/his own grid of regularization parameters with the argument \code{lambda1}. The latter acts like the arguments \code{VARlseq} of the original function \code{sparseVARMA}. Other arguments of this function can be passed via the \code{...} structure of the \code{psychNET} function. If you want the function to estimate the lag of the process simultaneously you must set \code{lag=NULL}. For details on the original function and additional parameters that can be used you can look at \code{?bigtime::sparseVARMA}.


\strong{\emph{Dynamic factor model} (DFM)}

A DFM is fitted by setting \code{model="DFM"} in the \code{psychNET} function. This model is available with the R package \pkg{dynfactoR} (available at github) in the function \code{dfm}. The model parameters are estimated by an EM algorithm. The argument \code{nFact} of the \code{psychNET} function is used to specify the number of static factors and it is analogous to the argument \code{r} of the original \code{dfm} function. The argument \code{lag} of the \code{psychNET} function is used to specify the lag order of the factor process analogously to the argument \code{p} of the \code{dfm} function. The DFM that we implement here is restricted to have an identity system state covariance matrix and zeros at the upper diagonal elements of the factor loading matrix for identifiability purposes. The user of the \pkg{psychNET} package can specify her/his own threshold for algorithm convergence by using an additional argument \code{threshold}.


\strong{\emph{Multi-level VAR model} for group level dynamics (MLVAR)}

An MLVAR can be fitted by setting \code{model="MLVAR"}. This model is available with the R package \pkg{mlVAR} in the function \code{mlVAR}. 

\strong{!!WARNING!!} Do not use the \code{...} to provide the following arguments of the \code{mlVAR} function: \code{vars}, \code{scale}, \code{idvar}, \code{beepvar}, \code{dayvar}, \code{estimator}, and \code{scaleWithin} since they are already specified internally. If you want to use \code{idvar}, \code{beepvar}, \code{dayvar} arguments of the original function \code{mlVAR}, these variables must be provided as columns in the \code{data} object explicitly named as "ID", "BEEP", and "DAY" respectively. For the \code{scaleWithin} argument you can use the argument \code{transform} of the \code{psychNET} function. For details on the original function and additional parameters that can be used you can look at \code{?mlVAR::mlVAR}. The only estimation possible is "lmer".

\strong{\emph{Sparse VAR model} with simultaneous covariance estimation for group level dynamics (GGVAR)}

A sparse VAR with simultaneous covariance estimation for multiple individuals (also known as graphical VAR) is fitted by setting \code{model="GGVAR"}. This model is available with the R package \pkg{SparseTSCGM} in the function \code{sparse.tscgm}. The model parameters are estimated by optimizing a penalized maximum likelihood. Penalization options when \code{model="GGVAR"} is LASSO and SCAD, which means that \code{penalty} argument must be equal to \code{"LASSO"} for LASSO penalization and equal to \code{"SCAD"} for SCAD penalization. The tuning of the penalty parameters is done by AIC, BIC, EBIC, GIC (generalized information criterion), and MBIC (modified Bayesian information criterion) by using the \code{criterion} argument which is the same as the \code{optimality} argument in the \code{sparse.tscgm} function. This model uses two distinct penalties, one for the autoregressive parameters and one for the covariance parameters. Typically, the function calculates its own grid of penalty parameters for the regularization, however, the user can provide her/his own grid by using the parameters \code{lambda1} (for covariance penalization) and \code{lambda2} (for autoregressive penalization) respectively in the function \code{psychNET}. Other arguments of the original \code{graphicalVAR} function can be passed via the \code{...} structure of the \code{psychNET} function. Be aware that only VAR models with up to 2 lags are possible.

\strong{Additional details}
The \code{psychNET} function will estimate (based on the lag of the process) one or more temporal conditional independence graphs unless the \code{model} equals to \code{"DFM"}. For \code{model="DFM"} and \code{lag=1} the function will estimate the equivalent VAR(1), whereas for \code{lag} greater than one, the temporal network of the factor interactions will be estimated. Contemporaneous conditional independence graphs are available only when the \code{model} argument equals to: \code{"GVAR"}, \code{"GGVAR"}, \code{"MLVAR"}, or \code{"DFM"} with one lag.

}


\value{
The function \code{psychNET} returns a list object of class \code{pnt}.
 
}

\author{Spyros E. Balafas (author, creator), Sanne Booij, Marco A. Grzegorczyk, Hanneke Wardenaar-Wigman, Ernst C. Wit 

Maintainer: Spyros E. Balafas (\email{s.balafas@rug.nl})
}

\references{
Lutkepohl, H. (2005). \emph{New Introduction to Multiple Time Series Analysis.} Springer, New York.

Basu, S., Michailidis, G. (2015). \emph{Regularized estimation in sparse high-dimensional time series models.} Ann. Statist. 43, no. 4, 1535-1567.

Abegaz, F., Wit, E. (2013). \emph{Sparse time series chain graphical models for reconstructing genetic networks.} Biostatistics. 14, 3: 586-599.

Haslbeck, J., Waldorp, L. J. (2015). \emph{mgm: Structure Estimation for time-varying Mixed Graphical Models in high-dimensional Data.} 

Nicholson, W. B., Bien, J., Matteson, D. S. (2017). \emph{High Dimensional Forecasting via Interpretable Vector Autoregression.}.

Wilms, I., Basu, S., Bien, J., Matteson D. S. (2017). \emph{Sparse Identification and Estimation of High-Dimensional Vector AutoRegressive Moving Averages.} 

Epskamp, S., Waldorp, L. J., Mottus, R.,  Borsboom, D. (2016). \emph{The Gaussian Graphical Model in Cross-sectional and Time-series Data.} 

}


\examples{
## load the psychNET library
library(psychNET)

## load the Canada dataset from the 'vars' package
data("Canada", package = "vars")
Canada_data_frame <- data.frame(Canada)

## fit a VAR model
VAR_model <- psychNET(Canada_data_frame, model = "VAR", lag = 1, type = "const")
# print the result
VAR_model
# summarize the resulting network
summary(VAR_model)
# summarize the VAR model using the original summary method
vars:::summary.varest(VAR_model$fit)
# plot the VAR model results using the original plot method
vars:::plot.varest(VAR_model$fit)
# plot the resulting network
plot(VAR_model)

## fit a sparse VAR model
sparse_VAR_model <- psychNET(Canada_data_frame, model = "SVAR", lag = 1)
# print the result
sparse_VAR_model
# summarize the resulting network
summary(sparse_VAR_model)
# plot the sparse VAR model results using the original plot method
psychNET:::plotVAR(sparse_VAR_model$fit)
# plot the resulting network
plot(sparse_VAR_model)

## fit a sparse VAR model as the one in the 'bigtime' package
sparse_lassVAR_model <- psychNET(Canada_data_frame, 
model = "SVARHL",penalty = "LASSO", lag = 1, VARgran=c(500,1000))
# print the result
sparse_lassVAR_model
# summarize the resulting network
summary(sparse_lassVAR_model)
# plot the resulting network
plot(sparse_lassVAR_model)

\donttest{
## Load the psychNET package
library(psychNET)
#################################### N=1 models ######################################
######################################################################################


####################################################################
######## REPRODUCE EXAMPLE OF A VAR FROM THE 'vars' PACKAGE ######## 
####################################################################

## load the 'vars' package
library(vars)

## Load the Canada dataset from the vars package
data(Canada)

## check the structure of the data
str(Canada)

## The data is in time series format. It needs to be transformed into
## a matrix, data.frame or longitudinal object for the psychNET package
Canada_data_frame <- data.frame(Canada)

## fitting a VAR using the vars package 
varmod <- vars::VAR(Canada, p = 2, type = "none")

## fitting the same VAR using the psychNET package 
psychvar <- psychNET(Canada_data_frame, model = "VAR", lag = 2, type = "none")

## Check if the results are the same 
all.equal(Acoef(varmod), psychvar$results$A, check.attributes = FALSE)

################################################################################
####### Fit A DYNAMIC FACTOR TO THE CANADA DATA FROM THE 'vars' PACKAGE  ####### 
################################################################################

## install and load the 'devtools' package and the 
## by using install.packages("devtools") and then
## library(devtools)
## install the 'dynfactoR' package available on github
## by using devtools::install_github("rbagd/dynfactoR")
## library(dynfactoR)

## Load the Canada dataset from the vars package
data(Canada, package = "vars")

## check the structure of the data
str(Canada)

## The data is in time series format. It needs to be transformed into
## a matrix, data.frame or longitudinal object for the psychNET package
Canada_data_frame <- data.frame(Canada)

## fitting a DFM using the dynfactoR package 
dfmmod <- psychNET:::dfm(Canada_data_frame, r=2, p = 1,
                         rQ= "identity", rC= "upper",max_iter = 100000)

## fitting the same DFM using the psychNET package 
psychdfm <- psychNET(Canada_data_frame, model = "DFM", nFact = 2,  lag = 1)

## Check if the results are the same 
all.equal(dfmmod$A, psychdfm$results$A_fact[[1]], check.attributes = FALSE)
all.equal(dfmmod$C, psychdfm$results$B_fac_symptoms, check.attributes = FALSE)
all.equal(dfmmod$Q, psychdfm$results$System_Covariance, check.attributes = FALSE)
all.equal(dfmmod$R, psychdfm$results$Obs_Covariance, check.attributes = FALSE)


#####################################################################################
########### FIT A SPARSE VAR TO THE CANADA DATA FROM THE 'vars' PACKAGE  ############
#####################################################################################

## Load the Canada dataset from the vars package
data(Canada, package = "vars")

## The data is in time series format. It needs to be transformed into
## a matrix, data.frame or longitudinal object for the psychNET package
Canada_data_frame <- data.frame(Canada)

## fitting a SVAR using the sparsevar package
set.seed(1)
svarmod <- psychNET:::fitVAR(Canada, p = 1, penalty="SCAD", method="cv", logScale=FALSE)

## fitting the same SVAR using the psychNET package
set.seed(1)
psychsvar <- psychNET(Canada_data_frame, model = "SVAR",penalty = "SCAD", lag = 1, criterion="CV")

## Check if the results are the same
all.equal(svarmod$A, psychsvar$results$A, check.attributes = FALSE)


########################################################################################
####### EXAMPLE OF A SPARSE VAR WITH HIERARCHICAL LAGS FROM THE 'bigtime' PACKAGE ######
########################################################################################

## load the 'bigtime' package
library(bigtime)
## Load the Y dataset from the bigtime package
data(Y, package = "bigtime")

## fitting a SVAR with hierarchical lags
## using the bigtime package
svarhlmod <- sparseVAR(Y, VARpen ="HLag")

## fitting the same model using the psychNET package
psychsvarhl <- psychNET(Y, model = "SVARHL", penalty = "HLag", criterion="CV")

## Check if the results are the same
all.equal(svarhlmod$Phihat, psychsvarhl$fit$Phihat, check.attributes = FALSE)



######################################################################################
########### REPRODUCE EXAMPLE OF A SPARSE MIXED VAR FROM THE 'mgm' PACKAGE ########### 
######################################################################################

## load the 'mgm' package
library(mgm)

# 1) Define mVAR model as in the mgm manual
p <- 6 # Six variables
type <- c("c", "c", "c", "c", "g", "g") # 4 categorical, 2 gaussians
level <- c(2, 2, 4, 4, 1, 1) # 2 categoricals with m=2, 2 categoricals with m=4, two continuous
max_level <- max(level)
lags <- 1:3 # include lagged effects of order 1-3
n_lags <- length(lags)

# Specify thresholds
thresholds <- list()
thresholds[[1]] <- rep(0, level[1])
thresholds[[2]] <- rep(0, level[2])
thresholds[[3]] <- rep(0, level[3])
thresholds[[4]] <- rep(0, level[4])
thresholds[[5]] <- rep(0, level[5])
thresholds[[6]] <- rep(0, level[6])

# Specify standard deviations for the Gaussians
sds <- rep(NULL, p)
sds[5:6] <- 1

# Create coefficient array
coefarray <- array(0, dim=c(p, p, max_level, max_level, n_lags))

# a.1) interaction between continuous 5<-6, lag=3
coefarray[5, 6, 1, 1, 2] <- .4
# a.2) interaction between 1<-3, lag=1
m1 <- matrix(0, nrow=level[2], ncol=level[4])
m1[1,1:2] <- 1
m1[2,3:4] <- 1
coefarray[1, 3, 1:level[2], 1:level[4], 1] <- m1
# a.3) interaction between 1<-5, lag=9
coefarray[1, 5, 1:level[1], 1:level[5], 3] <- c(0, 1)


# 2) Sample
set.seed(1)
dlist <- mvarsampler(coefarray = coefarray,
                     lags = lags,
                     thresholds = thresholds,
                     sds = sds,
                     type = type,
                     level = level,
                     N = 200,
                     pbar = TRUE)

# 3) Transform data into a data.frame for psychoNET suitability
# each categorical variable is coded as factor, each poisson as integer, gaussian as numeric
d1 <- as.data.frame(dlist$data)
d1$V1 <- as.factor(d1$V1)
d1$V2 <- as.factor(d1$V2)
d1$V3 <- as.factor(d1$V3)
d1$V4 <- as.factor(d1$V4)
dat <- d1

## fitting the SMVAR model using the mgm package 
smvarmod <- mvar(data = dlist$data,
                 type = type,
                 level = level,
                 lambdaSel = "EBIC",
                 lags = 1:3,
                 scale = FALSE,
                 signInfo = FALSE,
                 overparameterize = FALSE)

## fitting the same model using the psychNET package 
psychsmvar <- psychNET(dat, model = "SMVAR",
                       lag = 3,
                       criterion = "EBIC",
                       signInfo = FALSE,
                       overparameterize = FALSE)


all.equal(smvarmod$wadj,psychsmvar$fit$wadj)


################################### N>1 models ########################################
#######################################################################################

######################################################################################
########### REPRODUCE EXAMPLE OF A MULTILEVEL VAR FROM THE 'mlVAR' PACKAGE ########### 
######################################################################################

## load the 'mlVAR' package
library(mlVAR)

## Simulate data:
Model <- mlVARsim(nPerson = 50, nNode = 3, nTime = 50, lag=1)

# Estimate an MLVAR with correlated random effects using the mlVAR package
mlvarmod <- mlVAR(Model$Data, vars = Model$vars, 
          idvar = Model$idvar, lags = 1, temporal = "correlated")

# Estimate the same MLVAR using the psychNET package
psychmlvar <- psychNET(Model$Data, model="MLVAR",  lag = 1, temporal = "correlated")

# Check if the results are equal
all.equal(mlvarmod$results, psychmlvar$fit$results)

#######################################################################################
########### REPRODUCE EXAMPLE OF A GGVAR VAR FROM THE 'SparseTSCGM' PACKAGE ########### 
#######################################################################################

## load the 'SparseTSCGM' package
library(SparseTSCGM)

## Simulate data:
seed = 321
datas <- sim.data(model="ar1", time=10,n.obs=10, n.var=5,seed=seed,
                  prob0=0.35, network="random")
data.fit <-  datas$data1

# Estimate a group graphical VAR (also known as time series chain graphical model)  
# using the SparseTSCGM package
ggvarmod <- sparse.tscgm(data=data.fit, 
                         lam1=NULL, lam2=NULL, nlambda=NULL, 
                         model="ar1", penalty="scad",optimality="bic",
                         control=list(maxit.out = 10, maxit.in = 100))

# Estimate the same model using the psychNET package
psychggvar <- psychNET(data.fit, model="GGVAR", lag=1, penalty="SCAD", criterion="BIC",
                       control=list(maxit.out = 10, maxit.in = 100))

# Check if the results are equal
all.equal(ggvarmod$theta, psychggvar$fit$theta, check.attributes = FALSE)
all.equal(ggvarmod$gamma, psychggvar$fit$gamma, check.attributes = FALSE)

}
}

