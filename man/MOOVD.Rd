\name{MOOVD}
\alias{MOOVD}
\docType{data}
\title{
Multivariate symptom time series from an individual in the Mood and Movement in Daily life (MOOVD) study.
}
\description{
Discrete valued multivariate symptom time series from a depressed individual in the MOOVD study. The responses of the participant are recorded for 30 days, three times per day (one vector of responses obtained at each beep).
}
\usage{data("MOOVD")}
\format{
  A data frame with 90 incomplete observations on the following 30 variables.
  \describe{
    \item{\code{Appetite}}{: ordered factor. "My apetite is..." [-3,3]}
    \item{\code{Psychomotor retardation}}{: ordered factor. " My movements are... " [-3,3]}
    \item{\code{Interest}}{: ordered factor. "My interest in the things around me is..." [-3,3]}
    \item{\code{Talkative}}{: ordered factor. "I feel talkative " [1,7]}
    \item{\code{Energetic}}{: ordered factor. "I feel energetic..." [1,7]}
    \item{\code{Tense}}{: ordered factor. "I feel tense..." [1,7]}
    \item{\code{Anxious}}{: ordered factor. "I feel anxious..." [1,7]}
    \item{\code{Enthusiastic}}{: ordered factor. "I feel enthousiastic..." [1,7]}
    \item{\code{Self assured}}{: ordered factor. "I feel self assured..." [1,7]}
    \item{\code{Distracted}}{: ordered factor. "I feel self distracted..." [1,7]}
    \item{\code{Resteless}}{: ordered factor. "I feel self restless..." [1,7]}
    \item{\code{Irritated}}{: ordered factor. "I feel self irritated..." [1,7]}
    \item{\code{Content}}{: ordered factor. "I feel content..." [1,7]}
    \item{\code{Happy}}{: ordered factor. "I feel cheerful..." [1,7]}
    \item{\code{Sad}}{: ordered factor. "I feel sad..." [1,7]}
    \item{\code{Cheerful}}{: ordered factor. "I feel excited..." [1,7]}
    \item{\code{Guilty}}{: ordered factor. "I feel guilty..." [1,7]}
    \item{\code{Defeated}}{: ordered factor. "I feel defeated..." [1,7]}
    \item{\code{Tired}}{: ordered factor. "I feel tired..." [1,7]}
    \item{\code{Relaxed}}{: ordered factor. "I feel relaxed..." [1,7]}
    \item{\code{Hopeless}}{: ordered factor. "I feel hopeless..." [1,7]}
    \item{\code{Appreciated}}{: ordered factor. "I feel appreciated..." [1,7]}
    \item{\code{Thoughts about death}}{: ordered factor. "I think about death a lot..." [1,7]}
    \item{\code{Can handle what comes my way}}{: ordered factor. "I can handle what comes my way..." [1,7]}
    \item{\code{Rumination}}{: ordered factor. "I worry a lot..." [1,7]}
    \item{\code{Negative thoughts}}{: ordered factor. "Many negative thoughts go through my mind..." [1,7]}
    \item{\code{Concentrate}}{: ordered factor. "I am able to concentrate very well..." [1,7]}
    \item{\code{Slowed thoughts}}{: ordered factor. " My thoughts are slow..." [1,7]}
    \item{\code{DAY}}{: numeric. Day variable }
    \item{\code{BEEP}}{: integer. Beep variable}
  }
}


\references{
Booij,  S.  H.,  Bos,  E.  H.,  Bouwmans,  M.  E.  J.,  van  Faassen,  M.,  Kema,  I.  P.,  Oldehinkel,A.  J.,  and  de  Jonge,  P.   (2015).   Cortisol  and a-amylase  secretion  patterns  between and within depressed and non-depressed individuals.PlOS  ONE,10(7), 1-15.  Retrieved from https://doi.org/10.1371/journal.pone.0131002
}
\examples{
data(MOOVD)
str(MOOVD)
}
\keyword{datasets}
