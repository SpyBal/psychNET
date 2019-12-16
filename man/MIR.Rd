\name{MIR}
\alias{MIR}
\docType{data}
\title{
Multivariate symptom time series from the Mapping Individual Routes Of Risk and Resilience (MIRORR) study.
}
\description{
The sample comprises 51 individuals of 18-35 years from two subsamples (n1=25 and n2=26 respectively). Each
subgroup, has an increasingly more severe psychopathological level and thus representing subsequent
clinical stages. 

}
\usage{data("MIR")}
\format{
  A data frame with 4488 observations on the following 36 variables.
  \describe{
    \item{\code{ID}}{: An integer vector with identification numbers for each participant.}
    \item{\code{GROUP}}{: An integer vector with identification numbers for each subgroup of participants.}
    \item{\code{TIME}}{: An integer vector with the time point when each measurement has been taken..}
    \item{\code{relx }}{: numeric. "I felt relaxed today." [0,100]}
    \item{\code{calm }}{: numeric. "I felt calm today." [0,100]}
    \item{\code{satf }}{: numeric. "I felt satisfied today." [0,100]}
    \item{\code{ener }}{: numeric. "I felt energetic today." [0,100]}
    \item{\code{enth }}{: numeric. "I felt enthousiastic today." [0,100]}
    \item{\code{chrf }}{: numeric. "I felt cheerful today." [0,100]}
    \item{\code{apat }}{: numeric. "I felt apathetic today." [0,100]}
    \item{\code{tire }}{: numeric. "I felt tired today." [0,100]}
    \item{\code{down }}{: numeric. "I felt down today." [0,100]}
    \item{\code{anxi }}{: numeric. "I felt anxious today." [0,100]}
    \item{\code{rstl }}{: numeric. "I felt restless today." [0,100]}
    \item{\code{irtb }}{: numeric. "I felt irritable today." [0,100]}
    \item{\code{irtt }}{: numeric. "I felt irritated today." [0,100]}
    \item{\code{talk }}{: numeric. "I felt talkative today." [0,100]}
    \item{\code{conf }}{: numeric. "I felt confident today." [0,100]}
    \item{\code{empt }}{: numeric. "I felt empty today." [0,100]}
    \item{\code{worr }}{: numeric. "I felt worried today." [0,100]}
    \item{\code{susp }}{: numeric. "I felt suspicious today." [0,100]}
    \item{\code{incm }}{: numeric. "Today I had the feeling of falling short." [0,100]}
    \item{\code{hand }}{: numeric. "Today I could handle what came my way." [0,100]}
    \item{\code{conc }}{: numeric. "I could concentrate well today." [0,100]}
    \item{\code{life }}{: numeric. "I found my life was worthwhile today." [0,100]}
    \item{\code{phys }}{: numeric. "I was bothered by physical symptoms today." [0,100]}
    \item{\code{unrs }}{: numeric. "Today I had the tendency to do something wild." [0,100]}
    \item{\code{thal }}{: numeric. "My thoughts wouldn't leave me alone today." [0,100]}
    \item{\code{thrc }}{: numeric. "My thoughts were racing today." [0,100]}
    \item{\code{thex }}{: numeric. "My thoughts were difficult to express today." [0,100]}
    \item{\code{otnl }}{: numeric. "Today I had the feeling that others did not like me." [0,100]}
    \item{\code{otrt }}{: numeric. "I felt that others could read my thoughts today." [0,100]}
    \item{\code{unrl }}{: numeric. "I felt unreal today." [0,100]}
    \item{\code{cont }}{: numeric. "I felt that other could control me today." [0,100]}
    \item{\code{enjy}}{: numeric. "I could experience pleasure when nice things happened today." [0,100]}
    \item{\code{undt }}{: numeric. "I felt like undertaking something today." [0,100]}

  }
}

\references{
Booij, S. H., Wichers, M., de Jonge, P., Sytema, S., van Os, J., Wunderink, L., and Wigman,J. T. W.  (2018).  Study protocol for a prospective cohort study examining the predictive potential  of  dynamic symptom networks for the onset and progression of  psychosis: the mapping individual routes of risk and resilience (mirorr) study.BMJ Open,8(1). Retrieved from https://bmjopen.bmj.com/content/8/1/e019059
}
\examples{
data(MIR)
str(MIR)
}
\keyword{datasets}
