\name{birthwt}
\alias{birthwt}
\docType{data}
\title{
Birth weight data
}
\description{
Birth weight data collected at Baystate Medical Center, Springfield, MA during 1986.
}
\usage{data(birthwt)}
\format{
  A data frame with 189 observations on the following 10 variables.
  \describe{
    \item{\code{LOW}}{Indicator of birth weight less than 2.5 kg}
    \item{\code{AGE}}{Mother's age (in years)}
    \item{\code{LWT}}{Mother's weight (in pounds) at last menstrual period}
    \item{\code{RACE}}{Mother's race ('1' = white, '2' = black, '3' = other)}
    \item{\code{SMOKE}}{Smoking status during pregnancy}
    \item{\code{PTL}}{Number of previous premature labours}
    \item{\code{HT}}{History of hypertension}
    \item{\code{UI}}{Presence of uterine irritability}
    \item{\code{FTV}}{Number of physician visits during the first trimester}
    \item{\code{BWT}}{Recorded birth weight (in grams)}
  }
}

\source{
Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied Logistic Regression, third edition, Wiley, New York, Table 1.6.
}

\examples{
data(birthwt)
head(birthwt)
}
\keyword{datasets}
