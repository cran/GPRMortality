\name{data.mortality}
\alias{data.mortality}
\docType{data}
\title{
age-sex specific mortality rates data set
}
\description{
A data frame of mortality rates using as data likelihood of the GPR model.

This data set for child mortality rates should contain the year, sex, age_cat, location, completeness of each mortality rate and population of each raw data.
}
\usage{data("data.mortality")}
\format{
  A data frame with 395 observations on the following 7 variables.
  \describe{
    \item{\code{sex}}{a numeric vector}
    \item{\code{age_cat}}{a numeric vector}
    \item{\code{location}}{a numeric vector}
    \item{\code{year}}{a numeric vector}
    \item{\code{pop}}{a numeric vector}
    \item{\code{mortality}}{a numeric vector}
    \item{\code{completeness}}{a numeric vector}
  }
}
\details{
This is a real data gathered in Iran including 31 provinces and 26 years for age-sex specific and adult mortality.


}
\source{
This data set was gathered in Non-Communicable Diseases Research Center affiliated to Endocrinology and Metabolism Research Institute,


Tehran University of Medical Sciences. ("http://www.ncdrc.info/")
}
\references{
1-Mehdipour P, Navidi I, Parsaeian M, Mohammadi Y, Moradi Lakeh M, Rezaei Darzi E, Nourijelyani K, Farzadfar F. . Application of Gaussian Process Regression (GPR) in estimating under-five mortality levels and trends in Iran 1990-2013, study protocol. Archives of Iranian medicine. 2014;17(3):189.

2-Mohammadi Y, Parsaeian M, Mehdipour P, Khosravi A, Larijani B, Sheidaei A, et al. Measuring Iran's success in achieving Millennium Development Goal 4: a systematic analysis of under-5 mortality at national and subnational levels from 1990 to 2015. The Lancet Global Health. 2017;5(5):e537-e44.
}
\examples{
data(data.mortality)
## maybe str(data.mortality) ; plot(data.mortality) ...
}
\keyword{datasets}
