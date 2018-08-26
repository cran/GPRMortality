\name{data.mean}
\alias{data.mean}
\docType{data}
\title{
mean of adult mortality rates
}
\description{
A data frame for GPR mean function and must contain an estimate for each year.

This data must contain year, sex, age_cat and location of each mortality rate. It can be different estimation models of mortality rates by year.
}
\usage{data("data.mean")}
\format{
  A data frame with 45942 observations on the following 5 variables.
  \describe{
    \item{\code{location}}{a numeric vector}
    \item{\code{age_cat}}{a numeric vector}
    \item{\code{sex}}{a numeric vector}
    \item{\code{year}}{a numeric vector}
    \item{\code{mean}}{a numeric vector}
  }
}
\details{
This data set is the result of Spatio-temporal Model that used in age-sex specific mortality rates as a part of NASBOD project in Iran from 1990 to 2015 for 31 provinces.
}
\source{
This data set was gathered in Non-Communicable Diseases Research Center affiliated to Endocrinology and Metabolism Research Institute,


Tehran University of Medical Sciences. ("http://www.ncdrc.info/")
}
\references{
1- Mohammadi Y, Parsaeian M, Farzadfar F, et al. Levels and trends of child and adult mortality rates in the Islamic Republic of Iran, 1990-2013; protocol of the NASBOD study. Arch Iran Med 2014; 17: 176–81.

2-Parsaeian M, Farzadfar F, Zeraati H, et al. Application of spatiotemporal model to estimate burden of diseases, injuries and risk factors in Iran 1990–2013. Arch Iran Med 2014; 17: 28–33.
}
\examples{
data(data.mean)
## maybe str(data.mean) ; plot(data.mean) ...
}
\keyword{datasets}
