\name{data.mean.child}
\alias{data.mean.child}
\docType{data}
\title{
mean of child mortality rates
}
\description{
A data frame for GPR mean function and must contain an estimate for each year.

This data must contain year and location of each mortality rate. It can be different models for mortality rate estimations by year.
}
\usage{data("data.mean.child")}
\format{
  A data frame with 1728 observations on the following 3 variables.
  \describe{
    \item{\code{year}}{a numeric vector}
    \item{\code{location}}{a numeric vector}
    \item{\code{mean}}{a numeric vector}
  }
}
\details{
This data set is the result of Spatio-temporal Model that used in child mortality project in Iran from 1990 to 2015 for 31 provinces.
}
\source{
This data set was gathered in Non-Communicable Diseases Research Center affiliated to Endocrinology and Metabolism Research Institute,


Tehran University of Medical Sciences. ("http://www.ncdrc.info/")

the result of this data published in "http://dx.doi.org/10.1016/S2214-109X(17)30105-5" and is available on "https://data.mendeley.com/datasets/9z3pzd6rmd/1" .

}
\references{
1- Mohammadi Y, Parsaeian M, Farzadfar F, et al. Levels and trends of child and adult mortality rates in the Islamic Republic of Iran, 1990-2013; protocol of the NASBOD study. Arch Iran Med 2014; 17: 176–81.

2-Parsaeian M, Farzadfar F, Zeraati H, et al. Application of spatiotemporal model to estimate burden of diseases, injuries and risk factors in Iran 1990–2013. Arch Iran Med 2014; 17: 28–33.
}
\examples{
data(mean.rates.child)
## maybe str(mean.rates.child) ; plot(mean.rates.child) ...
}
\keyword{datasets}
