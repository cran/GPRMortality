\name{data.mortality.child}
\alias{data.mortality.child}
\docType{data}
\title{
child mortality rates data
}
\description{
A data frame of mortality rates using as data likelihood of the GPR model.

This data set for child mortality rates should contain the year, location, name of each source of data, the population of each raw data and whether it comes from a Death Registration System or not.
}
\usage{
data("data.mortality.child")
}
\format{
  A data frame with 4107 observations on the following 6 variables.
  \describe{
    \item{\code{year}}{a numeric vector}
    \item{\code{location}}{a numeric vector}
    \item{\code{type}}{a character vector}
    \item{\code{mortality}}{a numeric vector}
    \item{\code{isDR}}{a numeric vector}
    \item{\code{pop}}{a numeric vector}
  }
}
\details{
This is a real data gathered in Iran including 31 provinces and 52 years for child mortality.

There are three types of data that include the summary (SBH) and complete birth history (CBH) data,

Census data and the Demographic and Health Survey (DHS) contain SBH questions only; the DHS contains CBH questions as well.

Subsequently, the census data of 1986, 1996, 2006, and 2011, the DHS data of 2000 and 2010 were determined as the sources of data.
}
\source{
This data set was gathered in Non-Communicable Diseases Research Center affiliated to Endocrinology and Metabolism Research Institute,


Tehran University of Medical Sciences. ("http://www.ncdrc.info/")

the result of this data published in "http://dx.doi.org/10.1016/S2214-109X(17)30105-5" and is available on "https://data.mendeley.com/datasets/9z3pzd6rmd/1" .
}
\references{
1-Mehdipour P, Navidi I, Parsaeian M, Mohammadi Y, Moradi Lakeh M, Rezaei Darzi E, Nourijelyani K, Farzadfar F. . Application of Gaussian Process Regression (GPR) in estimating under-five mortality levels and trends in Iran 1990-2013, study protocol. Archives of Iranian medicine. 2014;17(3):189.

2-Mohammadi Y, Parsaeian M, Mehdipour P, Khosravi A, Larijani B, Sheidaei A, et al. Measuring Iran's success in achieving Millennium Development Goal 4: a systematic analysis of under-5 mortality at national and subnational levels from 1990 to 2015. The Lancet Global Health. 2017;5(5):e537-e44.
}
\examples{
data(mortality.rate.child)
## maybe str(mortality.rate.child) ; plot(mortality.rate.child) ...
}
\keyword{datasets}
