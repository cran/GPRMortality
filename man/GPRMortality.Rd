\name{GPRMortality}
\alias{GPRMortality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gaussian Process Regression for estimation and projection each age-sex combination of mortality rates
}
\description{
Searching for the latest methods of estimating mortality rates is a major concern. The need to have accurate and valid estimation of age and sex mortality rate led to apply more powerful and reliable methods. The main challenge is how to combine and integrate these different time series and how to produce unified estimates of mortality rates during a specified time span.

GPR is a Bayesian statistical model for estimating adult (15-60 age group or age-specific) mortality which its data likelihood could be mortality rates from different data sources such as: Death Registration System, Censuses or surveys. This function produces a unique estimation and 95\% (or any desirable percentiles) uncertainty based on sampling and non-sampling errors due to variation in data sources.

The GP model utilizes Bayesian inference to update predicted mortality rates as a posterior in Bayes rule by combining data and a prior probability distribution over parameters in mean, covariance function, and the regression model. This package uses Markov Chain Monte Carlo (MCMC) to sample from posterior probability distribution by rstan package in R.


}
\usage{
GPRMortality(data.mortality,data.mean,minYear,maxYear,
             nu,rho_,product , n.itr=4000,n.warm=3000,verbose=FALSE)
}
\arguments{
  \item{data.mortality}{
a data frame on mortality rates, which needs seven compulsory basic variables: year, mortality rates, the completeness of corresponding mortality rates and their corresponding population; different levels of data such as age_cat, sex, location(countries and so on ) can be changed .
}
  \item{data.mean}{
a data frame of mean data, with the same levels including sex, age_cat, year and location.
}
  \item{minYear}{
min year will predicted
}
  \item{maxYear}{
max year will predicted
}
  \item{nu}{
the degree of differentiability parameter that is ranged from 0.2 to 2 and it controls the smoothness of samples driven from GP model.
}
  \item{rho_}{
a scale parameter that ranges from 0 to 1 and it defines the amount of correlation between years of data.
}
  \item{product}{
the value for calculating variance of Registration System, the default is 0.1 for high quality data.
}
  \item{n.itr}{
the number of iterations to use for running model in Rstan. If not specified, then the default is 4000.
}
  \item{n.warm}{
the number of samples for warm-up in MCMC. If not specified, the default is 3000.
}
 \item{verbose}{
TRUE or FALSE: flag indicating whether to print intermediate output from Stan on the console, which might be helpful for model debugging.
}
}
\details{
This package intra and extrapolates the mortality rates excluded from different data sources and various demographic models and produces uncertainty.

The algorithm  for GPR package was developed in rstan.

The mortality rates will change to per 100,000 population.

non-sampling variance is defined as 0.1 of completeness for countries/locations with high quality of registration more than 90\%.
}
\value{
	a matrix of GPR result.
}
\references{
1.	Wang H, Dwyer-Lindgren L, Lofgren KT, et al. Age-specific and sex-specific mortality in 187 countries, 1970-2010: a systematic analysis for the Global Burden of Disease Study 2010. Lancet 2012; 380: 2071–94.

2.	Rajaratnam JK, Marcus JR, Levin-Rector A, et al. Worldwide mortality in men and women aged 15–59 years from 1970 to 2010: a systematic analysis. Lancet 2010; published online April 30. DOI:10.1016/S0140-6736(10)60517-X.

3.	Mohammadi Y, Parsaeian M, Farzadfar F, et al. Levels and trends of child and adult mortality rates in the Islamic Republic of Iran, 1990-2013; protocol of the NASBOD study. Arch Iran Med 2014; 17: 176–81.

4.	Mehdipour P, Navidi I, Parsaeian M, Mohammadi Y, Moradi Lakeh M, Rezaei Darzi E, Nourijelyani K, Farzadfar F. . Application of Gaussian Process Regression (GPR) in estimating under-five mortality levels and trends in Iran 1990-2013, study protocol. Archives of Iranian medicine. 2014;17(3):189.

5. Williams CK, Rasmussen CE. Gaussian processes for machine learning. the MIT Press. 2006;2(3):4.
}
\author{
Parinaz Mehdipour, Ali Ghanbari, Iman Navidi
}
\note{
Each data source must contin at least two data points in the time span between min and max years to be predicted.

This package need Rtools to be installed. The version should match with R version.

min and max years determined the number of years in GPR results.
}


\seealso{
\code{\link{GPRMortalityChild}},\code{\link{GPRMortalitySummary}}
}
\examples{

 library("rstan")
 library("GPRMortality")
 head(data.mortality)
 head(data.mean)
 mortality <- data.mortality[data.mortality$location\%in\%c(0,5) &
           data.mortality$age_cat\%in\%c(1,10) & data.mortality$sex\%in\%c(0,1),]
 mean <- data.mean[data.mean$location\%in\%c(0,5) &
      data.mean$age_cat\%in\%c(1,10) & data.mean$sex\%in\%c(0,1),]
# WARNING: The following code will take a long time to run
\donttest{
 fit = GPRMortality(mortality,mean,minYear=1990,maxYear=2015,
                    nu = 2,rho_ =0.4 ,product =  0.1 ,verbose=TRUE)
 fit_sum = GPRMortalitySummary(fit)
 fit_sum
}

}
