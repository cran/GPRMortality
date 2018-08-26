\name{GPRMortalityChild}
\alias{GPRMortalityChild}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gaussian Process Regression for estimation and projection of child (under-5) mortality rates
}
\description{
 Searching for the latest methods of estimating mortality rates is a major concern. The need to have accurate and valid estimation of under-5 mortality rate led to apply more powerful and reliable methods. The main challenge is how to combine and integrate these different time series and how to produce unified estimates of mortality rates during a specified time span.

GPR is a Bayesian statistical model for estimating child (under-five age group) mortality which its data likelihood could be mortality rates from different data sources such as: Death Registration System, Censuses or surveys. This function produces estimations and uncertainty (95\% or any desirable percentiles) based on sampling and non-sampling errors due to variation in data sources.

The GP model utilizes Bayesian inference to update predicted mortality rates as a posterior in Bayes rule by combining data and a prior probability distribution over parameters in mean, covariance function, and the regression model. This package uses Markov Chain Monte Carlo (MCMC) to sample from posterior probability distribution by rstan package in R.

}
\usage{
GPRMortalityChild( data.mortality , data.mean ,minYear,maxYear,
                   nu,rho_,n.itr=4000,n.warm=3000,verbose=FALSE)
}
\arguments{
  \item{data.mortality}{
a data frame on mortality rates and population. It needs five compulsory basic variables: year, type of data sources, mortality rates and their corresponding population, location (countries and so on ).
}
  \item{data.mean}{
a data frame of data mean, with the same levels (year and location) as data.mortality.
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
  \item{n.itr}{
the number of iterations to use for running model in rstan. If not specified, then the default is 4000.
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

The algorithm for GPR package was developed in rstan.

The mortality rates will change to per 1000 population.
}
\value{
a list including two parts:


1.        a matrix of GPR result.

2.        a matrix of completeness bias information of Dearth Registration System by location including coefficients, Standard error of coefficient, and p-value of coefficients whether it can assess that DR system is biased in each location.

3.        a matrix of different data sources variances.
}
\references{
1.	Wang H, Liddell CA, Coates MM, et al. Global, regional, and national levels of neonatal, infant, and under-5 mortality during 1990-2013: a systematic analysis for the Global Burden of Disease Study 2013. Lancet 2014; published online May 2. http://dx.doi.org/10.1016/S0140-6736(14)60497-9.

2.	Wang H, Dwyer-Lindgren L, Lofgren KT, et al. Age-specific and sex-specific mortality in 187 countries, 1970-2010: a systematic analysis for the Global Burden of Disease Study 2010. Lancet 2012; 380: 2071-94.

3.	Mohammadi Y, Parsaeian M, Farzadfar F, et al. Levels and trends of child and adult mortality rates in the Islamic Republic of Iran, 1990-2013; protocol of the NASBOD study. Arch Iran Med 2014; 17: 176–81.

4. 	Parsaeian M, Farzadfar F, Zeraati H, et al. Application of spatiotemporal model to estimate burden of diseases, injuries and risk factors in Iran 1990-2013. Arch Iran Med 2014; 17: 28-33.

5.	Mehdipour P, Navidi I, Parsaeian M, Mohammadi Y, Moradi Lakeh M, Rezaei Darzi E, Nourijelyani K, Farzadfar F. . Application of Gaussian Process Regression (GPR) in estimating under-five mortality levels and trends in Iran 1990-2013, study protocol. Archives of Iranian medicine. 2014;17(3):189.

6.	Mohammadi Y, Parsaeian M, Mehdipour P, Khosravi A, Larijani B, Sheidaei A, et al. Measuring Iran's success in achieving Millennium Development Goal 4: a systematic analysis of under-5 mortality at national and subnational levels from 1990 to 2015. The Lancet Global Health. 2017;5(5):e537-e44.

7. Williams CK, Rasmussen CE. Gaussian processes for machine learning. the MIT Press. 2006;2(3):4.

}
\author{
Parinaz Mehdipour, Ali Ghanbari, Iman Navidi
}
\note{
Each data source must contin at least two data points in the time span between min and max years to be predicted.


min and max years determined the number of years in GPR results.


This package need Rtools to be installed. The version should match with R version.
}


\seealso{
\code{\link{GPRMortality}},\code{\link{GPRMortalitySummary}}
}
\examples{

 library("rstan")
 library("GPRMortality")
 head(data.mortality.child)
 head(data.mean.child)
 mortality <- data.mortality.child[data.mortality.child$location\%in\%c(0,5) ,]
 mean <- data.mean.child[data.mean.child$location\%in\%c(0,5) ,]
# WARNING: The following code will take a long time to run
\donttest{
 fit <- GPRMortalityChild(mortality,mean,minYear=1990,maxYear=2015,
     nu = 2,rho_ = 0.4 , n.itr=2000,n.warm=1000,verbose=TRUE)
 fit$simulation
 fit$variance
 fit_sum = GPRMortalitySummary(fit)
 fit_sum
}

}
