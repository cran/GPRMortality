\name{GPRMortalitySummary}
\alias{GPRMortalitySummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary table for Gaussian Process Regression results
}
\description{
This function summarize the results by percentiles.

This function works for both child and age-sex specific mortality rates.
}
\usage{
GPRMortalitySummary(model,percentile=c(0.025,0.5,0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
a model object.
}
  \item{percentile}{
a vector of uncertainty interval percentile. The default is 0.025, 0.5 and 0.975.
}


}

\value{
a matrix of GPR result including percentile.
}

\author{
Parinaz Mehdipour, Ali Ghanbari, Iman Navidi
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


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
 fit = GPRMortality(mortality,mean,minYear =  1990,maxYear =  2015,
                    nu = 2,rho_ =0.4 ,product =  0.1 ,verbose=TRUE)
####### summary
 fit_sum = GPRMortalitySummary(fit)
 fit_sum
}
}
