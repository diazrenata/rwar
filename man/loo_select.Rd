% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fxns.R
\name{loo_select}
\alias{loo_select}
\title{Select the simplest model within 1 se of the best model}
\usage{
loo_select(some_compares)
}
\arguments{
\item{some_compares}{from compare_both_stanarms}
}
\value{
df with model names for best model for each currency
}
\description{
Select the simplest model within 1 se of the best model
}
